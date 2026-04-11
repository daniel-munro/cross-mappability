#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests


SERVER = "https://rest.ensembl.org"
HOMOLOGY_TEMPLATE = (
    SERVER
    + "/homology/id/human/{gene_id}?target_species=rattus_norvegicus;"
    + "type=orthologues;sequence=none;format=condensed"
)
LOOKUP_URL = SERVER + "/lookup/id"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--human-crossmap", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args()


def read_human_gene_ids(path: str) -> list[str]:
    gene_ids: set[str] = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            gene_a, gene_b, _ = line.rstrip("\n").split("\t")
            gene_ids.add(gene_a.split(".", 1)[0])
            gene_ids.add(gene_b.split(".", 1)[0])
    return sorted(gene_ids)


def fetch_homology(session: requests.Session, gene_id: str, cache_dir: Path) -> list[dict]:
    cache_path = cache_dir / f"{gene_id}.json"
    if cache_path.exists():
        payload = json.loads(cache_path.read_text())
    else:
        url = HOMOLOGY_TEMPLATE.format(gene_id=gene_id)
        last_error = None
        for attempt in range(5):
            try:
                response = session.get(
                    url,
                    headers={"Content-Type": "application/json", "Accept": "application/json"},
                    timeout=60,
                )
                response.raise_for_status()
                payload = response.json()
                cache_path.write_text(json.dumps(payload))
                break
            except Exception as exc:  # pragma: no cover - network errors are environment-dependent
                last_error = exc
                time.sleep(1 + attempt)
        else:
            raise RuntimeError(f"failed to fetch homology for {gene_id}") from last_error

    data = payload.get("data", [])
    if not data:
        return []
    homologies = data[0].get("homologies", [])
    return [
        {
            "human_ensembl_gene_id": gene_id,
            "rat_ensembl_gene_id": homology["id"],
            "orthology_type": homology["type"],
        }
        for homology in homologies
        if homology.get("species") == "rattus_norvegicus"
    ]


def lookup_display_names(session: requests.Session, ids: list[str]) -> dict[str, str]:
    display_names: dict[str, str] = {}
    for start in range(0, len(ids), 1000):
        chunk = ids[start : start + 1000]
        response = session.post(
            LOOKUP_URL,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            json={"ids": chunk},
            timeout=120,
        )
        response.raise_for_status()
        payload = response.json()
        for gene_id, record in payload.items():
            if record is not None:
                display_names[gene_id] = record.get("display_name")
    return display_names


def main() -> None:
    args = parse_args()
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    human_gene_ids = read_human_gene_ids(args.human_crossmap)

    session = requests.Session()
    rows: list[dict] = []
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(fetch_homology, session, gene_id, cache_dir): gene_id
            for gene_id in human_gene_ids
        }
        for future in as_completed(futures):
            rows.extend(future.result())

    rat_ids = sorted({row["rat_ensembl_gene_id"] for row in rows})
    rat_symbols = lookup_display_names(session, rat_ids)
    for row in rows:
        row["rat_gene_symbol"] = rat_symbols.get(row["rat_ensembl_gene_id"])

    rows.sort(
        key=lambda row: (
            row["human_ensembl_gene_id"],
            row["orthology_type"],
            row["rat_ensembl_gene_id"],
        )
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "human_ensembl_gene_id",
                "rat_ensembl_gene_id",
                "rat_gene_symbol",
                "orthology_type",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
