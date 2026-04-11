#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rat-crossmap", required=True)
    parser.add_argument("--human-crossmap", required=True)
    parser.add_argument("--orthologs", required=True)
    parser.add_argument("--joined-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--report-md", required=True)
    return parser.parse_args()


def summarize_crossmap(path: str, strip_version: bool) -> dict[str, dict[str, int]]:
    out_degree: Counter[str] = Counter()
    in_degree: Counter[str] = Counter()
    out_strength: Counter[str] = Counter()
    in_strength: Counter[str] = Counter()

    with gzip.open(path, "rt") as handle:
        for line in handle:
            gene_a, gene_b, score_text = line.rstrip("\n").split("\t")
            if strip_version:
                gene_a = gene_a.split(".", 1)[0]
                gene_b = gene_b.split(".", 1)[0]
            score = int(score_text)
            out_degree[gene_a] += 1
            in_degree[gene_b] += 1
            out_strength[gene_a] += score
            in_strength[gene_b] += score

    genes = sorted(set(out_degree) | set(in_degree))
    return {
        gene: {
            "out_degree": out_degree[gene],
            "in_degree": in_degree[gene],
            "out_strength": out_strength[gene],
            "in_strength": in_strength[gene],
            "total_degree": out_degree[gene] + in_degree[gene],
            "total_strength": out_strength[gene] + in_strength[gene],
        }
        for gene in genes
    }


def rankdata(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda idx: values[idx])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        avg_rank = (i + j - 1) / 2 + 1
        for k in range(i, j):
            ranks[order[k]] = avg_rank
        i = j
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    num = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    den_x = math.sqrt(sum((a - mean_x) ** 2 for a in x))
    den_y = math.sqrt(sum((b - mean_y) ** 2 for b in y))
    return num / (den_x * den_y) if den_x and den_y else float("nan")


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(rankdata(x), rankdata(y))


def read_orthologs(path: str) -> list[dict[str, str]]:
    with gzip.open(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def one_to_one_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows = [row for row in rows if row["orthology_type"] == "ortholog_one2one" and row["rat_gene_symbol"]]
    human_counts = Counter(row["human_ensembl_gene_id"] for row in rows)
    rat_counts = Counter(row["rat_gene_symbol"] for row in rows)
    return [
        row
        for row in rows
        if human_counts[row["human_ensembl_gene_id"]] == 1 and rat_counts[row["rat_gene_symbol"]] == 1
    ]


def write_joined(path: str, rows: list[dict]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    return f"{value:.3f}"


def main() -> None:
    args = parse_args()
    rat = summarize_crossmap(args.rat_crossmap, strip_version=False)
    human = summarize_crossmap(args.human_crossmap, strip_version=True)
    ortholog_rows = one_to_one_rows(read_orthologs(args.orthologs))

    joined_rows = []
    for row in ortholog_rows:
        human_id = row["human_ensembl_gene_id"]
        rat_symbol = row["rat_gene_symbol"]
        if human_id not in human or rat_symbol not in rat:
            continue
        human_stats = human[human_id]
        rat_stats = rat[rat_symbol]
        joined_rows.append(
            {
                "human_ensembl_gene_id": human_id,
                "rat_gene_symbol": rat_symbol,
                "rat_ensembl_gene_id": row["rat_ensembl_gene_id"],
                "human_out_degree": human_stats["out_degree"],
                "rat_out_degree": rat_stats["out_degree"],
                "human_in_degree": human_stats["in_degree"],
                "rat_in_degree": rat_stats["in_degree"],
                "human_out_strength": human_stats["out_strength"],
                "rat_out_strength": rat_stats["out_strength"],
                "human_in_strength": human_stats["in_strength"],
                "rat_in_strength": rat_stats["in_strength"],
                "human_total_strength": human_stats["total_strength"],
                "rat_total_strength": rat_stats["total_strength"],
            }
        )

    if not joined_rows:
        raise ValueError("no ortholog pairs overlap the human and rat crossmap outputs")

    write_joined(args.joined_tsv, joined_rows)

    human_out_strength = [math.log1p(row["human_out_strength"]) for row in joined_rows]
    rat_out_strength = [math.log1p(row["rat_out_strength"]) for row in joined_rows]
    human_in_strength = [math.log1p(row["human_in_strength"]) for row in joined_rows]
    rat_in_strength = [math.log1p(row["rat_in_strength"]) for row in joined_rows]
    human_total_strength = [math.log1p(row["human_total_strength"]) for row in joined_rows]
    rat_total_strength = [math.log1p(row["rat_total_strength"]) for row in joined_rows]

    summary = {
        "ortholog_pairs_total": len(ortholog_rows),
        "ortholog_pairs_compared": len(joined_rows),
        "ortholog_pairs_overlap_fraction": len(joined_rows) / len(ortholog_rows),
        "log1p_out_strength_pearson": pearson(human_out_strength, rat_out_strength),
        "log1p_out_strength_spearman": spearman(human_out_strength, rat_out_strength),
        "log1p_in_strength_pearson": pearson(human_in_strength, rat_in_strength),
        "log1p_in_strength_spearman": spearman(human_in_strength, rat_in_strength),
        "log1p_total_strength_pearson": pearson(human_total_strength, rat_total_strength),
        "log1p_total_strength_spearman": spearman(human_total_strength, rat_total_strength),
    }

    summary_path = Path(args.summary_json)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    top_rows = sorted(
        joined_rows,
        key=lambda row: row["human_total_strength"] + row["rat_total_strength"],
        reverse=True,
    )[:10]

    report_lines = [
        "# Ortholog comparison",
        "",
        "## Summary",
        "",
        f"- One-to-one human-rat ortholog pairs retrieved: {summary['ortholog_pairs_total']}",
        f"- One-to-one ortholog pairs present in both crossmap outputs: {summary['ortholog_pairs_compared']}",
        f"- Overlap fraction among retrieved one-to-one orthologs: {fmt(summary['ortholog_pairs_overlap_fraction'])}",
        f"- log1p(outgoing strength) Pearson: {fmt(summary['log1p_out_strength_pearson'])}",
        f"- log1p(outgoing strength) Spearman: {fmt(summary['log1p_out_strength_spearman'])}",
        f"- log1p(incoming strength) Pearson: {fmt(summary['log1p_in_strength_pearson'])}",
        f"- log1p(incoming strength) Spearman: {fmt(summary['log1p_in_strength_spearman'])}",
        f"- log1p(total strength) Pearson: {fmt(summary['log1p_total_strength_pearson'])}",
        f"- log1p(total strength) Spearman: {fmt(summary['log1p_total_strength_spearman'])}",
        "",
        "## Interpretation",
        "",
        "The ortholog-level correspondence is positive but weak. That is enough to argue against a completely species-random result, but it is not strong evidence that rat and human cross-mappability burdens are conserved gene by gene.",
        "The main use of this comparison is reassurance at the level of broad tendency and overlap, not precise cross-species prediction for individual genes.",
        "",
        "## Top ortholog pairs by combined cross-mappability burden",
        "",
        "| Human Ensembl | Rat symbol | Human total strength | Rat total strength |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in top_rows:
        report_lines.append(
            f"| {row['human_ensembl_gene_id']} | {row['rat_gene_symbol']} | "
            f"{row['human_total_strength']} | {row['rat_total_strength']} |"
        )
    report_lines.append("")

    Path(args.report_md).write_text("\n".join(report_lines))


if __name__ == "__main__":
    main()
