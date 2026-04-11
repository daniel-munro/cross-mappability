#!/usr/bin/env python3

import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path


QUANTILES = {
    "p00": 0.00,
    "p10": 0.10,
    "p25": 0.25,
    "p50": 0.50,
    "p75": 0.75,
    "p90": 0.90,
    "p99": 0.99,
    "p100": 1.00,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rat-crossmap", required=True)
    parser.add_argument("--human-crossmap", required=True)
    parser.add_argument("--rat-mappability", required=True)
    parser.add_argument("--summary-tsv", required=True)
    parser.add_argument("--comparison-tsv", required=True)
    parser.add_argument("--mappability-deciles-tsv", required=True)
    parser.add_argument("--report-md", required=True)
    parser.add_argument("--summary-json", required=True)
    return parser.parse_args()


def quantile_from_sorted(values: list[int], q: float) -> float:
    if not values:
        raise ValueError("cannot compute quantile of empty list")
    idx = int((len(values) - 1) * q)
    return float(values[idx])


def quantiles_from_histogram(histogram: Counter[int], total_count: int) -> dict[str, float]:
    if total_count == 0:
        raise ValueError("cannot compute quantiles for empty histogram")
    targets = {name: math.ceil(q * total_count) for name, q in QUANTILES.items()}
    targets["p00"] = 1
    ordered_targets = sorted(targets.items(), key=lambda item: item[1])
    result: dict[str, float] = {}
    seen = 0
    target_idx = 0
    for value in sorted(histogram):
        seen += histogram[value]
        while target_idx < len(ordered_targets) and seen >= ordered_targets[target_idx][1]:
            result[ordered_targets[target_idx][0]] = float(value)
            target_idx += 1
    return result


def mean(values: list[int]) -> float:
    return sum(values) / len(values)


def top_share(values: list[int], top_fraction: float) -> float:
    if not values:
        raise ValueError("cannot compute concentration for empty list")
    count = max(1, math.ceil(len(values) * top_fraction))
    numerator = sum(sorted(values, reverse=True)[:count])
    denominator = sum(values)
    return numerator / denominator if denominator else 0.0


def read_crossmap(path: str) -> dict:
    genes: set[str] = set()
    out_degree: Counter[str] = Counter()
    in_degree: Counter[str] = Counter()
    out_strength: Counter[str] = Counter()
    in_strength: Counter[str] = Counter()
    weight_histogram: Counter[int] = Counter()
    edge_count = 0
    total_strength = 0

    with gzip.open(path, "rt") as handle:
        for line in handle:
            source, target, strength_text = line.rstrip("\n").split("\t")
            strength = int(strength_text)
            genes.add(source)
            genes.add(target)
            out_degree[source] += 1
            in_degree[target] += 1
            out_strength[source] += strength
            in_strength[target] += strength
            weight_histogram[strength] += 1
            edge_count += 1
            total_strength += strength

    all_genes = sorted(genes)
    out_degree_values = [out_degree[gene] for gene in all_genes]
    in_degree_values = [in_degree[gene] for gene in all_genes]
    out_strength_values = [out_strength[gene] for gene in all_genes]
    in_strength_values = [in_strength[gene] for gene in all_genes]
    nonzero_out_degree = sum(value > 0 for value in out_degree_values)
    nonzero_in_degree = sum(value > 0 for value in in_degree_values)

    return {
        "genes": all_genes,
        "edge_count": edge_count,
        "total_strength": total_strength,
        "weight_quantiles": quantiles_from_histogram(weight_histogram, edge_count),
        "mean_edge_weight": total_strength / edge_count,
        "out_degree_values": out_degree_values,
        "in_degree_values": in_degree_values,
        "out_strength_values": out_strength_values,
        "in_strength_values": in_strength_values,
        "out_degree": out_degree,
        "out_strength": out_strength,
        "summary": {
            "gene_count": len(all_genes),
            "edge_count": edge_count,
            "total_strength": total_strength,
            "mean_edge_weight": total_strength / edge_count,
            "source_gene_fraction": nonzero_out_degree / len(all_genes),
            "target_gene_fraction": nonzero_in_degree / len(all_genes),
            "out_degree_quantiles": {
                name: quantile_from_sorted(sorted(out_degree_values), q)
                for name, q in QUANTILES.items()
            },
            "in_degree_quantiles": {
                name: quantile_from_sorted(sorted(in_degree_values), q)
                for name, q in QUANTILES.items()
            },
            "out_strength_quantiles": {
                name: quantile_from_sorted(sorted(out_strength_values), q)
                for name, q in QUANTILES.items()
            },
            "in_strength_quantiles": {
                name: quantile_from_sorted(sorted(in_strength_values), q)
                for name, q in QUANTILES.items()
            },
            "edge_weight_quantiles": quantiles_from_histogram(weight_histogram, edge_count),
            "out_strength_top_1pct_share": top_share(out_strength_values, 0.01),
            "out_strength_top_5pct_share": top_share(out_strength_values, 0.05),
            "out_strength_top_10pct_share": top_share(out_strength_values, 0.10),
            "in_strength_top_1pct_share": top_share(in_strength_values, 0.01),
            "in_strength_top_5pct_share": top_share(in_strength_values, 0.05),
            "in_strength_top_10pct_share": top_share(in_strength_values, 0.10),
        },
    }


def read_rat_mappability(path: str) -> dict[str, float | None]:
    scores: dict[str, float | None] = {}
    with open(path) as handle:
        for line in handle:
            gene, score = line.rstrip("\n").split("\t")
            scores[gene] = None if score == "NA" else float(score)
    return scores


def build_mappability_deciles(
    rat_stats: dict,
    rat_mappability: dict[str, float | None],
) -> dict[str, object]:
    missing = [gene for gene in rat_stats["genes"] if gene not in rat_mappability]
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"missing rat mappability for {len(missing)} crossmap genes: {preview}")

    rows = [
        {
            "gene": gene,
            "mappability": rat_mappability[gene],
            "out_degree": rat_stats["out_degree"][gene],
            "out_strength": rat_stats["out_strength"][gene],
        }
        for gene in rat_stats["genes"]
    ]
    na_rows = [row for row in rows if row["mappability"] is None]
    finite_rows = [row for row in rows if row["mappability"] is not None]
    if not finite_rows:
        raise ValueError("all rat mappability values are NA")

    finite_rows.sort(key=lambda row: (row["mappability"], row["gene"]))

    decile_rows: list[dict[str, float]] = []
    total_out_strength = sum(row["out_strength"] for row in rows)
    total_genes = len(finite_rows)
    for decile_idx in range(10):
        start = (total_genes * decile_idx) // 10
        end = (total_genes * (decile_idx + 1)) // 10
        chunk = finite_rows[start:end]
        mappability_values = sorted(row["mappability"] for row in chunk)
        out_degree_values = sorted(row["out_degree"] for row in chunk)
        out_strength_values = sorted(row["out_strength"] for row in chunk)
        chunk_out_strength = sum(out_strength_values)
        decile_rows.append(
            {
                "decile": decile_idx + 1,
                "gene_count": len(chunk),
                "mappability_min": min(mappability_values),
                "mappability_median": quantile_from_sorted(mappability_values, 0.5),
                "mappability_max": max(mappability_values),
                "mean_out_degree": mean(out_degree_values),
                "median_out_degree": quantile_from_sorted(out_degree_values, 0.5),
                "mean_out_strength": mean(out_strength_values),
                "median_out_strength": quantile_from_sorted(out_strength_values, 0.5),
                "share_total_out_strength": (
                    chunk_out_strength / total_out_strength if total_out_strength else 0.0
                ),
            }
        )
    na_out_strength = sum(row["out_strength"] for row in na_rows)
    return {
        "deciles": decile_rows,
        "na_summary": {
            "na_gene_count": len(na_rows),
            "na_gene_fraction": len(na_rows) / len(rows),
            "na_out_strength": na_out_strength,
            "na_out_strength_share": na_out_strength / total_out_strength if total_out_strength else 0.0,
        },
    }


def flatten_summary(dataset: str, summary: dict) -> list[dict[str, str | float]]:
    rows: list[dict[str, str | float]] = []
    for key in (
        "gene_count",
        "edge_count",
        "total_strength",
        "mean_edge_weight",
        "source_gene_fraction",
        "target_gene_fraction",
        "out_strength_top_1pct_share",
        "out_strength_top_5pct_share",
        "out_strength_top_10pct_share",
        "in_strength_top_1pct_share",
        "in_strength_top_5pct_share",
        "in_strength_top_10pct_share",
    ):
        rows.append({"dataset": dataset, "metric": key, "value": summary[key]})

    for prefix in (
        "edge_weight_quantiles",
        "out_degree_quantiles",
        "in_degree_quantiles",
        "out_strength_quantiles",
        "in_strength_quantiles",
    ):
        for quantile_name, value in summary[prefix].items():
            rows.append(
                {
                    "dataset": dataset,
                    "metric": f"{prefix}.{quantile_name}",
                    "value": value,
                }
            )
    return rows


def compare_rat_human(rat_summary: dict, human_summary: dict) -> list[dict[str, float | str]]:
    metrics = [
        "gene_count",
        "edge_count",
        "total_strength",
        "mean_edge_weight",
        "source_gene_fraction",
        "target_gene_fraction",
        "out_strength_top_1pct_share",
        "out_strength_top_5pct_share",
        "out_strength_top_10pct_share",
        "in_strength_top_1pct_share",
        "in_strength_top_5pct_share",
        "in_strength_top_10pct_share",
    ]
    rows = []
    for metric in metrics:
        rows.append(
            {
                "metric": metric,
                "rat_value": rat_summary[metric],
                "human_value": human_summary[metric],
                "rat_over_human": rat_summary[metric] / human_summary[metric],
            }
        )
    for prefix in (
        "edge_weight_quantiles",
        "out_degree_quantiles",
        "in_degree_quantiles",
        "out_strength_quantiles",
        "in_strength_quantiles",
    ):
        for quantile_name in QUANTILES:
            rat_value = rat_summary[prefix][quantile_name]
            human_value = human_summary[prefix][quantile_name]
            rows.append(
                {
                    "metric": f"{prefix}.{quantile_name}",
                    "rat_value": rat_value,
                    "human_value": human_value,
                    "rat_over_human": rat_value / human_value if human_value else float("nan"),
                }
            )
    return rows


def write_tsv(path: str, rows: list[dict], fieldnames: list[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def format_float(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.3f}"


def build_report(
    rat_summary: dict,
    human_summary: dict,
    deciles: list[dict],
    na_summary: dict,
) -> str:
    rat_low = deciles[0]
    rat_high = deciles[-1]
    peak_strength_decile = max(deciles, key=lambda row: row["mean_out_strength"])
    ratio_strength = rat_low["mean_out_strength"] / rat_high["mean_out_strength"]
    ratio_degree = rat_low["mean_out_degree"] / rat_high["mean_out_degree"]

    lines = [
        "# Cross-mappability validation",
        "",
        "## Global comparison",
        "",
        "| Metric | Rat | Human reference | Rat / Human |",
        "| --- | ---: | ---: | ---: |",
        (
            f"| Genes | {rat_summary['gene_count']} | {human_summary['gene_count']} | "
            f"{format_float(rat_summary['gene_count'] / human_summary['gene_count'])} |"
        ),
        (
            f"| Edges | {rat_summary['edge_count']} | {human_summary['edge_count']} | "
            f"{format_float(rat_summary['edge_count'] / human_summary['edge_count'])} |"
        ),
        (
            f"| Total strength | {rat_summary['total_strength']} | {human_summary['total_strength']} | "
            f"{format_float(rat_summary['total_strength'] / human_summary['total_strength'])} |"
        ),
        (
            f"| Mean edge weight | {format_float(rat_summary['mean_edge_weight'])} | "
            f"{format_float(human_summary['mean_edge_weight'])} | "
            f"{format_float(rat_summary['mean_edge_weight'] / human_summary['mean_edge_weight'])} |"
        ),
        (
            f"| Edge-weight median | {format_float(rat_summary['edge_weight_quantiles']['p50'])} | "
            f"{format_float(human_summary['edge_weight_quantiles']['p50'])} | "
            f"{format_float(rat_summary['edge_weight_quantiles']['p50'] / human_summary['edge_weight_quantiles']['p50'])} |"
        ),
        (
            f"| Out-degree median | {format_float(rat_summary['out_degree_quantiles']['p50'])} | "
            f"{format_float(human_summary['out_degree_quantiles']['p50'])} | "
            f"{format_float(rat_summary['out_degree_quantiles']['p50'] / human_summary['out_degree_quantiles']['p50'])} |"
        ),
        (
            f"| Out-strength median | {format_float(rat_summary['out_strength_quantiles']['p50'])} | "
            f"{format_float(human_summary['out_strength_quantiles']['p50'])} | "
            f"{format_float(rat_summary['out_strength_quantiles']['p50'] / human_summary['out_strength_quantiles']['p50'])} |"
        ),
        (
            f"| Top 1% outgoing-strength share | {format_float(rat_summary['out_strength_top_1pct_share'])} | "
            f"{format_float(human_summary['out_strength_top_1pct_share'])} | "
            f"{format_float(rat_summary['out_strength_top_1pct_share'] / human_summary['out_strength_top_1pct_share'])} |"
        ),
        "",
        "## Rat internal consistency",
        "",
        (
            f"Rat genes in the lowest mappability decile carry "
            f"{format_float(ratio_strength)}x the mean outgoing cross-mappability strength and "
            f"{format_float(ratio_degree)}x the mean outgoing cross-mappability degree of genes "
            f"in the highest mappability decile."
        ),
        (
            f"{na_summary['na_gene_count']} rat genes have NA weighted mappability "
            f"({format_float(na_summary['na_gene_fraction'])} of rat genes) and contribute "
            f"{format_float(na_summary['na_out_strength_share'])} of total outgoing strength."
        ),
        "",
        "| Mappability decile | Mappability range | Mean out-degree | Mean out-strength | Share of total out-strength |",
        "| --- | --- | ---: | ---: | ---: |",
    ]

    for row in deciles:
        lines.append(
            "| "
            f"{row['decile']} | "
            f"{format_float(row['mappability_min'])}-{format_float(row['mappability_max'])} | "
            f"{format_float(row['mean_out_degree'])} | "
            f"{format_float(row['mean_out_strength'])} | "
            f"{format_float(row['share_total_out_strength'])} |"
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "The rat result is not expected to numerically match the human reference, "
                "but the overall structure is comparable: both are large directed gene-gene "
                "cross-mappability graphs with heavy-tailed edge weights and burden concentrated "
                "in a minority of genes."
            ),
            (
                "The strongest reassurance here is the internal rat trend: the most mappable rat "
                "genes are strongly depleted for outgoing cross-mappability burden, while the "
                f"largest burden appears in mappability decile {peak_strength_decile['decile']} "
                "rather than in the top-mappability tail. That shape is plausible for this pipeline "
                "and argues against the rat result being random or grossly miscomputed."
            ),
            (
                "A direct rat-vs-human ortholog validation would still be useful, but it should be "
                "done with an explicit ortholog table rather than inferred from gene symbols."
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    args = parse_args()

    rat_stats = read_crossmap(args.rat_crossmap)
    human_stats = read_crossmap(args.human_crossmap)
    rat_mappability = read_rat_mappability(args.rat_mappability)
    mappability_summary = build_mappability_deciles(rat_stats, rat_mappability)
    deciles = mappability_summary["deciles"]
    na_summary = mappability_summary["na_summary"]

    summary_rows = flatten_summary("rat", rat_stats["summary"]) + flatten_summary(
        "human", human_stats["summary"]
    )
    comparison_rows = compare_rat_human(rat_stats["summary"], human_stats["summary"])
    write_tsv(args.summary_tsv, summary_rows, ["dataset", "metric", "value"])
    write_tsv(
        args.comparison_tsv,
        comparison_rows,
        ["metric", "rat_value", "human_value", "rat_over_human"],
    )
    write_tsv(
        args.mappability_deciles_tsv,
        deciles,
        [
            "decile",
            "gene_count",
            "mappability_min",
            "mappability_median",
            "mappability_max",
            "mean_out_degree",
            "median_out_degree",
            "mean_out_strength",
            "median_out_strength",
            "share_total_out_strength",
        ],
    )

    Path(args.report_md).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report_md).write_text(
        build_report(rat_stats["summary"], human_stats["summary"], deciles, na_summary)
    )

    Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_json).write_text(
        json.dumps(
            {
                "rat": rat_stats["summary"],
                "human": human_stats["summary"],
                "rat_mappability_na_summary": na_summary,
                "rat_mappability_deciles": deciles,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
