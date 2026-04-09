#!/usr/bin/env python3

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass
class TranscriptRecord:
    gene_id: str
    transcript_id: str
    gene_name: str
    gene_type: str
    transcript_type: str
    chrom: str
    strand: str
    source: str
    exons: list
    cds: list


def parse_attributes(attr_field: str) -> dict[str, list[str] | str]:
    attributes: dict[str, list[str] | str] = {}
    for raw_item in attr_field.strip().split(";"):
        item = raw_item.strip()
        if not item:
            continue
        if " " not in item:
            key = item
            value = ""
        else:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        if key == "tag":
            attributes.setdefault("tag", []).append(value)
        else:
            attributes[key] = value
    if "gene_name" not in attributes:
        attributes["gene_name"] = attributes.get("gene", attributes.get("gene_id", ""))
    if "gene_type" not in attributes:
        attributes["gene_type"] = attributes.get("gene_biotype", "unknown")
    if "transcript_type" not in attributes:
        attributes["transcript_type"] = attributes.get(
            "transcript_biotype", attributes["gene_type"]
        )
    return attributes


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def subtract_intervals(
    intervals: list[tuple[int, int]], subtractors: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    if not subtractors:
        return intervals[:]
    result: list[tuple[int, int]] = []
    for start, end in intervals:
        segments = [(start, end)]
        for sub_start, sub_end in subtractors:
            next_segments = []
            for seg_start, seg_end in segments:
                if sub_end < seg_start or sub_start > seg_end:
                    next_segments.append((seg_start, seg_end))
                    continue
                if sub_start > seg_start:
                    next_segments.append((seg_start, sub_start - 1))
                if sub_end < seg_end:
                    next_segments.append((sub_end + 1, seg_end))
            segments = next_segments
            if not segments:
                break
        result.extend(segments)
    return result


def load_gene_filter(path: str | None) -> set[str] | None:
    if path is None:
        return None
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def load_transcripts(gtf_path: str, gene_filter: set[str] | None) -> dict[str, TranscriptRecord]:
    transcripts: dict[str, TranscriptRecord] = {}

    with open(gtf_path) as handle:
        for raw_line in handle:
            if raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            chrom, source, feature, start, end, _, strand, _, attrs = fields
            attributes = parse_attributes(attrs)
            gene_id = attributes.get("gene_id", "")
            if gene_filter is not None and gene_id not in gene_filter:
                continue

            if feature not in {"transcript", "exon", "CDS"}:
                continue

            transcript_id = str(attributes.get("transcript_id", "")).strip()
            if not transcript_id:
                continue

            if transcript_id not in transcripts:
                transcripts[transcript_id] = TranscriptRecord(
                    gene_id=gene_id,
                    transcript_id=transcript_id,
                    gene_name=str(attributes.get("gene_name", gene_id)),
                    gene_type=str(attributes.get("gene_type", "unknown")),
                    transcript_type=str(
                        attributes.get(
                            "transcript_type", attributes.get("gene_type", "unknown")
                        )
                    ),
                    chrom=chrom,
                    strand=strand,
                    source=source,
                    exons=[],
                    cds=[],
                )

            record = transcripts[transcript_id]
            interval = (int(start), int(end))
            if feature == "exon":
                record.exons.append(interval)
            elif feature == "CDS":
                record.cds.append(interval)

    return transcripts


def transcript_feature_rows(
    transcripts: dict[str, TranscriptRecord],
) -> tuple[list[dict], dict]:
    rows: list[dict] = []
    summary = {
        "transcripts": 0,
        "genes_with_exons": set(),
        "genes_with_utrs": set(),
        "exon_rows": 0,
        "utr_rows": 0,
    }

    for transcript in transcripts.values():
        if not transcript.exons:
            continue
        summary["transcripts"] += 1
        merged_exons = merge_intervals(transcript.exons)
        merged_cds = merge_intervals(transcript.cds)
        utr_intervals = subtract_intervals(merged_exons, merged_cds) if merged_cds else []

        summary["genes_with_exons"].add(transcript.gene_id)
        base_attrs = {
            "gene_id": transcript.gene_id,
            "transcript_id": transcript.transcript_id,
            "gene_name": transcript.gene_name,
            "gene_type": transcript.gene_type,
            "transcript_type": transcript.transcript_type,
        }
        for start, end in merged_exons:
            rows.append(
                {
                    "chrom": transcript.chrom,
                    "source": transcript.source,
                    "feature": "exon",
                    "start": start,
                    "end": end,
                    "strand": transcript.strand,
                    "attributes": base_attrs.copy(),
                }
            )
            summary["exon_rows"] += 1

        for start, end in utr_intervals:
            rows.append(
                {
                    "chrom": transcript.chrom,
                    "source": transcript.source,
                    "feature": "UTR",
                    "start": start,
                    "end": end,
                    "strand": transcript.strand,
                    "attributes": base_attrs.copy(),
                }
            )
            summary["utr_rows"] += 1
            summary["genes_with_utrs"].add(transcript.gene_id)

    summary["genes_with_exons"] = len(summary["genes_with_exons"])
    summary["genes_with_utrs"] = len(summary["genes_with_utrs"])
    return rows, summary


def format_attributes(attributes: dict[str, str]) -> str:
    ordered_keys = [
        "gene_id",
        "transcript_id",
        "gene_name",
        "gene_type",
        "transcript_type",
    ]
    return " ".join(f'{key} "{attributes[key]}";' for key in ordered_keys if key in attributes)


def write_gtf(rows: list[dict], output_path: str) -> None:
    with open(output_path, "w") as handle:
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        row["chrom"],
                        row["source"],
                        row["feature"],
                        str(row["start"]),
                        str(row["end"]),
                        ".",
                        row["strand"],
                        ".",
                        format_attributes(row["attributes"]),
                    ]
                )
                + "\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out-crossmap-gtf", required=True)
    parser.add_argument("--gene-filter")
    parser.add_argument("--summary-json")
    args = parser.parse_args()

    gene_filter = load_gene_filter(args.gene_filter)
    transcripts = load_transcripts(args.gtf, gene_filter)
    rows, summary = transcript_feature_rows(transcripts)
    write_gtf(rows, args.out_crossmap_gtf)

    if args.summary_json:
        Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
