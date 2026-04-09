#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.config) as handle:
        config = yaml.safe_load(handle)

    metadata = {
        "reference_fasta": config["reference"]["fasta"],
        "reference_gtf": config["reference"]["gtf"],
        "genome_split_dir": config["reference"]["split_dir"],
        "bowtie_index_prefix": config["reference"]["bowtie_index_prefix"],
        "exon_k": config["params"]["exon_k"],
        "utr_k": config["params"]["utr_k"],
        "mismatch": config["params"]["mismatch"],
        "exon_kmer_mappability": config["mappability"]["exon_bedgraph"],
        "utr_kmer_mappability": config["mappability"]["utr_bedgraph"],
        "output_root": config["output_root"],
    }
    Path(args.output).write_text(json.dumps(metadata, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
