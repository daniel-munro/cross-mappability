#!/usr/bin/env python3

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kmer-dir", required=True)
    parser.add_argument("--stamp", required=True)
    args = parser.parse_args()

    kmer_dir = Path(args.kmer_dir)
    for txt_path in sorted(kmer_dir.rglob("*.kmer.txt")):
        fasta_path = txt_path.with_suffix(".fa")
        with txt_path.open() as in_handle, fasta_path.open("w") as out_handle:
            for i, line in enumerate(in_handle):
                kmer = line.rstrip("\n")
                out_handle.write(f">{i}\n{kmer}\n")
    Path(args.stamp).write_text("complete\n")


if __name__ == "__main__":
    main()
