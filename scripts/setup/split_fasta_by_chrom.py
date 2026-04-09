#!/usr/bin/env python3

import argparse
from pathlib import Path


def write_record(out_dir: Path, header: str | None, seq_lines: list[str]) -> None:
    if header is None:
        return
    out_path = out_dir / f"{header}.fa"
    with out_path.open("w") as handle:
        handle.write(f">{header}\n")
        sequence = "".join(seq_lines)
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i : i + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--manifest", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    header = None
    seq_lines: list[str] = []
    chroms: list[str] = []
    with open(args.fasta) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                write_record(out_dir, header, seq_lines)
                header = line[1:].split()[0]
                chroms.append(header)
                seq_lines = []
            else:
                seq_lines.append(line)
    write_record(out_dir, header, seq_lines)
    Path(args.manifest).write_text("".join(f"{chrom}\n" for chrom in chroms))


if __name__ == "__main__":
    main()
