#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--crossmap-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    crossmap_dir = Path(args.crossmap_dir)
    paths = sorted(crossmap_dir.rglob("*.crossmap.txt"))
    with gzip.open(args.output, "wt") as out_handle:
        # out_handle.write("G1\tG2\tS1\n")
        for path in paths:
            with path.open() as in_handle:
                for line in in_handle:
                    out_handle.write(line)


if __name__ == "__main__":
    main()
