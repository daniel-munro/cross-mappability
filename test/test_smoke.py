import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


class WorkflowSmokeTest(unittest.TestCase):
    def test_toy_workflow_completes(self):
        repo = Path.cwd()
        crossmap_bin = repo.parent / "miniforge3" / "envs" / "crossmap" / "bin"
        with tempfile.TemporaryDirectory() as tmpdir:
            workdir = Path(tmpdir)
            output_root = workdir / "results"
            config_path = workdir / "config.yaml"
            config_path.write_text(
                textwrap.dedent(
                    f"""\
                    reference:
                      fasta: {repo / "test/data/toy.fa"}
                      gtf: {repo / "test/data/toy.gtf"}
                      split_dir: {output_root / "reference/genome_split"}
                      split_manifest: {output_root / "reference/genome_split.chromosomes.txt"}
                      bowtie_index_prefix: {output_root / "reference/bowtie_index/toy"}

                    mappability:
                      genmap_bin: {crossmap_bin / "genmap"}
                      index_dir: {output_root / "reference/genmap_index"}
                      threads: 1

                    params:
                      exon_k: 5
                      utr_k: 3
                      mismatch: 1
                      max_chr: 1
                      max_gene_alignment: 10
                      dir_name_len: 5
                      verbose: 0

                    output_root: {output_root}
                    """
                )
            )

            env = os.environ.copy()
            env["PATH"] = f"{crossmap_bin}:{env['PATH']}"
            env["TMPDIR"] = "/tmp"
            env["HOME"] = "/tmp"
            env["XDG_CACHE_HOME"] = "/tmp"
            env["XDG_CONFIG_HOME"] = "/tmp"

            subprocess.run(
                [
                    "snakemake",
                    "--snakefile",
                    str(repo / "Snakefile"),
                    "--configfile",
                    str(config_path),
                    "--cores",
                    "1",
                ],
                cwd=repo,
                env=env,
                check=True,
            )

            combined = output_root / "cross_mappability.tsv"
            exon_bedgraph = output_root / "reference/mappability/rat_5mer_mappability.bedgraph"
            utr_bedgraph = output_root / "reference/mappability/rat_3mer_mappability.bedgraph"

            self.assertTrue(exon_bedgraph.exists())
            self.assertTrue(utr_bedgraph.exists())
            self.assertTrue(combined.exists())

            crossmap_lines = combined.read_text().splitlines()
            self.assertEqual(crossmap_lines[0], "G1\tG2\tS1")
            self.assertGreater(len(crossmap_lines), 1)


if __name__ == "__main__":
    unittest.main()
