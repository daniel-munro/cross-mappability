import unittest
from pathlib import Path
import tempfile

from scripts.setup.prepare_refseq_crossmap import (
    load_transcripts,
    transcript_feature_rows,
    write_gtf,
)


class PrepareRefseqCrossmapTest(unittest.TestCase):
    def test_derives_utrs_from_exon_minus_cds(self):
        transcripts = load_transcripts("test/data/toy.gtf", gene_filter=None)
        rows, summary = transcript_feature_rows(transcripts)
        gene_a_utrs = [
            (row["start"], row["end"])
            for row in rows
            if row["attributes"]["gene_id"] == "geneA" and row["feature"] == "UTR"
        ]
        self.assertEqual(gene_a_utrs, [(1, 10), (36, 40), (50, 54), (76, 80)])
        self.assertEqual(summary["genes_with_utrs"], 2)

    def test_noncoding_transcript_has_no_utrs(self):
        transcripts = load_transcripts("test/data/toy.gtf", gene_filter=None)
        rows, _ = transcript_feature_rows(transcripts)
        gene_c_utrs = [
            row
            for row in rows
            if row["attributes"]["gene_id"] == "geneC" and row["feature"] == "UTR"
        ]
        gene_c_exons = [
            row
            for row in rows
            if row["attributes"]["gene_id"] == "geneC" and row["feature"] == "exon"
        ]
        self.assertEqual(gene_c_utrs, [])
        self.assertEqual(len(gene_c_exons), 1)

    def test_gene_filter_limits_output(self):
        transcripts = load_transcripts("test/data/toy.gtf", gene_filter={"geneB"})
        rows, _ = transcript_feature_rows(transcripts)
        self.assertEqual({row["attributes"]["gene_id"] for row in rows}, {"geneB"})

    def test_writes_crossmap_compatible_gtf(self):
        transcripts = load_transcripts("test/data/toy.gtf", gene_filter=None)
        rows, _ = transcript_feature_rows(transcripts)
        with tempfile.NamedTemporaryFile("w", delete=False) as handle:
            out_path = Path(handle.name)
        try:
            write_gtf(rows, str(out_path))
            content = out_path.read_text()
            self.assertIn('gene_name "geneA";', content)
            self.assertIn('gene_type "protein_coding";', content)
            self.assertIn('\tUTR\t', content)
        finally:
            out_path.unlink(missing_ok=True)


if __name__ == "__main__":
    unittest.main()
