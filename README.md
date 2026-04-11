# cross-mappability

Rat cross-mappability workflow built around the original Battle lab `crossmap`
method for trans-eQTL alignment-confounding analyses.

## What this repo does

This workflow takes:

- a whole-genome reference FASTA
- a RefSeq-style GTF

and produces:

- a rat exon/UTR annotation table
- exon and UTR k-mer mappability bedGraphs
- gene-level mappability scores
- ambiguous k-mer intermediates
- final gene-to-gene cross-mappability calls

The upstream `battle-lab/crossmap` logic is included as a git submodule under
[`scripts/crossmap`](scripts/crossmap)
and is used essentially as-is. Rat-specific logic is confined to preprocessing
and workflow orchestration.

## Rat-specific adaptation

The current rat RefSeq GTF does not contain explicit `UTR` features. This repo
derives transcript UTR intervals as `exon - CDS`, writes an augmented
crossmap-compatible GTF, and then runs the original upstream `gtf_to_txt.R`
against that preprocessed GTF.

The default rat inputs in [`config.yaml`](config.yaml)
point to:

- `../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.fa`
- `../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf`

## Workflow

Run the full workflow with:

```bash
TMPDIR=/tmp XDG_CACHE_HOME=/tmp snakemake --cores 1
```

The default targets are:

- `results/rat/annot/annot.crossmap.gtf`
- `results/rat/annot/annot.exon_utr.txt`
- `results/rat/gene_mappability/gene_mappability.txt`
- `results/rat/cross_mappability.tsv.gz`

## Configuration

Edit [`config.yaml`](config.yaml) to set:

- reference FASTA and GTF
- split-genome output directory
- Bowtie index prefix
- GenMap binary, index directory, and thread count
- k values and mismatch setting
- output root

This repo now generates the exon and UTR mappability bedGraphs with `GenMap`
from the reference FASTA. The default config expects `genmap` to be available
on `PATH`, but you can change `mappability.genmap_bin` if it is installed
elsewhere. The generated bedGraphs are written under
`results/rat/reference/mappability/` using the pattern
`rat_{k}mer_mappability.bedgraph`.

## Dependencies

- `snakemake`
- `python3`
- `Rscript`
- `bowtie`
- `bowtie-build`
- `genmap`
- R packages: `argparser=0.4`, `data.table`, `intervals`, `seqinr`, `stringr`

## Tests

Python unit and smoke tests:

```bash
python -m unittest test.test_prepare_refseq_crossmap
python -m unittest test.test_smoke
```

## Validation

If `data/validate/hg38_cross_mappability_strength.txt.gz` is present, the
workflow can generate a validation report at
`data/validate/crossmap_validation.md` with:

```bash
snakemake --cores 1 data/validate/crossmap_validation.md
```

An ortholog-based gene-level comparison can be generated with:

```bash
snakemake --cores 1 data/validate/crossmap_ortholog_comparison.md
```

That report does two things:

- compares the rat cross-mappability graph against the downloaded human
  reference using graph-level summaries that do not require ortholog mapping
- checks internal consistency in rat by confirming that low-mappability genes
  carry much more outgoing cross-mappability burden than high-mappability genes

This is intended as a first-pass reassurance step. A stricter species-to-species
validation would require an explicit ortholog table.

The ortholog-based report fetches human-rat orthologs from Ensembl REST,
restricts to one-to-one orthologs, and compares rat and human gene-level
cross-mappability burdens on the shared set.
