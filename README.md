# cross-mappability

Rat cross-mappability workflow built around the original Battle lab `crossmap`
method for trans-eQTL alignment-confounding analyses.

## What this repo does

This workflow takes:

- a whole-genome reference FASTA
- a RefSeq-style GTF
- exon and UTR k-mer mappability bedGraphs

and produces:

- a rat exon/UTR annotation table
- gene-level mappability scores
- ambiguous k-mer intermediates
- final gene-to-gene cross-mappability calls
- run metadata

The upstream `battle-lab/crossmap` logic is included as a git submodule under
[`scripts/crossmap`](/data/hps/home/dmunro/gdml/cross-mappability/scripts/crossmap)
and is used essentially as-is. Rat-specific logic is confined to preprocessing
and workflow orchestration.

## Rat-specific adaptation

The current rat RefSeq GTF does not contain explicit `UTR` features. This repo
derives transcript UTR intervals as `exon - CDS`, writes an augmented
crossmap-compatible GTF, and then runs the original upstream `gtf_to_txt.R`
against that preprocessed GTF.

The default rat inputs in [`config/config.yaml`](/data/hps/home/dmunro/gdml/cross-mappability/config/config.yaml)
point to:

- [`../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.fa`](/data/hps/home/dmunro/gdml/ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.fa)
- [`../ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf`](/data/hps/home/dmunro/gdml/ratgtex/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf)

## Workflow

Run the full workflow with:

```bash
TMPDIR=/tmp XDG_CACHE_HOME=/tmp snakemake --cores 1
```

The default targets are:

- `results/rat/annot/annot.crossmap.gtf`
- `results/rat/annot/annot.exon_utr.txt`
- `results/rat/gene_mappability/gene_mappability.txt`
- `results/rat/cross_mappability.tsv`
- `results/rat/run_metadata.json`

## Configuration

Edit [`config/config.yaml`](/data/hps/home/dmunro/gdml/cross-mappability/config/config.yaml) to set:

- reference FASTA and GTF
- split-genome output directory
- Bowtie index prefix
- exon and UTR mappability bedGraphs
- k values and mismatch setting
- output root

The mappability bedGraphs are required inputs. This repo does not generate them
automatically. If rat bedGraphs are unavailable, generate them separately with
GEM or another equivalent workflow, then point the config at those files.

## Dependencies

- `snakemake`
- `python3`
- `Rscript`
- `bowtie`
- `bowtie-build`
- R packages: `argparser`, `data.table`, `intervals`, `seqinr`, `stringr`

## Tests

Python unit and smoke tests:

```bash
python -m unittest test.test_prepare_refseq_crossmap
python -m unittest test.test_smoke
```

The smoke test is dependency-gated and skips unless the original crossmap R
packages and Bowtie executables are installed.
