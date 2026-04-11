import json
from pathlib import Path

configfile: "config.yaml"


OUTPUT_ROOT = config["output_root"]
EXON_K = config["params"]["exon_k"]
UTR_K = config["params"]["utr_k"]
MISMATCH = config["params"]["mismatch"]
MAX_CHR = config["params"]["max_chr"]
MAX_GENE_ALIGNMENT = config["params"]["max_gene_alignment"]
DIR_NAME_LEN = config["params"]["dir_name_len"]
VERBOSE = config["params"]["verbose"]
CROSSMAP_GTF = f"{OUTPUT_ROOT}/annot/annot.crossmap.gtf"
ANNOT_TXT = f"{OUTPUT_ROOT}/annot/annot.exon_utr.txt"
ANNOT_SUMMARY = f"{OUTPUT_ROOT}/annot/annot.summary.json"
GENOME_SPLIT_DIR = f"{OUTPUT_ROOT}/reference/genome_split"
GENOME_SPLIT_MANIFEST = f"{OUTPUT_ROOT}/reference/genome_split.chromosomes.txt"
BOWTIE_PREFIX = f"{OUTPUT_ROOT}/reference/bowtie_index/bowtie_index"
GENMAP_BIN = config["mappability"]["genmap_bin"]
GENMAP_INDEX_DIR = f"{OUTPUT_ROOT}/reference/genmap_index"
GENMAP_THREADS = config["mappability"]["threads"]
MAPPABILITY_DIR = f"{OUTPUT_ROOT}/reference/mappability"
BEDGRAPH_TEMPLATE = f"{MAPPABILITY_DIR}/rat_{{k}}mer_mappability.bedgraph"
EXON_BEDGRAPH = BEDGRAPH_TEMPLATE.format(k=EXON_K)
UTR_BEDGRAPH = BEDGRAPH_TEMPLATE.format(k=UTR_K)
GENE_MAPPABILITY = f"{OUTPUT_ROOT}/gene_mappability/gene_mappability.txt"
GENE_MAPPABILITY_INFO = f"{GENE_MAPPABILITY}.info"
AMBIGUOUS_KMERS_DIR = f"{OUTPUT_ROOT}/ambiguous_kmers"
KMER_FASTA_READY = f"{OUTPUT_ROOT}/flags/ambiguous_kmers_fasta_complete.txt"
ALIGNMENT_DIR = f"{OUTPUT_ROOT}/ambiguous_kmers_alignment"
CROSSMAP_DIR = f"{OUTPUT_ROOT}/cross_mappability"
CROSSMAP_BATCH_DIR = f"{CROSSMAP_DIR}/batches"
COMBINED_CROSSMAP = f"{OUTPUT_ROOT}/cross_mappability.tsv.gz"
POS2GENE_READY = f"{CROSSMAP_DIR}/pos2gene/.init_complete"

VALIDATE_DIR = f"{OUTPUT_ROOT}/validate"
HUMAN_VALIDATE_CROSSMAP = f"{VALIDATE_DIR}/hg38_cross_mappability_strength.txt.gz"
VALIDATION_SUMMARY_TSV = f"{VALIDATE_DIR}/crossmap_validation_summary.tsv"
VALIDATION_COMPARISON_TSV = f"{VALIDATE_DIR}/crossmap_validation_comparison.tsv"
VALIDATION_MAPPABILITY_DECILES_TSV = f"{VALIDATE_DIR}/crossmap_validation_mappability_deciles.tsv"
VALIDATION_REPORT = f"{VALIDATE_DIR}/crossmap_validation.md"
VALIDATION_SUMMARY_JSON = f"{VALIDATE_DIR}/crossmap_validation_summary.json"
ORTHOLOG_CACHE_DIR = f"{VALIDATE_DIR}/ensembl_human_rat_ortholog_cache"
ORTHOLOG_TABLE = f"{VALIDATE_DIR}/ensembl_human_rat_orthologs.tsv.gz"
ORTHOLOG_COMPARISON_TSV = f"{VALIDATE_DIR}/crossmap_ortholog_comparison.tsv"
ORTHOLOG_COMPARISON_JSON = f"{VALIDATE_DIR}/crossmap_ortholog_comparison.json"
ORTHOLOG_COMPARISON_REPORT = f"{VALIDATE_DIR}/crossmap_ortholog_comparison.md"

BOWTIE_INDEX_FILES = [
    f"{BOWTIE_PREFIX}.1.ebwt",
    f"{BOWTIE_PREFIX}.2.ebwt",
    f"{BOWTIE_PREFIX}.3.ebwt",
    f"{BOWTIE_PREFIX}.4.ebwt",
    f"{BOWTIE_PREFIX}.rev.1.ebwt",
    f"{BOWTIE_PREFIX}.rev.2.ebwt",
]


def crossmap_batch_ranges_from_gene_count(gene_count):
    return [
        (n1, min(n1 + MAX_GENE_ALIGNMENT - 1, gene_count))
        for n1 in range(1, gene_count + 1, MAX_GENE_ALIGNMENT)
    ]


def crossmap_batch_stamp(n1, n2):
    return f"{CROSSMAP_BATCH_DIR}/n1_{n1}_n2_{n2}.complete"


def crossmap_batch_outputs(wildcards):
    checkpoint_output = checkpoints.preprocess_annotation.get(**wildcards).output
    with open(checkpoint_output.summary) as handle:
        gene_count = json.load(handle)["genes_with_exons"]
    return [
        crossmap_batch_stamp(n1, n2)
        for n1, n2 in crossmap_batch_ranges_from_gene_count(gene_count)
    ]

rule all:
    input:
        CROSSMAP_GTF,
        ANNOT_TXT,
        EXON_BEDGRAPH,
        UTR_BEDGRAPH,
        GENE_MAPPABILITY,
        GENE_MAPPABILITY_INFO,
        COMBINED_CROSSMAP,


rule split_reference_fasta:
    input:
        fasta=config["reference"]["fasta"]
    output:
        directory(GENOME_SPLIT_DIR),
        GENOME_SPLIT_MANIFEST
    shell:
        """
        mkdir -p {GENOME_SPLIT_DIR}
        python3 scripts/setup/split_fasta_by_chrom.py \
            --fasta {input.fasta} \
            --out-dir {GENOME_SPLIT_DIR} \
            --manifest {GENOME_SPLIT_MANIFEST}
        """


rule build_bowtie_index:
    input:
        fasta=config["reference"]["fasta"]
    output:
        BOWTIE_INDEX_FILES
    params:
        prefix=BOWTIE_PREFIX
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        bowtie-build {input.fasta} {params.prefix}
        """


rule build_genmap_index:
    input:
        fasta=config["reference"]["fasta"]
    output:
        directory(GENMAP_INDEX_DIR)
    params:
        index_parent=lambda wildcards, output: str(Path(str(output[0])).parent)
    resources:
        mem_mb=32000,
    shell:
        """
        mkdir -p {params.index_parent}
        {GENMAP_BIN} index \
            -F {input.fasta} \
            -I {output}
        """


checkpoint preprocess_annotation:
    input:
        gtf=config["reference"]["gtf"]
    output:
        gtf=CROSSMAP_GTF,
        summary=ANNOT_SUMMARY
    shell:
        """
        mkdir -p $(dirname {output.gtf})
        python3 scripts/setup/prepare_refseq_crossmap.py \
            --gtf {input.gtf} \
            --out-crossmap-gtf {output.gtf} \
            --summary-json {output.summary}
        """


rule gtf_to_txt:
    input:
        gtf=CROSSMAP_GTF
    output:
        txt=ANNOT_TXT
    shell:
        """
        mkdir -p $(dirname {output.txt})
        Rscript scripts/crossmap/gtf_to_txt.R \
            -gtf {input.gtf} \
            -f exon,UTR \
            -o {output.txt}
        """


rule genmap_bedgraph:
    input:
        index=GENMAP_INDEX_DIR
    output:
        bedgraph=BEDGRAPH_TEMPLATE
    params:
        prefix=lambda wildcards, output: str(output.bedgraph).removesuffix(".bedgraph")
    resources:
        mem_mb=32000,
    threads: GENMAP_THREADS
    shell:
        """
        mkdir -p $(dirname {output.bedgraph})
        {GENMAP_BIN} map \
            -I {GENMAP_INDEX_DIR} \
            -O {params.prefix} \
            -E {MISMATCH} \
            -K {wildcards.k} \
            -bg \
            -T {GENMAP_THREADS}
        """


rule compute_mappability:
    input:
        annot=ANNOT_TXT,
        exon_bedgraph=EXON_BEDGRAPH,
        utr_bedgraph=UTR_BEDGRAPH
    output:
        txt=GENE_MAPPABILITY,
        info=GENE_MAPPABILITY_INFO
    resources:
        mem_mb=64000,
    shell:
        """
        mkdir -p $(dirname {output.txt})
        Rscript scripts/crossmap/compute_mappability.R \
            -annot {input.annot} \
            -k_exon {EXON_K} \
            -k_utr {UTR_K} \
            -kmap_exon {input.exon_bedgraph} \
            -kmap_utr {input.utr_bedgraph} \
            -verbose {VERBOSE} \
            -o {output.txt}
        """


rule generate_ambiguous_kmers:
    input:
        annot=ANNOT_TXT,
        mappability=GENE_MAPPABILITY,
        genome_dir=GENOME_SPLIT_DIR,
        exon_bedgraph=EXON_BEDGRAPH,
        utr_bedgraph=UTR_BEDGRAPH
    output:
        directory(AMBIGUOUS_KMERS_DIR)
    resources:
        mem_mb=64000,
    shell:
        """
        mkdir -p {AMBIGUOUS_KMERS_DIR}
        Rscript scripts/crossmap/generate_ambiguous_kmers.R \
            -mappability {input.mappability} \
            -genome {input.genome_dir} \
            -annot {input.annot} \
            -k_exon {EXON_K} \
            -k_utr {UTR_K} \
            -kmap_exon {input.exon_bedgraph} \
            -kmap_utr {input.utr_bedgraph} \
            -th1 0 \
            -th2 1 \
            -dir_name_len {DIR_NAME_LEN} \
            -verbose {VERBOSE} \
            -o {AMBIGUOUS_KMERS_DIR}
        """


rule init_crossmap_resources:
    input:
        annot=ANNOT_TXT,
        mappability=GENE_MAPPABILITY,
        kmers=KMER_FASTA_READY,
        index=BOWTIE_INDEX_FILES
    output:
        POS2GENE_READY
    shell:
        """
        mkdir -p {ALIGNMENT_DIR} {CROSSMAP_DIR} $(dirname {POS2GENE_READY})
        Rscript scripts/crossmap/compute_cross_mappability.R \
            -annot {input.annot} \
            -mappability {input.mappability} \
            -kmer {AMBIGUOUS_KMERS_DIR} \
            -align {ALIGNMENT_DIR} \
            -index {BOWTIE_PREFIX} \
            -mismatch {MISMATCH} \
            -max_chr {MAX_CHR} \
            -max_gene {MAX_GENE_ALIGNMENT} \
            -initonly TRUE \
            -dir_name_len {DIR_NAME_LEN} \
            -verbose {VERBOSE} \
            -o {CROSSMAP_DIR}
        touch {POS2GENE_READY}
        """


rule make_kmer_fastas:
    input:
        AMBIGUOUS_KMERS_DIR
    output:
        KMER_FASTA_READY
    shell:
        """
        python3 scripts/setup/make_kmer_fastas.py \
            --kmer-dir {input} \
            --stamp {output}
        """


rule compute_crossmap:
    input:
        ready=POS2GENE_READY,
        annot=ANNOT_TXT,
        mappability=GENE_MAPPABILITY,
        kmers=KMER_FASTA_READY
    output:
        touch(f"{CROSSMAP_BATCH_DIR}/n1_{{n1}}_n2_{{n2}}.complete")
    resources:
        mem_mb=64000,
        runtime='8h',
    shell:
        """
        mkdir -p {ALIGNMENT_DIR} {CROSSMAP_DIR} {CROSSMAP_BATCH_DIR}
        Rscript scripts/crossmap/compute_cross_mappability.R \
            -annot {input.annot} \
            -mappability {input.mappability} \
            -kmer {AMBIGUOUS_KMERS_DIR} \
            -align {ALIGNMENT_DIR} \
            -index {BOWTIE_PREFIX} \
            -n1 {wildcards.n1} \
            -n2 {wildcards.n2} \
            -mismatch {MISMATCH} \
            -max_chr {MAX_CHR} \
            -max_gene {MAX_GENE_ALIGNMENT} \
            -initonly FALSE \
            -dir_name_len {DIR_NAME_LEN} \
            -verbose {VERBOSE} \
            -o {CROSSMAP_DIR}
        touch {output}
        """


rule combine_crossmap:
    input:
        crossmap_batch_outputs
    output:
        COMBINED_CROSSMAP
    shell:
        """
        python3 scripts/setup/combine_crossmap.py \
            --crossmap-dir {CROSSMAP_DIR} \
            --output {output}
        """


rule validate_crossmap:
    input:
        rat_crossmap=COMBINED_CROSSMAP,
        human_crossmap=HUMAN_VALIDATE_CROSSMAP,
        rat_mappability=GENE_MAPPABILITY
    output:
        summary_tsv=VALIDATION_SUMMARY_TSV,
        comparison_tsv=VALIDATION_COMPARISON_TSV,
        mappability_deciles_tsv=VALIDATION_MAPPABILITY_DECILES_TSV,
        report=VALIDATION_REPORT,
        summary_json=VALIDATION_SUMMARY_JSON
    shell:
        """
        python3 scripts/validate/summarize_crossmap_validation.py \
            --rat-crossmap {input.rat_crossmap} \
            --human-crossmap {input.human_crossmap} \
            --rat-mappability {input.rat_mappability} \
            --summary-tsv {output.summary_tsv} \
            --comparison-tsv {output.comparison_tsv} \
            --mappability-deciles-tsv {output.mappability_deciles_tsv} \
            --report-md {output.report} \
            --summary-json {output.summary_json}
        """


rule fetch_ensembl_human_rat_orthologs:
    input:
        human_crossmap=HUMAN_VALIDATE_CROSSMAP
    output:
        ORTHOLOG_TABLE
    shell:
        """
        python3 scripts/validate/fetch_ensembl_human_rat_orthologs.py \
            --human-crossmap {input.human_crossmap} \
            --output {output} \
            --cache-dir {ORTHOLOG_CACHE_DIR}
        """


rule compare_crossmap_orthologs:
    input:
        rat_crossmap=COMBINED_CROSSMAP,
        human_crossmap=HUMAN_VALIDATE_CROSSMAP,
        orthologs=ORTHOLOG_TABLE
    output:
        joined_tsv=ORTHOLOG_COMPARISON_TSV,
        summary_json=ORTHOLOG_COMPARISON_JSON,
        report=ORTHOLOG_COMPARISON_REPORT
    shell:
        """
        python3 scripts/validate/compare_crossmap_orthologs.py \
            --rat-crossmap {input.rat_crossmap} \
            --human-crossmap {input.human_crossmap} \
            --orthologs {input.orthologs} \
            --joined-tsv {output.joined_tsv} \
            --summary-json {output.summary_json} \
            --report-md {output.report}
        """
