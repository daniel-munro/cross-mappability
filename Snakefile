configfile: "config/config.yaml"


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
GENOME_SPLIT_DIR = config["reference"]["split_dir"]
GENOME_SPLIT_MANIFEST = config["reference"]["split_manifest"]
BOWTIE_PREFIX = config["reference"]["bowtie_index_prefix"]
GENE_MAPPABILITY = f"{OUTPUT_ROOT}/gene_mappability/gene_mappability.txt"
GENE_MAPPABILITY_INFO = f"{GENE_MAPPABILITY}.info"
AMBIGUOUS_KMERS_DIR = f"{OUTPUT_ROOT}/ambiguous_kmers"
KMER_FASTA_READY = f"{OUTPUT_ROOT}/ambiguous_kmers/.fasta_complete"
ALIGNMENT_DIR = f"{OUTPUT_ROOT}/ambiguous_kmers_alignment"
CROSSMAP_DIR = f"{OUTPUT_ROOT}/cross_mappability"
COMBINED_CROSSMAP = f"{OUTPUT_ROOT}/cross_mappability.tsv"
RUN_METADATA = f"{OUTPUT_ROOT}/run_metadata.json"
POS2GENE_READY = f"{CROSSMAP_DIR}/pos2gene/.init_complete"
CROSSMAP_COMPLETE = f"{CROSSMAP_DIR}/.complete"
CONFIG_PATH = workflow.overwrite_configfiles[0] if workflow.overwrite_configfiles else "config/config.yaml"

BOWTIE_INDEX_FILES = [
    f"{BOWTIE_PREFIX}.1.ebwt",
    f"{BOWTIE_PREFIX}.2.ebwt",
    f"{BOWTIE_PREFIX}.3.ebwt",
    f"{BOWTIE_PREFIX}.4.ebwt",
    f"{BOWTIE_PREFIX}.rev.1.ebwt",
    f"{BOWTIE_PREFIX}.rev.2.ebwt",
]


rule all:
    input:
        CROSSMAP_GTF,
        ANNOT_TXT,
        GENE_MAPPABILITY,
        GENE_MAPPABILITY_INFO,
        COMBINED_CROSSMAP,
        RUN_METADATA,


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


rule preprocess_annotation:
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


rule compute_mappability:
    input:
        annot=ANNOT_TXT,
        exon_bedgraph=config["mappability"]["exon_bedgraph"],
        utr_bedgraph=config["mappability"]["utr_bedgraph"]
    output:
        txt=GENE_MAPPABILITY,
        info=GENE_MAPPABILITY_INFO
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
        exon_bedgraph=config["mappability"]["exon_bedgraph"],
        utr_bedgraph=config["mappability"]["utr_bedgraph"]
    output:
        directory(AMBIGUOUS_KMERS_DIR)
    shell:
        """
        mkdir -p {output}
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
            -o {output}
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
        mkdir -p {ALIGNMENT_DIR}
        mkdir -p {CROSSMAP_DIR}
        mkdir -p $(dirname {output})
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
        touch {output}
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
        CROSSMAP_COMPLETE
    shell:
        """
        mkdir -p {ALIGNMENT_DIR}
        mkdir -p {CROSSMAP_DIR}
        Rscript scripts/crossmap/compute_cross_mappability.R \
            -annot {input.annot} \
            -mappability {input.mappability} \
            -kmer {AMBIGUOUS_KMERS_DIR} \
            -align {ALIGNMENT_DIR} \
            -index {BOWTIE_PREFIX} \
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
        CROSSMAP_COMPLETE
    output:
        COMBINED_CROSSMAP
    shell:
        """
        python3 scripts/setup/combine_crossmap.py \
            --crossmap-dir {CROSSMAP_DIR} \
            --output {output}
        """


rule write_run_metadata:
    output:
        RUN_METADATA
    shell:
        """
        mkdir -p $(dirname {output})
        python3 scripts/setup/write_run_metadata.py \
            --config {CONFIG_PATH} \
            --output {output}
        """
