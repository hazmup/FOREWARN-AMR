# Directories------------------------------------------------------------------
configfile: "config.yaml"

import glob

# Define input and output directories
input_dir = config["base_dir"] + config["input_dir"]
output_dir = config["base_dir"] + config["output_dir"]
samples_file = config["input_list"]

# Get a list of samples from a samples file
with open(samples_file) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]
READS = ["1", "2"]


# Rules------------------------------------------------------------------------
rule all:
    input:
        # expand(output_dir + "/fastqc/{sample}_R{read}_fastqc.zip", sample=SAMPLES, read=READS),
        # expand(output_dir + "/trimmed/{sample}_R{read}_clean.fastq.gz", sample=SAMPLES, read=READS),
        # expand(output_dir + "/kraken2/{sample}_classification.tsv", sample=SAMPLES),
        # expand(output_dir + "/bracken/{sample}_bracken_" + config["bracken"]["level"] + ".tsv", sample=SAMPLES),
        # expand(output_dir + "/rgi/{sample}.txt", sample=SAMPLES),
        expand(output_dir + "/resfinder/{sample}/{sample}.json", sample=SAMPLES),


rule test:
    input:
        "/users/pa22/sgkionis/test.txt"
    output:
        "/users/pa22/sgkionis/test2.txt"
    shell:
        """
        cat {input} > {output}
        """


# Run FastQC on the fastq.gz files
rule fastqc_fastq_gz:
    input:
        input_dir + "/{sample}_{read}_001.fastq.gz",
    output:
        html = output_dir + "/fastqc/{sample}_{read}_fastqc.html",
        zip = output_dir + "/fastqc/{sample}_{read}_fastqc.zip",
    params: 
        "--quiet",
    log:
        output_dir+"/logs/fastqc/{sample}_{read}.log",
    params:
        threads=config["fastqc_params"],
    wrapper:
        "v1.25.0/bio/fastqc"


# Run FastQC on the fastq files
rule fastqc_fastq:
    input:
        input_dir + "/{sample}_{read}_001.fastq",
    output:
        html = output_dir + "/fastqc/{sample}_{read}_fastqc.html",
        zip = output_dir + "/fastqc/{sample}_{read}_fastqc.zip",
    params: 
        "--quiet",
    log:
        output_dir+"/logs/fastqc/{sample}_{read}.log",
    params:
        threads=config["fastqc_params"],
    wrapper:
        "v1.25.0/bio/fastqc"


# Run MultiQC on the FastQC reports
rule multiqc:
    input:
        expand(output_dir + "/fastqc/{sample}_R{read}_fastqc.zip", sample=SAMPLES, read=READS),
    output:
        output_dir + "/multiqc/multiqc_report.html",
    params:
        extra = "",
        use_input_files_only = True, # Optional, use only a.txt and don't search folder samtools_stats for files
    log:
        output_dir + "/logs/multiqc.log",
    wrapper:
        "v1.25.0/bio/multiqc"


# Run bbduk to remove adapters
rule remove_adapters:
    input:
        sample = [input_dir + "/{sample}_R1_001.fastq.gz", input_dir + "/{sample}_R2_001.fastq.gz"],
        adapters = config["bbduk"]["ref"]["remove_adapters"],
    output:
        trimmed = temp([output_dir+"/trimmed/{sample}_R1_adaptered.fastq.gz", output_dir + "/trimmed/{sample}_R2_adaptered.fastq.gz"])
    log:
        output_dir + "/logs/bbduk/pe/{sample}_remove_adapters.log",
    params:
        extra = lambda w, input: "ref={adapters} ktrim={ktrim} k={k} mink={mink} hdist={hdist} tpe tbo trimpolygright={trimpolygright} ftr={ftr}".format(
            adapters = input.adapters,
            ktrim = config["bbduk"]["ktrim"],
            k = config["bbduk"]["k"]["remove_adapters"],
            mink = config["bbduk"]["mink"],
            hdist = config["bbduk"]["hdist"],
            ftr = config["bbduk"]["ftr"],
            trimpolygright = config["bbduk"]["trimpolygright"],
        ),
        java_opts = "",
    resources:
        mem_mb = config["bbduk"]["memory"],
    threads:
        config["bbduk"]["threads"],
    wrapper:
        "v1.25.0/bio/bbtools/bbduk"


# Run bbduk to remove sequencing artifacts
rule remove_artifacts:
    input:
        sample = [output_dir + "/trimmed/{sample}_R1_adaptered.fastq.gz", output_dir + "/trimmed/{sample}_R2_adaptered.fastq.gz"],
        adapters = config["bbduk"]["ref"]["remove_artifacts"],
    output:
        trimmed = [output_dir + "/trimmed/{sample}_R1_clean.fastq.gz", output_dir + "/trimmed/{sample}_R2_clean.fastq.gz"],
        stats = output_dir + "/trimmed/{sample}.stats.txt",
    log:
        output_dir + "/logs/bbduk/pe/{sample}_remove_artifacts.log",
    params:
        extra = lambda w, input: "ref={ref} k={k} hdist={hdist} ".format(
            ref = input.adapters,
            ktrim = config["bbduk"]["ktrim"],
            k = config["bbduk"]["k"]["remove_artifacts"],
            hdist = config["bbduk"]["hdist"],
        ),
        java_opts = "",
    resources:
        mem_mb = config["bbduk"]["memory"],
    threads:
        config["bbduk"]["threads"],
    wrapper:
        "v1.25.0/bio/bbtools/bbduk"


# Run kraken2 to classify reads
rule kraken2:
    input:
        readF = output_dir + "/trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "/trimmed/{sample}_R2_clean.fastq.gz",
        db = config["kraken2"]["database"]
    output:
        rep = output_dir + "/kraken2/{sample}_kraken.tsv",
        classif = output_dir + "/kraken2/{sample}_classification.tsv",
    log:
        output_dir + "/logs/kraken2/{sample}.log",
    params:
        threads = config["kraken2"]["threads"],
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {input.db} --threads {params.threads} --paired {input.readF} {input.readR}  --report {output.rep} > {output.classif}

        """


# Run bracken on kraken2 reports to get clasification to specific taxonomic level
rule bracken:
    input:
        rep = output_dir + "/kraken2/{sample}_kraken.tsv",
        db = config["bracken"]["database"],
    output:
        brack_rep = output_dir + "/bracken/{sample}_bracken_" + config["bracken"]["level"] + ".tsv",
    log:
        output_dir + "/logs/bracken/{sample}.log",
    params:
        read_len = config["bracken"]["read_length"],
        threshold = config["bracken"]["threshold"],
        level = config["bracken"]["level"],
    conda:
        "envs/bracken.yaml",
    shell:
        """
        bracken -d {input.db} -i {input.rep} -o {output.brack_rep} -r {params.read_len} -l {params.level} -t {params.threshold}

        """


# Run rgi with CARD database to get resistome
rule rgi:
    input:
        readF = output_dir + "/trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "/trimmed/{sample}_R2_clean.fastq.gz",
        cardDb = config["rgi"]["cardDB"],
    output:
        cardTxt = output_dir + "/rgi/{sample}.txt",
    log:
        output_dir + "/logs/rgi/{sample}.log",
    params:
        aligner = config["rgi"]["aligner"],
        threads = config["rgi"]["threads"],
    conda:
        "envs/rgi.yaml"
    shell:
        """
        rgi load --card_json {input.cardDb}/card.json --local --card_annotation {input.cardDb}/card_database_v3.2.6.fasta \
--card_annotation_all_models {input.cardDb}/card_database_v3.2.6_all.fasta --wildcard_annotation {input.cardDb}/wildcard_database_v4.0.0.fasta \
--wildcard_annotation_all_models {input.cardDb}/wildcard_database_v4.0.0_all.fasta --wildcard_index {input.cardDb}/wildcard/index-for-model-sequences.txt \
--wildcard_version 4.0.0 --amr_kmers {input.cardDb}/wildcard/all_amr_61mers.txt --kmer_database {input.cardDb}/wildcard/61_kmer_db.json --kmer_size 61

        rgi bwt --read_one {input.readF} --read_two {input.readR} --aligner {params.aligner} --threads {params.threads} --output_file {output.cardTxt} --include_wildcard --include_other_models --local

        """


# Run resfinder to get resistome
rule resfinder:
    input:
        readF = output_dir + "/trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "/trimmed/{sample}_R2_clean.fastq.gz",
    output:
        output_dir + "/resfinder/{sample}/{sample}.json",
    log:
        output_dir + "/logs/rgi/{sample}.log",
    params:
        out_dir = output_dir + "/{sample}/",
        db_res = config["resfinder"]["db_res"],
    conda:
        "envs/resfinder.yaml",
    shell:
        """
        run_resfinder.py -ifq {input.readF} {input.readR} -o {params.out_dir} --acquired \
-db_res {params.db_res}
        """


# # Run resfinderfg to get resistome
# rule resfinderfg:
#     input:
#     output:
