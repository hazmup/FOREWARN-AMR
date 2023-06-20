# Directories------------------------------------------------------------------
configfile: "config.yaml"

import os
import pandas as pd
import glob

# Define input and output directories
input_dir = config["base_dir"] + config["input_dir"]
output_dir = config["base_dir"] + config["output_dir"]
samples_file = config["input_list"]

# Get a list of samples from a samples file
with open(samples_file) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]
READS = ["1", "2"]


# All rule - Modify accordigly-------------------------------------------------
rule all:
    input:
        output_dir + "multiqc/raw/multiqc_report.html",
        expand(output_dir + "bracken/{sample}_bracken_" + config["bracken"]["level"] + ".tsv", sample=SAMPLES),
        expand(output_dir + "rgi/{sample}.txt", sample=SAMPLES),
        expand(output_dir + "resfinderfg/{sample}_sorted.res", sample=SAMPLES),
        expand(output_dir + "resfinder/{sample}/resfinder_kma/{sample}_all.res", sample=SAMPLES),
        expand(output_dir + "resfinderfg/{sample}_filtered_normalized.res", sample=SAMPLES),
        expand(output_dir + "resfinder/{sample}_filtered_normalized.res", sample=SAMPLES),


# Rule definitions-------------------------------------------------------------

# Run FastQC on the fastq.gz files
rule fastqc_raw:
    input:
        input_dir + "{sample}_{read}_001.fastq.gz",
    output:
        html = output_dir + "fastqc/raw/{sample}_{read}_fastqc.html",
        zip = output_dir + "fastqc/raw/{sample}_{read}_fastqc.zip",
    params: 
        extra = "--quiet",
    log:
        output_dir + "logs/fastqc_raw/{sample}_{read}.log",
    threads: 1,
    resources:
        mem_mb = 1024,   
    wrapper:
        "v1.31.0/bio/fastqc"


# Run FastQC on the clean fastq.gz files
rule fastqc_clean:
    input:
        output_dir + "rqcfilter2/{sample}/{sample}_R1_001.anqdpht.fastq.gz",
    output:
        html = output_dir + "fastqc/clean/{sample}_clean_fastqc.html",
        zip = output_dir + "fastqc/clean/{sample}_clean_fastqc.zip",
    params: 
        extra = "--quiet",
    log:
        output_dir + "logs/fastqc_clean/{sample}.log",
    threads: 1,
    resources:
        mem_mb = 1024,   
    wrapper:
        "v1.31.0/bio/fastqc"


# Run MultiQC on the FastQC raw data reports
rule multiqc_raw:
    input:
        expand(output_dir + "fastqc/raw/{sample}_R{read}_fastqc.zip", sample=SAMPLES, read=READS),
    output:
        output_dir + "multiqc/raw/multiqc_report.html",
    params:
        extra = "",
        use_input_files_only = True, # Optional, use only a.txt and don't search folder samtools_stats for files
    log:
        output_dir + "logs/multiqc.log",
    wrapper:
        "v1.31.0/bio/multiqc"


# Run MultiQC on the FastQC clean data reports
rule multiqc_clean:
    input:
        expand(output_dir + "fastqc/clean/{sample}_clean_fastqc.zip", sample=SAMPLES),
    output:
        output_dir + "multiqc/clean/multiqc_report.html",
    params:
        extra = "",
        use_input_files_only = True, # Optional, use only a.txt and don't search folder samtools_stats for files
    log:
        output_dir + "logs/multiqc.log",
    wrapper:
        "v1.31.0/bio/multiqc"


# Performs standard QC per JGI metagenomics pipeline; does not remove bacterial contaminants, simply reports on them
rule rqcfilter2:
    input:
        readF = input_dir + "{sample}_R1_001.fastq.gz",
        readR = input_dir + "{sample}_R2_001.fastq.gz",
        ref_data = db_dir + config["rqcfilter2"]["database"],
    output:
        output_dir + "rqcfilter2/{sample}/{sample}_R1_001.anqdpht.fastq.gz",
    log:
        output_dir + "logs/rqcfilter2/{sample}.log",
    params:
        out_dir = output_dir + "rqcfilter2/{sample}/",
    resources:
        mem_gb = config["rqcfilter2"]["memory"],
    threads:
        config["rqcfilter2"]["threads"],
    conda:
        "envs/bbmap.yaml",
    shell:
        """
        rqcfilter2.sh -Xmx{resources.mem_gb}g chastityfilter=f jni=t in={input.readF} in2={input.readR} rqcfilterdata={input.ref_data} \
        path={params.out_dir} rna=f trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 \
        mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t \
        detectmicrobes=t sketch kapa=t clumpify=t tmpdir= barcodefilter=f trimpolyg=5 usejni=f \
        """


# Run kraken2 to classify reads
rule kraken2:
    input:
        readF = output_dir + "trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "trimmed/{sample}_R2_clean.fastq.gz",
        db = config["kraken2"]["database"]
    output:
        rep = output_dir + "kraken2/{sample}_kraken.tsv",
        classif = output_dir + "kraken2/{sample}_classification.tsv",
    log:
        output_dir + "logs/kraken2/{sample}.log",
    params:
        threads = config["kraken2"]["threads"],
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        kraken2 --db {input.db} --threads {params.threads} --paired {input.readF} {input.readR}  --report {output.rep} > {output.classif}
        """


# Filter kraken2 results to only keep bacteria
rule filter_kraken2:
    input:
        rep = output_dir + "kraken2/{sample}_kraken.tsv",
    output:
        bacteria_rep = output_dir + "kraken2/{sample}_kraken_bacteria.tsv",
    shell:
        """
        awk '$6=="D"{{if(++found==2)exit}} {{print $0}}' {input.rep} > {output.bacteria_rep}
        """


# Run bracken on the filtered kraken2 reports to get clasification to specific taxonomic level
rule bracken:
    input:
        rep = output_dir + "kraken2/{sample}_kraken_bacteria.tsv",
        db = db_dir + config["bracken"]["database"],
    output:
        brack_rep_temp = temp(output_dir + "bracken/{sample}_bracken_" + config["bracken"]["level"] + ".unsorted.tsv"),
        brack_rep = output_dir + "bracken/{sample}_bracken_" + config["bracken"]["level"] + ".tsv",
    log:
        output_dir + "logs/bracken/{sample}.log",
    params:
        read_len = config["bracken"]["read_length"],
        threshold = config["bracken"]["threshold"],
        level = config["bracken"]["level"],
    conda:
        "envs/bracken.yaml",
    shell:
        """
        bracken -d {input.db} -i {input.rep} -o {output.brack_rep_temp} -r {params.read_len} -l {params.level} -t {params.threshold}

        (head -n 1 {output.brack_rep_temp} && tail -n +2 {output.brack_rep_temp}  | sort -nr -k6 -t$'\t') > {output.brack_rep}
        """


# Run RGI with CARD database to get resistome
rule rgi:
    input:
        readF = output_dir + "trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "trimmed/{sample}_R2_clean.fastq.gz",
        cardDb = config["rgi"]["cardDB"],
    output:
        cardTxt = output_dir + "rgi/{sample}.txt",
    log:
        output_dir + "logs/rgi/{sample}.log",
    params:
        aligner = config["rgi"]["aligner"],
        threads = config["rgi"]["threads"],
    conda:
        "envs/rgi.yaml",
    shell:
        """
        rgi load --card_json {input.cardDb}/card.json --local --card_annotation {input.cardDb}/card_database_v3.2.6.fasta \
--card_annotation_all_models {input.cardDb}/card_database_v3.2.6_all.fasta --wildcard_annotation {input.cardDb}/wildcard_database_v4.0.0.fasta \
--wildcard_annotation_all_models {input.cardDb}/wildcard_database_v4.0.0_all.fasta --wildcard_index {input.cardDb}/wildcard/index-for-model-sequences.txt \
--wildcard_version 4.0.0 --amr_kmers {input.cardDb}/wildcard/all_amr_61mers.txt --kmer_database {input.cardDb}/wildcard/61_kmer_db.json --kmer_size 61

        rgi bwt --read_one {input.readF} --read_two {input.readR} --aligner {params.aligner} --threads {params.threads} --output_file {output.cardTxt} --include_wildcard --include_other_models --local
        """


# Run Resfinder to get resistome
rule resfinder:
    input:
        readF = output_dir + "trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "trimmed/{sample}_R2_clean.fastq.gz",
    output:
        res_temp = temp(output_dir + "resfinder/{sample}/resfinder_kma/{sample}_all.unsorted.res"),
        res = output_dir + "resfinder/{sample}/resfinder_kma/{sample}_all.res",
    log:
        output_dir + "logs/resfinder/{sample}.log",
    params:
        out_dir = output_dir + "resfinder/{sample}/",
        db_res = db_dir + config["resfinder"]["db_res"],
    conda:
        "envs/resfinder.yaml",
    shell:
        """
        # Checks if databases have been built and builds them otherwise
        if [ ! -e {params.db_res}/all.comp.b ]
        then
            {params.db_res}/INSTALL.py
        fi

        # Runs resfinder
        run_resfinder.py -ifq {input.readF} {input.readR} -o {params.out_dir} --acquired -db_res {params.db_res}

        # Create a summary results file and sorts it
        array=( {params.out_dir}resfinder_kma/*.res )
        {{ cat ${{array[@]:0:1}}; grep -vh "^#" ${{array[@]:1}}; }} > {output.res_temp}
        (head -n 1 {output.res_temp} && tail -n +2 {output.res_temp}  | sort -nr -k2 -t$'\t') > {output.res}
        """


# Run ResfinderFG to get resistome
rule resfinderfg:
    input:
        readF = output_dir + "trimmed/{sample}_R1_clean.fastq.gz",
        readR = output_dir + "trimmed/{sample}_R2_clean.fastq.gz",
    output:
        res_sorted = output_dir + "resfinderfg/{sample}_sorted.res",
        res = output_dir + "resfinderfg/{sample}.res",
    log:
        output_dir + "logs/resfinderfg/{sample}.log",
    params:
        out_prefix = output_dir + "resfinderfg/{sample}/",
        db_res = db_dir + config["resfinderfg"]["db_res"],
        threads = config["resfinderfg"]["threads"],
    threads:
        config["resfinderfg"]["threads"],
    conda:
        "envs/kma.yaml",
    shell:
        """
        # Checks if databases have been built and builds them otherwise
        if [ ! -e {params.db_res}/all.comp.b ]
        then
            {params.db_res}/INSTALL.py
        fi

        # Runs kma with resfinderfg db
        kma -ipe {input.readF} {input.readR} -o {params.out_prefix} -t_db {params.db_res}/all -t {params.threads}

        # Sorts the summary results file
        (head -n 1 {output.res} && tail -n +2 {output.res}  | sort -nr -k2 -t$'\t') > {output.res_sorted}
        """


# Filters Resfinder results based on 90% identity, 60% coverage, calculates Transcripts per Million, and sorts based on that
rule filter_and_normalize_resfinder:
    input:
        res_sorted = output_dir + "resfinder/{sample}/resfinder_kma/{sample}_all.res",
        ref = output_dir + "multiqc/clean/multiqc_report_data/multiqc_fastqc.txt",
    output:
        res = output_dir + "resfinder/{sample}_filtered_normalized.res"
    run:
        # Reads files
        ref = pd.read_csv(input.ref, sep="\t")
        res = pd.read_csv(input.res_sorted, sep="\t")
        
        # Filter results using Strict criteria of 60% coverage and 90% identity
        res = res.query('Template_Coverage >= 60 & Query_Identity >= 90')

        # Calculates TPM by first normalizing for gene length, and then by total sequencing depth (not just ARG reads)
        new_res = res
        new_res['TPM'] = new_res['Depth'] / (new_res['Template_length']/1000) * ((ref.query('Sample == "'+wildcards.sample+'_R1_001.anqdpht"').iat[0, 4])/1000000)

        # Saves new results
        new_res.to_csv(output.res, sep="\t", index=False)


# Filters ResfinderFG results based on 90% identity, 60% coverage, calculates Transcripts per Million, and sorts based on that
rule filter_and_normalize_resfinder_fg:
    input:
        res_sorted = output_dir + "resfinderfg/{sample}_sorted.res",
        ref = output_dir + "multiqc/clean/multiqc_report_data/multiqc_fastqc.txt",
    output:
        res = output_dir + "resfinderfg/{sample}_filtered_normalized.res"
    run:
        # Reads files
        ref = pd.read_csv(input.ref, sep="\t")
        res = pd.read_csv(input.res_sorted, sep="\t")
        
        # Filter results using Strict criteria of 60% coverage and 90% identity
        res = res.query('Template_Coverage >= 60 & Query_Identity >= 90')

        # Calculates TPM by first normalizing for gene length, and then by total sequencing depth (not just ARG reads)
        new_res = res
        new_res['TPM'] = new_res['Depth'] / (new_res['Template_length']/1000) * ((ref.query('Sample == "'+wildcards.sample+'_R1_001.anqdpht"').iat[0, 4])/1000000)

        # Saves new results
        new_res.to_csv(output.res, sep="\t", index=False)
