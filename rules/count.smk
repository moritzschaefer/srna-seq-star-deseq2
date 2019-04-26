from os import path
import numpy as np


def feature_counts_extra(wildcards):
    extra = config["params"]["feature_counts"]
    # if is_paired:
    #     extra += " -p"
    return extra


rule index_bam:
    input:
        bam="star/{sample}-{unit}/Aligned.out.bam"
    output:
        bam="star/{sample}-{unit}/Aligned.sorted.bam",
        bai="star/{sample}-{unit}/Aligned.sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """samtools sort -o {output.bam} {input.bam}
        samtools index {output.bam}"""


rule feature_counts:
    input:
        # see STAR manual for additional output files
        bam="star/{sample}-{unit}/Aligned.sorted.bam",
        bai="star/{sample}-{unit}/Aligned.sorted.bam.bai",
        annotation="ref/{annotation}.gff3"
    output:
        counts="counts/{sample}-{unit}/{annotation}.featureCounts.txt",
        summary="qc/{sample}-{unit}/{annotation}.featureCounts.txt"
    log:
        "logs/feature_counts/{sample}-{unit}-{annotation}.txt"
    params:
        extra=feature_counts_extra
    threads:
        4
    wrapper:
        "https://raw.githubusercontent.com/moritzschaefer/srna-seq-star-deseq2/master/wrappers/subread/feature_counts"
        # "file:///home/schamori/moritzsphd/lib/srna-seq-star-deseq2/wrappers/subread/feature_counts"


rule merge_counts:
    input:
        expand("counts/{unit.sample}-{unit.unit}/{{annotation}}.featureCounts.txt", unit=units.itertuples())
    output:
        "counts/{annotation}.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule extract_snsnornas:
    input:
        counts="counts/gencode.vM20.annotation.tsv",
        annotation="ref/gencode.vM20.annotation.gff3"
    output:
        snrnas="counts/snrnas.tsv",
        snornas="counts/snornas.tsv",
        mt_trnas="counts/mt_trnas.tsv",
    script:
        "../scripts/split-filter-counts.py"

rule normalize_counts:
    input:
        "counts/{annotation}.tsv"
    output:
        "counts/{annotation}.cpm.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/cpm-counts.py"
