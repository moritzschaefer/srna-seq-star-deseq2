# TODO adapt!! code from https://github.com/jrderuiter/snakemake-rnaseq
from os import path

import numpy as np

################################################################################
# Functions                                                                    #
################################################################################

def normalize_counts(counts):
    """Normalizes expression counts using DESeq's median-of-ratios approach."""

    with np.errstate(divide="ignore"):
        size_factors = estimate_size_factors(counts)
        return counts / size_factors


def estimate_size_factors(counts):
    """Calculate size factors for DESeq's median-of-ratios normalization."""

    def _estimate_size_factors_col(counts, log_geo_means):
        log_counts = np.log(counts)
        mask = np.isfinite(log_geo_means) & (counts > 0)
        return np.exp(np.median((log_counts - log_geo_means)[mask]))

    log_geo_means = np.mean(np.log(counts), axis=1)
    size_factors = np.apply_along_axis(
        _estimate_size_factors_col, axis=0,
        arr=counts, log_geo_means=log_geo_means)

    return size_factors


################################################################################
# Rules                                                                        #
################################################################################


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
    output:
        counts="counts/{sample}-{unit}/featureCounts.txt",
        summary="qc/{sample}-{unit}/featureCounts.txt"
    params:
        annotation=config["ref"]["annotation"],
        extra=feature_counts_extra
    threads:
        4
    log:
        "logs/feature_counts/{sample}-{unit}.txt"
    wrapper:
        "file://" + path.join(workflow.basedir, "wrappers/subread/feature_counts")


rule merge_counts:
    input:
        expand("counts/{unit.sample}-{unit.unit}/featureCounts.txt", unit=units.itertuples())
    output:
        "counts/all.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule normalize_counts:
    input:
        "counts/all.tsv"
    output:
        "counts/all.log2.tsv"
    run:
        counts = pd.read_csv(input[0], sep="\t", index_col=list(range(6)))
        norm_counts = np.log2(normalize_counts(counts) + 1)
        norm_counts.to_csv(output[0], sep="\t", index=True)
