from os import path


def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                      group=[1, 2], **wildcards)
    # single end sample
    return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


rule align:
    input:
        fq1=get_trimmed
    output:
        # see STAR manual for additional output files
        "star/{sample}-{unit}/Aligned.out.bam",
        # "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="{}".format(config["params"]["star"])

        # TODO requires --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon miRNA_primary_transcript

        # extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
        #       config["ref"]["annotation"], config["params"]["star"])
    threads: 10
    wrapper:
        # "file://" + path.join(workflow.basedir, "wrappers/star/align")
        "file:///home/schamori/moritzsphd/lib/srna-seq-star-deseq2/wrappers/star/align"  # TODO get this right (use github!!)
        # "https://github.com/moritzschaefer/srna-seq-star-deseq2/wrappers/star/align"
