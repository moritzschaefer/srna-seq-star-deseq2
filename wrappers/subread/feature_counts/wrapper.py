from snakemake.shell import shell
from pathlib import Path

# Run command.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("featureCounts"
      " {snakemake.params.extra}"
      " -a {snakemake.input.annotation}"
      " -o {snakemake.output.counts}"
      " -T {snakemake.threads}"
      " {snakemake.input.bam} {log}")

# fix column header
shell("sed -i 's/%s/{snakemake.wildcards.sample}/' {snakemake.output.counts}" %
      (snakemake.input.bam.replace('/', '\/'), ))

# Move summary to expected location.
summary_path = snakemake.output.counts + '.summary'

if summary_path != snakemake.output.summary:
    shell("mv {summary_path} {snakemake.output.summary}")
