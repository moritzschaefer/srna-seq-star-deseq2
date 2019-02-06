from snakemake.shell import shell
from pathlib import Path

# Run command.
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("featureCounts"
      " {snakemake.params.extra}"
      " -a {snakemake.params.annotation}"
      " -o {snakemake.output.counts}"
      " -T {snakemake.threads}"
      " -M -f -O -t miRNA -g Name"
      " {snakemake.input.bam} {log}")

# fix column header
shell("sed 's/%s/{snakemake.wildcards.sample}/' {snakemake.output.counts}" %
      (Path(snakemake.input.bam).name, ))

# Move summary to expected location.
summary_path = snakemake.output.counts + '.summary'

if summary_path != snakemake.output.summary:
    shell("mv {summary_path} {snakemake.output.summary}")
