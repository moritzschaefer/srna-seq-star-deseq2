import pandas as pd

counts = pd.read_csv(snakemake.input[0], sep="\t", index_col=list(range(6)))
norm_counts = 1e6 * counts / counts.sum(axis=0)
norm_counts.to_csv(snakemake.output[0], sep="\t", index=True)
