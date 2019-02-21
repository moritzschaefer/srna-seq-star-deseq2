import pandas as pd

frames = (pd.read_csv(fp, sep="\t", skiprows=1, index_col=list(range(6)))
          for fp in snakemake.input)
merged = pd.concat(frames, axis=1)

merged.to_csv(snakemake.output[0], sep="\t", index=True)
