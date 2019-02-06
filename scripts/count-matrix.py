import pandas as pd
from os import path

frames = (pd.read_csv(fp, sep="\t", skiprows=1, index_col=list(range(6)))
          for fp in snakemake.input)
merged = pd.concat(frames, axis=1)

# Extract sample names.
merged = merged.rename(columns=lambda c: path.splitext(path.basename(c))[0])

merged.to_csv(snakemake.output[0], sep="\t", index=True)
