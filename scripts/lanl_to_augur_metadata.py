# Dirty script to transform lanl metadata to augur format

import sys
import pandas as pd

lanl = sys.argv[1]

df = pd.read_csv(lanl, sep="\t", skiprows=2)
df = df[["Accession", "Organism", "Sampling Year"]]
df = df.rename(columns={"Accession":"strain", "Organism":"virus", "Sampling Year":"date"})
df.to_csv(lanl[:-4]+"_augur.tsv", sep="\t", index=False)
