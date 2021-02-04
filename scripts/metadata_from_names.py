# Dirty script to generate metadata from names

import sys
import pandas as pd
from Bio import SeqIO

def metadata_from_names(sequences):
    """
    Creates a metadata tsv file from the sequences names.
    """
    columns = ["subtype", "country", "date", "name"]
    raw_names = [seq.name for seq in sequences]
    df = pd.DataFrame(data=None, index=None, columns=columns)
    for raw_name in raw_names:
        tmp_col = raw_name.split(".")
        tmp_col = tmp_col[:3]
        tmp_col.append(raw_name)
        tmp = pd.DataFrame(data=[tmp_col], columns=columns)
        df = df.append(tmp)
    return df

sequences_file = sys.argv[1]
output = sys.argv[2]

sequences = list(SeqIO.parse(sequences_file, "fasta"))
metadata = metadata_from_names(sequences)
metadata.to_csv(output, index=False, sep="\t")
