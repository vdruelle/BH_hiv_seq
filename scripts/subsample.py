import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

sequences = sys.argv[1]
metadata = sys.argv[2]
number = sys.argv[3]
output_sequences = sys.argv[4]
output_metadata = sys.argv[5]

df = pd.read_csv(metadata, sep='\t')
df = df[df["subtype"] != "O"]  # removing subtype O
df["date"] = [int(date[:4]) for date in df["date"]]
df = df.sort_values("date")

# Computing number of sequences for each year
hist, bins = np.histogram(df["date"], bins=np.arange(1976, 2022))
hist = hist[:-1]
bins = bins[:-1]
average_per_bin = float(number) / len(bins)
nb_per_bin = np.zeros_like(hist)
nb_per_bin[hist < average_per_bin] = hist[hist < average_per_bin]
nb_per_bin[hist > average_per_bin] = average_per_bin

while np.sum(nb_per_bin) < int(number):
    nb_per_bin[hist >= nb_per_bin + 1] += 1

# adjust to have exactly the desired number of sequences
while np.sum(nb_per_bin) > int(number):
    ii = np.argmax(nb_per_bin)
    nb_per_bin[ii] -= 1


# Getting the sequences names
seq_names = []
for ii, year in enumerate(bins[:-1]):
    if nb_per_bin[ii] == hist[ii]:
        seq_names = seq_names + df[df["date"] == year]["strain"].to_list()
    else:
        permut = np.random.permutation(hist[ii])
        names = df[df["date"] == year]["strain"].to_numpy()
        names = names[permut[:nb_per_bin[ii]]].tolist()
        seq_names = seq_names + names

# Selecting the given sequences
for ii, seq_name in enumerate(seq_names):
    if ii == 0:
        output_df = df[df["strain"] == seq_name]
    else:
        output_df = output_df.append(df[df["strain"] == seq_name])

# Reformating the date field for TreeTime
output_df["date"] = [str(date)+"-XX-XX" for date in output_df["date"]]

# Creating the output sequences
sequences = list(SeqIO.parse(sequences, "fasta"))
sequences = [seq for seq in sequences if seq.name in seq_names]

# Cleaning the sequences (some characters are non ATGC-N sometimes)
for ii in range(len(sequences)):
    seq = np.array(sequences[ii])
    tmp1 = np.logical_and(seq!="a", seq!="t")
    tmp2= np.logical_and(seq!="g", seq!="c")
    tmp3= np.logical_and(seq!="-", seq!="n")
    tmp4 = np.logical_and(tmp1, tmp2)
    tmp5 = np.logical_and(tmp4, tmp3)
    seq[tmp5] = "n"
    sequences[ii].seq = Seq("".join(seq))

# Creating the output files
output_df.to_csv(output_metadata, index=False, sep="\t")
SeqIO.write(sequences, output_sequences, "fasta")
