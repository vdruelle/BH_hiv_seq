import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO
import json


def get_consensus_sequence(filename):
    """
    Loads the consensus sequence from the given file.
    """
    if ".json" in filename:
        with open(filename) as f:
            data = json.load(f)
            consensus_sequence = list(data["nodes"]["NODE_0000000"]["sequence"])
        consensus_sequence = np.array(consensus_sequence)
    else:
        consensus = AlignIO.read(filename, "fasta")
        consensus_sequence = np.array(consensus)[0]
    return consensus_sequence


def get_gap_mask(alignment_array, threshold=0.1):
    "Return a vector were true are the sites seen with less than threshold fraction of N."
    gaps = alignment_array == "N"
    gap_proportion = np.sum(gaps, axis=0, dtype=int) / gaps.shape[0]
    return gap_proportion < threshold


def get_mean_distance_in_time(alignment_file, consensus_sequence, subtype="B"):
    alignment = AlignIO.read(alignment_file, "fasta")
    names = [seq.id for seq in alignment]
    dates = np.array([int(name.split(".")[2]) for name in names])
    subtypes = np.array([name.split(".")[0] for name in names])
    alignment_array = np.array(alignment)
    gap_mask = get_gap_mask(alignment_array)

    # Selecting subtype
    alignment_array = alignment_array[subtypes == subtype]
    dates = dates[subtypes == subtype]

    # Distance to consensus sequence
    distance_matrix = (alignment_array != consensus_sequence)[:, gap_mask]
    distance = np.sum(distance_matrix, axis=1, dtype=int) / distance_matrix.shape[-1]

    # Distance average per year
    average_distance = []
    std_distance = []
    years = np.unique(dates)
    for year in years:
        average_distance += [np.mean(distance[dates == year])]
        std_distance += [np.std(distance[dates == year])]
        average_distance_in_time = np.array(average_distance)
        std_distance_in_time = np.array(std_distance)

    return years, average_distance_in_time, std_distance_in_time


alignment_file = "data/alignments/to_HXB2/pol_1000.fasta"
# consensus_file = "data/alignments/to_HXB2/pol_1000_consensus.fasta"
consensus_file = "intermediate_files/pol_1000_nt_muts.json"

plt.figure()
subtypes = ["B", "C"]
colors = ["C0", "C1"]
c = 0
for subtype in subtypes:
    consensus_sequence = get_consensus_sequence(consensus_file)
    years, dist, std = get_mean_distance_in_time(alignment_file, consensus_sequence, subtype)
    fit = np.polyfit(years[std != 0], dist[std != 0], deg=1, w=(1 / std[std != 0]))
    plt.errorbar(years, dist, yerr=std, fmt=".", label=subtype, color=colors[c])
    plt.plot(years, np.polyval(fit, years), "--",
             color=colors[c], label=f"{round(fit[0],5)}x + {round(fit[1],5)}")
    c += 1

plt.grid()
plt.xlabel("Time [years]")
plt.ylabel("Average fraction difference")
plt.legend()
plt.savefig("20042021-Distance_subtypes.png", format="png")
plt.show()
