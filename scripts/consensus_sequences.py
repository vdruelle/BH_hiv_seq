import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

alignment_file = sys.argv[1]
output_file = sys.argv[2]

alignment = AlignIO.read(alignment_file, "fasta")
alignment_array = np.array(alignment)

# Consensus computation
consensus_sequence = []
for ii in range(alignment_array.shape[-1]):
    values, counts = np.unique(alignment_array[:,ii], return_counts=True)
    consensus_nucleotide = values[np.argmax(counts)]
    consensus_sequence = consensus_sequence + [consensus_nucleotide]

consensus_sequence = np.array(consensus_sequence)
consensus_sequence = Seq("".join(consensus_sequence))
consensus_sequence = SeqRecord(seq=consensus_sequence, id=f"Consensus_{alignment_file}", name="", description="")

with open(output_file, "w") as handle:
    SeqIO.write([consensus_sequence], handle, "fasta")

# Distance to consensus sequence
# distance_matrix = (alignment_array != consensus_sequence)
# distance = np.sum(distance_matrix, axis=1, dtype=int) / distance_matrix.shape[-1]
