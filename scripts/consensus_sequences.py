import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment

alignment_file = sys.argv[1]

alignment = AlignIO.read(alignment_file, "fasta")
alignment_array = np.array(alignment)

distance_matrix = (alignment_array != alignment_array[0,:])
distance = np.sum(distance_matrix, axis=1, dtype=int) / distance_matrix.shape[-1]


hist, bins = np.histogram(distance, bins=20)
bins = 0.5*(bins[1:] + bins[:-1])

plt.figure()
plt.plot(bins, hist, '.-')
plt.grid()
plt.show()
