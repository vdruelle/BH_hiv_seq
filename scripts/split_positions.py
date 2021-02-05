# Quick script to subsample an alignment into 1st 2nd and 3rd positions

import sys
import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqRecord import SeqRecord

alignment_file = sys.argv[1]

alignment = AlignIO.read(alignment_file, "fasta")
alignment_array = np.array(alignment)

basename = alignment_file.replace(".fasta", "")

for position, name in zip([0, 1, 2], ["1st", "2nd", "3rd"]):
    sub_alignment = alignment_array[:, position::3]
    seq_list = []
    for ii in range(alignment_array.shape[0]):
        seq = "".join(sub_alignment[ii, :])
        seq_list += [SeqRecord(Seq(seq, SingleLetterAlphabet()), id=alignment[ii].id,
                               name=alignment[ii].name, description="")]

    sub_alignment = MultipleSeqAlignment(seq_list)
    AlignIO.write([sub_alignment], basename + "_" + name + ".fasta", "fasta")
