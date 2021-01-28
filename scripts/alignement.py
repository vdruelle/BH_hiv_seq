from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import SingleLetterAlphabet
import numpy as np


DATA_FOLDER = "data/"
ALIGNMENT_FOLDER = DATA_FOLDER + "alignments/"


def remove_gaps(alignment, ref_row=0):
    """
    Removes the column where the reference has gaps and return the obtained alignment.
    """

    alignment_array = np.array(alignment)
    idxs = np.where(alignment_array[ref_row, :] != "-")[0]
    alignment_array = alignment_array[:, idxs]

    seq_list = []
    for ii in range(alignment_array.shape[0]):
        seq = "".join(alignment_array[ii, :])
        seq_list += [SeqRecord(Seq(seq, SingleLetterAlphabet()), id=alignment[ii].id, name=alignment[ii].name, description="")]

    alignment = MultipleSeqAlignment(seq_list)
    return alignment

def get_pol_region(alignment):
    """
    Slices the alignment to get the pol region only. Uses HXB2 coordinatese for the selection.
    """
    HXB2_start = 2084
    HXB2_stop = 5095
    return alignment[:, HXB2_start:HXB2_stop+1]

if __name__ == '__main__':
    alignment = AlignIO.read(ALIGNMENT_FOLDER + "raw/mafft_alignment.fasta", 'fasta')
    alignment = remove_gaps(alignment, ref_row=0)
    alignment = get_pol_region(alignment)
    AlignIO.write([alignment], ALIGNMENT_FOLDER + "to_HXB2/pol.fasta","fasta")
