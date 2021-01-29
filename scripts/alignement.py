from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import SingleLetterAlphabet
import numpy as np
import pandas as pd


DATA_FOLDER = "data/"
ALIGNMENT_FOLDER = DATA_FOLDER + "alignments/"

# What I used to create / filter the data
def sub_sample_raw_data():
    """
    Subsample the raw fasta file for easier analysis. Saves the sub_sampled file in the same directory.
    Adds HXB2 sequence as the first sequence.
    """
    nb_sample = 1000
    record =  list(SeqIO.parse(DATA_FOLDER + "raw/pol.fasta", "fasta"))
    sub_sampled = record[:nb_sample]
    sumb_sampled = insert_sequence(sub_sampled, DATA_FOLDER + "reference/HXB2.fasta")
    SeqIO.write(sub_sampled, DATA_FOLDER + "raw/pol_subsampled.fasta", "fasta")


def MSA_pol_HXB2():
    """
    Uses the raw subsampled sequences, align them to HXB2 and remove regions that correspond to gap is HXB2.
    Saves the newly obtained MultiSequenceAlignement in fasta.
    """
    alignment = AlignIO.read(ALIGNMENT_FOLDER + "raw/mafft_alignment.fasta", 'fasta')
    alignment = remove_gaps(alignment, ref_row=0)
    alignment = get_pol_region(alignment)
    AlignIO.write([alignment], ALIGNMENT_FOLDER + "to_HXB2/pol.fasta", "fasta")


# Helper functions
def insert_sequence(record, sequence_file):
    "Insert the sequence sequence at the beginning of the file."
    record.insert(0, SeqIO.read(sequence_file, "fasta"))
    return record


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
        seq_list += [SeqRecord(Seq(seq, SingleLetterAlphabet()), id=alignment[ii].id,
                               name=alignment[ii].name, description="")]

    alignment = MultipleSeqAlignment(seq_list)
    return alignment


def get_pol_region(alignment):
    """
    Slices the alignment to get the pol region only. Uses HXB2 coordinatese for the selection.
    """
    HXB2_start = 2084
    HXB2_stop = 5095
    return alignment[:, HXB2_start:HXB2_stop + 1]


def metadata_from_names(alignment):
    """
    Creates a metadata tsv file from the MSA using the names of the sequences.
    """
    columns = ["subtype", "country", "date", "name", "accession"]
    raw_names = [seq.name for seq in alignment]
    df = pd.DataFrame(data=None, index=None, columns=columns)
    for raw_name in raw_names:
        tmp_col = raw_name.split(".")
        if len(tmp_col) > 5:
            tmp_col[4] = "".join(tmp_col[4:])
            tmp_col = tmp_col[:5]
        tmp = pd.DataFrame(data=[tmp_col], columns=columns)
        df = df.append(tmp)
    return df


if __name__ == '__main__':
    sub_sample_raw_data()
