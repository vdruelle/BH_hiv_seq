from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import SingleLetterAlphabet
import numpy as np
import pandas as pd
import os
import pickle
from treetime import TreeTime
from treetime.utils import parse_dates

########## File and directory location ##########

DATA_FOLDER = "data/"
INTERMEDIATE_FOLDER = "intermediate_files/"
ALIGNMENT_FOLDER = DATA_FOLDER + "alignments/"
RAW_DATA = DATA_FOLDER + "raw/pol.fasta"
RAW_SUB_SAMPLED = DATA_FOLDER + "raw/pol_subsampled.fasta"
REFERENCE_HXB2 = DATA_FOLDER + "reference/HXB2.fasta"
RAW_ALIGNMENT = ALIGNMENT_FOLDER + "raw/pol.fasta"
HXB2_ALIGNMENT = ALIGNMENT_FOLDER + "to_HXB2/pol.fasta"
HXB2_ALIGNMENT_META = ALIGNMENT_FOLDER + "to_HXB2/pol_metadata.csv"
HXB2_ALIGNMENT_TREE = ALIGNMENT_FOLDER + "to_HXB2/pol_tree.nwk"

########## What I used to create / filter the data ##########


def sub_sample_raw_data():
    """
    Subsample the raw fasta file for easier analysis. Saves the sub_sampled file in the same directory.
    Adds HXB2 sequence as the first sequence.
    """
    nb_sample = 1000

    # This is not random so it contains mainly sequences from the same time
    # record =  list(SeqIO.parse(DATA_FOLDER + "raw/pol.fasta", "fasta"))
    # sub_sampled = record[:nb_sample]

    # This is random
    os.system(f"seqtk sample -s100 {RAW_DATA} {nb_sample} > {RAW_SUB_SAMPLED}")
    sub_sampled = list(SeqIO.parse(RAW_SUB_SAMPLED, "fasta"))
    sub_sampled = insert_sequence(sub_sampled, REFERENCE_HXB2)
    SeqIO.write(sub_sampled, RAW_SUB_SAMPLED, "fasta")


def align_subsampled_data():
    """
    Uses mafft to align the sequences in the subsampled data.
    """
    os.system(f"mafft {RAW_SUB_SAMPLED} > {RAW_ALIGNMENT}")


def MSA_pol_HXB2():
    """
    Uses the raw subsampled sequences, align them to HXB2 and remove regions that correspond to gap in HXB2.
    Saves the newly obtained MultiSequenceAlignement in fasta.
    """
    alignment = AlignIO.read(RAW_ALIGNMENT, 'fasta')
    alignment = remove_gaps(alignment, ref_row=0)
    alignment = get_pol_region(alignment)
    AlignIO.write([alignment], HXB2_ALIGNMENT, "fasta")


def make_metadata():
    """
    Creates the metadata file from the names in the alignment. Saves it to the same folder as alignment.
    """
    alignment = AlignIO.read(HXB2_ALIGNMENT, "fasta")
    df = metadata_from_names(alignment)
    df.to_csv(HXB2_ALIGNMENT_META, index=False)


def make_FastTree():
    """
    Uses Treetime on created the alignment.
    """
    os.system(f"fasttree -nt {HXB2_ALIGNMENT} > {HXB2_ALIGNMENT_TREE}")


def make_TreeTime():
    """
    Runs treetime and saves the results.
    """
    dates = parse_dates(HXB2_ALIGNMENT_META)
    tree = TreeTime(gtr="Jukes-Cantor", tree=HXB2_ALIGNMENT_TREE,
                    precision=1, aln=HXB2_ALIGNMENT, verbose=2, dates=dates)
    result = tree.run(root='best', infer_gtr=True, relaxed_clock=False, max_iter=2,
                      branch_length_mode='input', n_iqd=3, resolve_polytomies=True,
                      Tc='skyline', time_marginal="assign", vary_rate=True)

    assert result, "Error while running the tree."

    return tree


########## Helper functions #########
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
    columns = ["subtype", "country", "date", "name"]
    raw_names = [seq.name for seq in alignment]
    df = pd.DataFrame(data=None, index=None, columns=columns)
    for raw_name in raw_names:
        tmp_col = raw_name.split(".")
        tmp_col = tmp_col[:3]
        tmp_col.append(raw_name)
        tmp = pd.DataFrame(data=[tmp_col], columns=columns)
        df = df.append(tmp)
    return df


def runTree():
    """
    Run analysis from console.
    """
    os.system(f"treetime --aln {HXB2_ALIGNMENT} --tree {HXB2_ALIGNMENT_TREE} --dates {HXB2_ALIGNMENT_META} --outdir {INTERMEDIATE_FOLDER}treetime")
    os.system(f"treetime clock --aln {HXB2_ALIGNMENT} --tree {HXB2_ALIGNMENT_TREE} --dates {HXB2_ALIGNMENT_META} --reroot least-squares --outdir {INTERMEDIATE_FOLDER}treetime_clock")
    os.system(f"treetime ancestral --aln {HXB2_ALIGNMENT} --tree {HXB2_ALIGNMENT_TREE} --outdir {INTERMEDIATE_FOLDER}treetime_ancestral")

if __name__ == '__main__':
    # sub_sample_raw_data()
    # align_subsampled_data()
    # MSA_pol_HXB2()
    # make_metadata()
    # make_FastTree()
    # tree = make_TreeTime()
    runTree()
