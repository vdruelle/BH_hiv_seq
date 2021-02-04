# Dirty script to add HXB2

import sys
from Bio import SeqIO

def insert_sequence(record, sequence_file):
    "Insert the sequence sequence at the beginning of the file."
    record.insert(0, SeqIO.read(sequence_file, "fasta"))
    return record

sequences_file = sys.argv[1]
reference_file = sys.argv[2]

sequences = list(SeqIO.parse(sequences_file, "fasta"))
sequences = insert_sequence(sequences, reference_file)
SeqIO.write(sequences, sequences_file, "fasta")
