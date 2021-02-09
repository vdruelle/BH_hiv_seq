import sys
import numpy as np
import json

def get_mutation_rate(gtr_file, refine_file):
    """
    Returns the mutation rate in mutation per site per year from the GTR model and its corresponding
    tree file from augur refine output.
    """
    with open(gtr_file, "r") as file:
        gtr = json.load(file)
    with open(refine_file, "r") as file:
        refine = json.load(file)

    mutation_rate = gtr["mu"] * refine["clock"]["rate"]

    return mutation_rate

nb_nucleotide = 3012
gtr_file = "gtr/pol.json"
refine_file = "intermediate_files/branch_lengths_pol.json"
mut_rate = get_mutation_rate(gtr_file, refine_file)
