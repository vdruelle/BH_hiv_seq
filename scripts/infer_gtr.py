# Quick script to infer a gtr model using treetime

import sys
import json
import numpy as np
from Bio import AlignIO
from treetime import TreeAnc

def output_JSON(gtr_model, output_file):
    gtr = gtr_model.__dict__
    save_dict = {}
    for key in ["_mu", "_Pi", "_W"]:
        if isinstance(gtr[key], np.ndarray):
            save_dict[key.replace("_", "")] = gtr[key].tolist()
        else:
            save_dict[key.replace("_", "")] = gtr[key]

    with open(output_file, 'w') as output:
        json.dump(save_dict, output)


alignment_file = sys.argv[1]
tree_file = sys.argv[2]
output_file = sys.argv[3]

alignment = AlignIO.read(alignment_file, "fasta")
tt = TreeAnc(tree=tree_file, aln=alignment_file)
gtr = tt.infer_gtr(marginal=True, normalized_rate=False)
output_JSON(gtr, output_file)
