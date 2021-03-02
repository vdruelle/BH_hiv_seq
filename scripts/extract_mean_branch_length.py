import sys
import numpy as np
import json


def mean_branch_length(refine_file):
    """
    Returns the mutation rate in mutation per site per year from the GTR model and its corresponding
    tree file from augur refine output.
    """
    with open(refine_file, "r") as file:
        refine = json.load(file)

    total_length = 0
    for key in refine["nodes"].keys():
        total_length += refine["nodes"][key]["branch_length"]
    mean_length = total_length/len(refine["nodes"])

    return mean_length


refine_file = sys.argv[1]
output_file = sys.argv[2]

mean_length = {"mean_branch_length": mean_branch_length(refine_file)}

with open(output_file, "w") as output:
    json.dump(mean_length, output)
