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


refine_file = sys.argv[1]
all_file = sys.argv[2]
first_file = sys.argv[3]
second_file = sys.argv[4]
third_file = sys.argv[5]
output_file = sys.argv[6]


rate_dict = {}
rate_dict["all"] = get_mutation_rate(all_file, refine_file)
rate_dict["first"] = get_mutation_rate(first_file, refine_file)
rate_dict["second"] = get_mutation_rate(second_file, refine_file)
rate_dict["third"] = get_mutation_rate(third_file, refine_file)

with open(output_file, "w") as output:
    json.dump(rate_dict, output)
