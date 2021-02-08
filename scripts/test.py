import sys
import numpy as np
import json

def get_gtr(file_name):
    "Opens the result of gtr model from treetime and returns the parameters."
    with open(file_name, 'rb') as file:
        gtr = json.load(file)

    mu = gtr["mu"]
    W = np.array(gtr["W"])
    Pi = np.array(gtr["Pi"])
    return mu, W, Pi

nb_nucleotide = 3012
file_name = "gtr/pol.json"
mu, W, Pi = get_gtr(file_name)
mut_rate_nucleotide = np.sum(W, axis=1)
mut_rate_weighted = mut_rate_nucleotide*Pi
average_mut_rate = np.sum(mut_rate_weighted) / nb_nucleotide
print(average_mut_rate / 365)
