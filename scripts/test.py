import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

sequences = sys.argv[1]
metadata = sys.argv[2]
output = sys.argv[3]

df = pd.read_csv(metadata, sep='\t')
tmp = df["date"]
tmp = []

df.sort_values(by=["date"])
print(df["date"])
