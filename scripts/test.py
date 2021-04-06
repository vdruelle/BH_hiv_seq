import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

filename = sys.argv[1]

df = pd.read_csv(filename, sep='\t')
dates = list(df["date"])
dates = [float(date[:4]) for date in dates]


plt.figure()
plt.hist(dates, bins=np.arange(1970, 2020))
# plt.yscale("log")
plt.xlabel("Year")
plt.ylabel("# sequences")
plt.show()
