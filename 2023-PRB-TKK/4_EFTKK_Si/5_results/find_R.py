import numpy as np

prefix = [1, 4, 7, 10, 13]
logfiles = [str(_) + 'log' for _ in prefix]
for file in logfiles:
    log = np.loadtxt(file, skiprows=3)
    print(min(log[:,-1]))