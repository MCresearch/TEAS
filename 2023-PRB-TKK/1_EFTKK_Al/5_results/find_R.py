import numpy as np

prefix = [1, 3, 4, 5, 8, 11]
logfiles = [str(_) + 'log' for _ in prefix]
for file in logfiles:
    log = np.loadtxt(file, skiprows=3)
    print(min(log[:,-1]))
