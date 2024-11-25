import numpy as np

MgAl = np.zeros(8)
Al = np.zeros(8)
Mg = np.zeros(8)

i = 0
with open('./0_origin_data.txt') as f:
    f.readline()
    while i < 8:
        x = f.readline()
        temp = np.array(x.split()[1:])
        MgAl[i] = temp[0]
        Al[i] = temp[1]
        Mg[i] = temp[2]
        i += 1
print((4*MgAl - (3 * Mg + Al))/4 * 1000) 
