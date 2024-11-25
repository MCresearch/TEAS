import numpy as np
from create_kernel import *
from settings import *

g = np.linspace(0, ng * dg, ng + 1)
r = np.linspace(dr, nr * dr, nr)

jle = get_real_bessel(dr, nr, rcut, n_bessel, l)
jG = get_recip_bessel(r, dr, g, jle, n_bessel)
jle_inte = np.array([cotes(r, each * r * r) for each in jle])
J = np.array([jG[i, 0] - jle_inte[i] / jle_inte[-1] * jG[-1, 0]
              for i in range(n_bessel - 1)])

coef = np.loadtxt("./RCUT8_8Jlq_GCUT10/DUMP_COEF.dat", skiprows=3).T[1]

print(np.dot(coef, jle_inte))#0.127314021207093
print(np.dot(coef, jG[:,0]))#1.5999288060457284