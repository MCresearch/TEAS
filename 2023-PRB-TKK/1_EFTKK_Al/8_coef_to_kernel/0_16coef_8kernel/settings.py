import numpy as np


class Settings:
    def __init__(self):
        self.dr = 0.01
        self.rcut = 8
        self.n_bessel = 8
        self.l = 0
        self.nr = int(self.rcut/self.dr)
        if self.nr % 2 == 0:
            self.nr += 1
        self.dg = 0.0001
        self.gcut = 20  # gcut=10时跑并行会报错
        self.ng = int(self.gcut/self.dg)
        if self.ng % 2 == 0:
            self.ng += 1
        
        self.g = np.linspace(0, self.ng * self.dg, self.ng + 1)
        self.r = np.linspace(self.dr, self.nr * self.dr, self.nr)

        self.coef = np.loadtxt("4optcoef.dat")
        # self.coef[0] = 0.1
        # self.coef[1] = 0.1
        # self.coef[2] = 0.1
        # self.coef[3] = 0.1
        # self.coef[4] = 0.1
        self.kinetic = 4

        # self.a = 3.968
        # self.a_list = np.linspace(0.9, 1.1, 11) * self.a
        # self.v_list = self.a_list ** 3 / 4
        self.structures = ['fcc100', 'fcc110', 'fcc111', 'vacancy111', 'vacancy211', 'vacancy221']
        self.ecut = 800
        self.pseudo = 'al_HC.lda.recpot'
        self.kernelfile = 'PROFESS_KERNEL.dat'

        # self.lowest = 1e9
        self.T = 1
        self.cooling_rate = 0.5
        self.nT = 10
        self.nsteps = 300
        self.rate = np.ones(self.n_bessel-1)*0.1
        self.bias = 1e-3
        self.ncheck = 100
        self.accept_high = 0.4
        self.accept_low = 0.2

        self.Ha2eV = 27.2114
        self.weight = np.array([1/20, 3/20, 1/10, 1/20])

        self.name = 'al'
        self.targetfile = 'al_target'
