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

        self.coef = np.loadtxt("7optcoef.dat")
        # self.coef[0] = 0.1
        # self.coef[1] = 0.1
        # self.coef[2] = 0.1
        self.kinetic = 4

        #C for cd si, B for beta-tin si
        #cd si
        self.aC = 5.40835206773660193136
        self.aC_list = np.linspace(0.9, 1.1, 11) * self.aC
        # self.vC_list = self.aC_list ** 3 / 8
        #beta-tin si
        self.aB = 4.750454639824956 
        self.cB = 2.5915638554921636
        self.aB_list = np.linspace(0.9, 1.1, 11) * self.aB
        self.cB_list = np.linspace(0.9, 1.1, 11) * self.cB
        # self.vB_list = self.aB_list ** 2 * self.cB_list / 4
        
        self.structures = ['cd100', 'vacancy111', 'vacancy211']

        self.ecut = 1000 # eV
        self.pseudo = 'si.lda.recpot'
        self.kernelfile = 'PROFESS_KERNEL.dat'

        # self.lowest = 1e9
        self.T = 1
        self.cooling_rate = 0.5
        self.nT = 8
        self.nsteps = 250
        self.rate = np.ones(self.n_bessel-1)*0.1
        self.bias = 1e-3
        self.ncheck = 100
        self.accept_high = 0.4
        self.accept_low = 0.2

        self.Ha2eV = 27.2114
        self.weight = np.array([1/25, 3/25, 1/10, 1/20])

        self.name = 'si'
        self.targetfile = 'si_target'
