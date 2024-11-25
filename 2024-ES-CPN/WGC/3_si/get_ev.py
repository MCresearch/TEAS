import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit


class Profess:
    def __init__(self):
        self.pPROFESS = "PROFESS"
        self.structures = ['cd','btin','hd', 'cbcc','fcc', 'bcc', 'hcp', 'sc']
        self.a = [5.408543163224102, 4.752412941002936, 3.7999674195474435, 6.57, 3.8587830799369325, 3.0802588345439395,2.729216814764073, 2.4922597732684624]
        self.covera = [1, 2.592632187455004/4.752412941002936, 1.656, 1, 1, 1, 1.59270728861, 1]
        self.cell = [np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
                     np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
                     np.array([[1, 0, 0], [-0.5, np.sqrt(3)/2, 0], [0, 0, 1]]),
                     np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
                     np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
                     np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
                     np.array([[1, 0, 0], [-0.5, np.sqrt(3)/2, 0], [0, 0, 1]]),
                     np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])]
        self.position = [np.array([[0, 0, 0], [1/4, 1/4, 1/4]]),
                         np.array([[0, 0, 0], [3/4, 1/4, 1/2]]),
                         np.array([[0.3333333333333333, 0.6666666666666667, 0.0631],
                                [0.6666666666666667, 0.3333333333333333, 0.5631],
                                [0.3333333333333333, 0.6666666666666667, 1 - 0.5631],
                                [0.6666666666666667, 0.3333333333333333, 1 - 0.0631]]),
                         np.array([
                                [0.1484430000000001, 0.8515569999999999, 0.3515569999999999],
                                [0.8515569999999999, 0.6484430000000001, 0.3515569999999999],
                                [0.3515569999999999, 0.8515569999999999, 0.6484430000000001],
                                [0.6484430000000001, 0.6484430000000001, 0.6484430000000001],
                                [0.8515569999999999, 0.1484430000000001, 0.6484430000000001],
                                [0.6484430000000001, 0.8515569999999999, 0.1484430000000001],
                                [0.6484430000000001, 0.1484430000000001, 0.3515569999999999],
                                [0.3515569999999999, 0.3515569999999999, 0.3515569999999999],
                                [0.6484430000000001, 0.3515569999999999, 0.8515569999999999],
                                [0.3515569999999999, 0.1484430000000001, 0.8515569999999999],
                                [0.8515569999999999, 0.3515569999999999, 0.1484430000000001],
                                [0.1484430000000001, 0.1484430000000001, 0.1484430000000001],
                                [0.3515569999999999, 0.6484430000000001, 0.1484430000000001],
                                [0.1484430000000001, 0.3515569999999999, 0.6484430000000001],
                                [0.1484430000000001, 0.6484430000000001, 0.8515569999999999],
                                [0.8515569999999999, 0.8515569999999999, 0.8515569999999999]
                         ]),
                         np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0], [2/3, 1/3, 1/2]]),
                         np.array([[0, 0, 0]])]
        self.natom = [len(each) for each in self.position]
        self.ecut = 3200
        self.pseudo = 'si.gga.recpot'
        self.kernelfile = 'PROFESS_KERNEL.dat'
        self.id = 'Si'
        self.N = 11
        self.create_files()
        self.cal_e_v()

    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WGC\nPARA GAMM 4.2\nEXCH   PBE\n'.format(
                self.ecut))

    def create_ion(self, ionfile, cell, position):
        # create .ion file for PROFESS
        with open(ionfile, 'w') as f:
            f.write('%BLOCK LATTICE_CART\n')
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('%END BLOCK LATTICE_CART\n')
            f.write('%BLOCK POSITIONS_FRAC\n')
            for eachline in position:
                f.write(self.id)
                f.writelines('    ' + '%.12f' % each for each in eachline)
                f.write('\n')
            f.write('%END BLOCK POSITIONS_FRAC\n')
            f.write(
                '%BLOCK SPECIES_POT\nSi {0}\n%END BLOCK SPECIES_POT\n'.format(self.pseudo))

    def create_files(self):
        coef = np.linspace(0.9, 1.1, self.N)
        with open('jobs', 'w') as job:
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i] * self.a[i] * coef[j]
                    cell[:, -1] *= self.covera[i]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_inpt(path+'/'+self.id+'.inpt')
                    self.create_ion(path+'/'+self.id+'.ion', cell, self.position[i])
                    job.write(self.pPROFESS + ' ' + path + '/' + self.id + ' > '+ path + '/log\n')

    def cal_e_v(self):
        # calculate e_v curve and return energy list
        e_list = np.zeros(self.N * len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}{1}/{2}.out'.format(self.structures[i], j, self.id)
                try:
                    get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    total_e = float(get_total_e.stdout.read())/self.natom[i]  # maybe wrong
                    # get_kinetic_e = subprocess.Popen(args="grep 'TOTAL KINETIC ENERGY' %s|tr -s ' '|cut -d ' ' -f6" % outfile,
                    #                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    # kinetic_e = float(get_kinetic_e.stdout.read())
                    # get_hartree_e = subprocess.Popen(args="grep 'Coulombic Energy' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                    #                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    # hartree_e = float(get_hartree_e.stdout.read())
                except:
                    total_e = 1e5
                    # kinetic_e = 1e5
                    # hartree_e = 1e5
                e_list[self.N * i + j] = total_e
                # e_list[1][i] = kinetic_e
                # e_list[2][i] = hartree_e
        v_list = np.zeros_like(e_list)
        coef = np.linspace(0.9, 1.1, self.N)
        for i in range(len(self.structures)):
            v = abs(np.dot(self.cell[i][0], np.cross(self.cell[i][1], self.cell[i][2])))
            v *= self.a[i] ** 3 * self.covera[i] / self.natom[i]
            for j in range(self.N):
                v_list[self.N * i + j] = v * coef[j] ** 3
        with open(self.id + '.curve', 'w') as f:
            f.write('Volume per atom(Ang^3)\tEnergy per atom(eV)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(self.N * len(self.structures)):
                f.write('%.12e\t%.12e\n' % (v_list[i], e_list[i]))

if __name__ == '__main__':
    profess = Profess()
    # def dichotomy(self, f, a, b, eta=1e-7):
    #     # Solve the equation by dichotomy
    #     if f(a) * f(b) > 0:
    #         raise ValueError
    #     if abs(f(a)) <= eta:
    #         return a
    #     if abs(f(b)) <= eta:
    #         return b
    #     c = (a+b)/2
    #     while abs(f(c)) > eta:
    #         if f(c) * f(a) > 0:
    #             a = c
    #         else:
    #             b = c
    #         c = (a+b)/2
    #     return c

    # def fit(self, volume, energy, v1, v2, v0=0):
    #     #volume in A^3, energy in eV
    #     def _energy(v, a, b, c, d, e, f):
    #         return a + b * v + c * v ** 2 + d * v ** 3 + e * v ** 4 + f * v ** 5

    #     def order1(v):
    #         return b + 2 * c * v + 3 * d * v ** 2 + 4 * e * v ** 3 + 5 * f * v ** 4
    #     try:
    #         a, b, c, d, e, f = curve_fit(_energy, volume, energy)[0]
    #     except:
    #         return 1e5, 1e5, 1e5
    #     try:
    #         v0 = self.dichotomy(order1, v1, v2)
    #     except:
    #         v0 = v0
    #     e0 = _energy(v0, a, b, c, d, e, f)
    #     B = (2 * c + 6 * d * v0 + 12 * e * v0 ** 2 + 20 * f * v0 ** 3) * \
    #         v0 * 160.2  # 160.2 for converting unit to GPa
    #     return e0, v0, B
