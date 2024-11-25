import numpy as np
import os
import subprocess
# from scipy.optimize import curve_fit


class Profess:
    def __init__(self, settings):
        self.pPROFESS = "/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/PROFESS3.0_modified/pPROFESS"
        self.structures = ['fcc', 'hcp', 'bcc', 'sc']
        self.a = [3.97010529048443894240, 2.80853464541268804666, 3.18038082103364594388, 2.65887848246827296350]
        self.covera = [1, 1.6409235512, 1, 1]
        self.cell = [np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
                     np.array([[1, 0, 0], [-0.5, np.sqrt(3)/2, 0], [0, 0, 1]]),
                     np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
                     np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])]
        self.position = [np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0], [2/3, 1/3, 0.5]]),
                         np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0]])]
        # self.natom = [len(each) for each in self.position]
        # self.ecut = 800
        # self.pseudo = 'al_HC.lda.recpot'
        # self.kernelfile = 'PROFESS_KERNEL.dat'
        self.N = 5
        self.create_files(settings)

    def create_inpt(self, inptfile, ecut, kernelfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WT\nKERN {1}\nMAXI DEN 50\n'.format(
                ecut, kernelfile))

    def create_ion(self, name, ionfile, cell, position, pseudo):
        # create .ion file for PROFESS
        with open(ionfile, 'w') as f:
            f.write('%BLOCK LATTICE_CART\n')
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('%END BLOCK LATTICE_CART\n')
            f.write('%BLOCK POSITIONS_FRAC\n')
            for eachline in position:
                f.write(name)
                f.writelines('    ' + '%.12f' % each for each in eachline)
                f.write('\n')
            f.write('%END BLOCK POSITIONS_FRAC\n')
            f.write(
                '%BLOCK SPECIES_POT\nAl {0}\n%END BLOCK SPECIES_POT\n'.format(pseudo))

    def create_files(self, settings):
        coef = np.linspace(0.9, 1.1, self.N)
        with open('jobs', 'w') as job:
            for structure in settings.structures[-3:]:
                job.write(self.pPROFESS + ' {0}/{1}\n'.format(structure, settings.name))
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i] * self.a[i] * coef[j]
                    cell[:, -1] *= self.covera[i]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_inpt(path+'/'+settings.name+'.inpt', settings.ecut, settings.kernelfile)
                    self.create_ion(settings.name, path+'/'+settings.name+'.ion', cell, self.position[i], settings.pseudo)
                    job.write(self.pPROFESS + ' ' + path + '/' + settings.name + '\n')

    def cal_e_v(self, name, structures=[]):
        # calculate e_v curve and return energy list
        N = self.N * len(self.structures) + len(structures)
        #e_list = [np.zeros(N), np.zeros(N), np.zeros(N)]
        e_list = np.zeros((2, N))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for each in structures[:3]:
            cal = subprocess.Popen(args='mpirun -n 20 {0} {1}/{2}'.format(self.pPROFESS, each, name), shell=True, universal_newlines=False)
            cal.wait()
        for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}{1}/{2}.out'.format(self.structures[i], j, self.id)
                try:
                    get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    total_e = float(get_total_e.stdout.read()) # maybe wrong
                    # get_kinetic_e = subprocess.Popen(args="grep 'TOTAL KINETIC ENERGY' %s|tr -s ' '|cut -d ' ' -f6" % outfile,
                    #                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    # kinetic_e = float(get_kinetic_e.stdout.read())
                    get_hartree_e = subprocess.Popen(args="grep 'Coulombic Energy' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    hartree_e = float(get_hartree_e.stdout.read())
                except:
                    total_e = 1e5
                    # kinetic_e = 1e5
                    hartree_e = 1e5
                e_list[0,self.N * i + j] = total_e
                # e_list[1][i] = kinetic_e
                e_list[1,self.N * i + j] = hartree_e
        for i, stru in enumerate(structures):
            outfile = './{0}/{1}.out'.format(stru, name)
            try:
                get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read()) # maybe wrong
                # get_kinetic_e = subprocess.Popen(args="grep 'TOTAL KINETIC ENERGY' %s|tr -s ' '|cut -d ' ' -f6" % outfile,
                #                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                # kinetic_e = float(get_kinetic_e.stdout.read())
                get_hartree_e = subprocess.Popen(args="grep 'Coulombic Energy' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                hartree_e = float(get_hartree_e.stdout.read())
            except:
                total_e = 1e5
                # kinetic_e = 1e5
                hartree_e = 1e5
            e_list[0,self.N * len(self.structures) + i] = total_e
            # e_list[1][i] = kinetic_e
            e_list[1,self.N * len(self.structures) + i] = hartree_e
        return e_list

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
