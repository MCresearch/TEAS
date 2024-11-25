import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit


class PROFESS:
    def __init__(self):
        self.pPROFESS = "pPROFESS"
        # self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.structures = ['bcc', 'fcc', 'sc', 'diamond']
        self.a = [3.431257319966733, 4.327318413514269, 2.737252293704602, 5.90960801642828]
        self.covera = [1, 1, 1, 1]
        self.cell = [np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
                     np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
                     np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
                     np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])]
        self.position = [np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0], [1/4, 1/4, 1/4]])]
        self.natom = [len(each) for each in self.position]
        self.ecut = 900
        self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/reci/'
        self.pseudo = 'li.lda.recpot'
        self.kernel_dir = "/home/dell/1_work/TKK_abacus/2_Li/1_kernels/"
        # self.kernelfiles = []
        self.kernelfiles = ['PROFESS_KERNEL.dat','1PROFESS_KERNEL.dat','3PROFESS_KERNEL.dat','4PROFESS_KERNEL.dat','5PROFESS_KERNEL.dat','8PROFESS_KERNEL.dat','11PROFESS_KERNEL.dat']
        self.id = 'Li'
        for kernelfile in self.kernelfiles:
            os.chdir("/home/dell/1_work/TKK_abacus/2_Li/2_bulk_properties/1_optimize_cell/")
            os.mkdir(kernelfile[:-4])
            os.chdir('./'+kernelfile[:-4])
            self.kernelfile = kernelfile
            os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
            os.system('cp {0} .'.format(self.pseudo_dir+self.pseudo))
            self.create_files()
            self.calculate()

    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('MINI CEL\nMETH ION NON\nECUT	{0}\nKINE	WT\nKERN {1}\nMAXI DEN 50\n'.format(
                self.ecut, self.kernelfile))

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
                '%BLOCK SPECIES_POT\n{0} {1}\n%END BLOCK SPECIES_POT\n'.format(self.id, self.pseudo))

    def create_files(self):
        with open('jobs', 'w') as job:
            for i in range(len(self.structures)):
                cell = self.cell[i] * self.a[i]
                cell[:, -1] *= self.covera[i]
                path = self.structures[i]
                os.mkdir(path)
                self.create_inpt(path+'/'+self.id+'.inpt')
                self.create_ion(path+'/'+self.id+'.ion', cell, self.position[i])
                job.write(self.pPROFESS + ' ' + path + '/' + self.id + '\n')

    def calculate(self):
        # perform structure optimization
        # e_list = np.zeros((2,self.N * len(self.structures)))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        # for i in range(len(self.structures)):
        #     for j in range(self.N):
        #         outfile = './{0}{1}/OUT.blpstest/running_scf.log'.format(self.structures[i], j)
#                 try:
#                     get_total_e = subprocess.Popen(args="grep 'final etot' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
#                                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#                     total_e = float(get_total_e.stdout.read())  # maybe wrong
# #                    get_kinetic_e = subprocess.Popen(args="grep 'E_one_elec' %s|tr -s ' '|cut -d ' ' -f4" % outfile,
#  #                                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#   #                  kinetic_e = float(get_kinetic_e.stdout.read())
#                     get_hartree_e = subprocess.Popen(args="grep 'E_Hartree' %s|tr -s ' '|cut -d ' ' -f4" % outfile,
#                                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
#                     hartree_e = float(get_hartree_e.stdout.read())
#                 except:
#                     total_e = 1e5
#    #                 kinetic_e = 1e5
#                     hartree_e = 1e5
#                 e_list[0,self.N * i + j] = total_e
# #                e_list[1,self.N * i + j] = kinetic_e
#                 e_list[1,self.N * i + j] = hartree_e
#         # v_list = np.zeros_like(e_list)
#         # coef = np.linspace(0.9, 1.1, self.N)
#         # for i in range(len(self.structures)):
#         #     v = abs(np.dot(self.cell[i][0], np.cross(self.cell[i][1], self.cell[i][2])))
#         #     v *= self.a[i] ** 3 * self.covera[i] / self.natom[i]
#         #     for j in range(self.N):
#         #         v_list[self.N * i + j] = v * coef[j] ** 3
#         with open(self.id + '_target', 'w') as f:
#             f.write('1-5:fcc al 6-10:hcp 11-15:bcc 15-20:sc\n\
# Energy(eV)                   Hartree             natoms\n')
#             for i in range(self.N * len(self.structures)):
#                 f.write('%.12e\t%.12e\t%d\n' % (e_list[0,i], e_list[1,i], self.natom[i//5]))

if __name__ == '__main__':
    profess = PROFESS()
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
