import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit


class Abacus:
    def __init__(self):
        self.abacus = "/home/dell/1_work/5_ABACUS_DEBUG/0_abacus_develop/abacus-develop/build/abacus"
        self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.structures = ['hd', 'cbcc']
        self.a = [3.7999674195474435, 
                #   5.194174936740699, 
                #   5.44,
                  6.57]
        self.covera = [1.656, 1.]
        self.cell = [np.array([[1, 0, 0], [-0.5, np.sqrt(3)/2, 0], [0, 0, 1]]),
                     np.array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])]
        self.position = [np.array([[0.3333333333333333, 0.6666666666666667, 0.0631],
                                   [0.6666666666666667, 0.3333333333333333, 0.5631],
                                   [0.3333333333333333, 0.6666666666666667, 1 - 0.5631],
                                   [0.6666666666666667, 0.3333333333333333, 1 - 0.0631]]),
                        #  np.array([[0, 0, 0],
                        #            [0.5, 0.5, 0.5],
                        #            [0.1003, 0.1003, 0.1003],
                        #            [0.5 + 0.1003, 0.5 + 0.1003, 0.5 + 0.1003],
                        #            [0.5 + 0.1003, 0.5 - 0.1003, 1.0 - 0.1003],
                        #            [1-0.1003, 1-0.1003, 1-0.1003],
                        #            [1-(0.5 + 0.1003), 1-(0.5 - 0.1003), 0.1003],
                        #            [1-0.1003, 0.5 + 0.1003, 0.5 - 0.1003],
                        #            [0.5 - 0.1003, 1-0.1003, 0.5 + 0.1003]]),
                        #  np.array([[0, 0, 0], #refer to mp-149
                        #            [0.5, 0.5, 0],
                        #            [0.5, 0, 0.5],
                        #            [0, 0.5, 0.5],
                        #            [0.1003, 0.1003, 0.1003],
                        #            [0.5 + 0.1003, 0.5 + 0.1003, 0.0 + 0.1003],
                        #            [0.5 + 0.1003, 0.0 + 0.1003, 0.5 + 0.1003],
                        #            [0.0 + 0.1003, 0.5 + 0.1003, 0.5 + 0.1003]]),
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
                         ])]
        self.natom = [len(each) for each in self.position]
        self.ecut = 60 # in Ry
        self.k = [20, 20, 12, 20]
        self.pseudo = 'si.lda.lps'
        self.id = 'Si'
        self.zatom = 13
        self.N = 11
        self.create_files()
        self.cal_e_v()

    def create_INPUT(self, path):
        # create .inpt file for PROFESS
        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			blpstest\n\
calculation     scf\n\
esolver_type    ofdft\n\
ntype			1\n\
nbands			15\n\
symmetry		1\n\
pseudo_dir      /home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/real\n\
pseudo_rcut     16\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{0}\n\
scf_nmax			100\n\
#dft_functional  XC_GGA_X_PBE+XC_GGA_C_PBE\n\
\n\
\n\
#Parameters (3.Basis)\n\
basis_type		pw\n\
\n\
#OFDFT\n\
#of_kinetic     tf\n\
#of_kinetic     vw\n\
of_kinetic      wt\n\
of_method       tn\n\
of_conv         energy\n\
of_read_kernel  1\n\
of_kernel_file  ../../4PROFESS_KERNEL.dat\n\
of_tole         2e-6\n\
of_full_pw      1\n\
of_full_pw_dim  0\n".format(
                self.ecut))

    def create_STRU(self, path, cell, position, a):
        # create .ion file for PROFESS
        with open(path+"/STRU", 'w') as f:
            f.write("ATOMIC_SPECIES\n\
{0} {1} {2} blps\n\
\n\
LATTICE_CONSTANT\n\
{3}  // add lattice constant\n\
\n\
LATTICE_VECTORS\n".format(self.id, self.zatom, self.pseudo, a/self.bohr2a))
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('\nATOMIC_POSITIONS\nDirect\n\n{0}\n0.0\n{1}\n'.format(self.id, len(position)))
            for eachline in position:
                f.writelines('    ' + '%.12f' % each for each in eachline)
                f.write(' 1 1 1\n')

    def create_KPT(self, path, k):
        with open(path+"/KPT", 'w') as f:
            f.write("K_POINTS\n\
0\n\
Gamma\n\
1 1 1 0 0 0\n".format(k))
# {0} {0} {0} 0 0 0\n".format(k))

    def create_files(self):
        coef = np.linspace(0.9, 1.1, self.N)
        with open('jobs', 'w') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i].copy()
                    cell[:, -1] *= self.covera[i]
                    a = self.a[i] * coef[j]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_INPUT(path)
                    self.create_STRU(path, cell, self.position[i], a)
                    self.create_KPT(path, self.k[i])
                    job.write("cd $temp_path\{0} && ".format(path) + self.abacus + ' > log\n')



    def cal_e_v(self):
        # calculate e_v curve and return energy list
        e_list = np.zeros(self.N * len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}{1}/OUT.blpstest/running_scf.log'.format(self.structures[i], j)
                try:
                    get_total_e = subprocess.Popen(args="grep '!FINAL_ETOT_IS' %s|cut -d ' ' -f3" % outfile,
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
    abacus = Abacus()
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
