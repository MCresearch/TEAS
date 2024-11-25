import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit


class Abacus:
    def __init__(self):
        self.abacus = "/home/dell/1_work/5_ABACUS_DEBUG/0_abacus_develop/abacus-develop/build/abacus"
        self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.a = 3.9850582555647662
        self.ground = -57.934369

        self.structures = ['fcc100', 'fcc110', 'fcc111']
        self.cell = [np.array([[1/2**0.5, 0, 0], [0, 1/2**0.5, 0], [0, 0, 6]]),
                     np.array([[1/2**0.5, 0, 0], [0, 1, 0], [0, 0, 2**0.5/4*16]]),
                     np.array([[1/2**0.5, 0, 0], [1/2**0.5/2, 1/2**0.5/2*3**0.5, 0], [0, 0, 3**0.5/3*12]])]
        self.position = [np.array([[0, 0, 0], [0.5, 0.5 ,1/12], [0, 0, 2/12], [0.5, 0.5, 3/12], [0, 0, 4/12], [0.5, 0.5, 5/12], [0, 0, 6/12]]),
                         np.array([[0, 0, 0], [0.5, 0.5 ,1/16], [0, 0, 2/16], [0.5, 0.5, 3/16], [0, 0, 4/16], [0.5, 0.5 ,5/16], [0, 0, 6/16], [0.5, 0.5 ,7/16], [0, 0, 8/16]]),
                         np.array([[0, 0, 0], [1/3, 1/3, 1/12], [2/3, 2/3, 2/12], [0, 0, 3/12], [1/3, 1/3, 4/12], [0, 0, 5/12], [1/3, 1/3, 6/12]])]
        self.natom = np.array([len(each) for each in self.position])
        self.area = np.array([np.linalg.norm(np.cross(cell[0], cell[1])) for cell in self.cell])

        self.ecut = 60 # in Ry
        self.pseudo = 'al.lda.lps'
        self.id = 'Al'
        self.zatom = 13
        self.full_pw = 1
        self.full_pw_dim = 0

        self.create_files()
        self.cal()

    def create_INPUT(self, path):
        # create .inpt file for PROFESS
        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			oftest\n\
calculation     scf\n\
esolver_type    ofdft\n\
symmetry		1\n\
pseudo_dir      /home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/real\n\
pseudo_rcut     16\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{0}\n\
scf_nmax        200\n\
# dft_functional  XC_GGA_X_PBE+XC_GGA_C_PBE\n\
\n\
\n\
#Parameters (3.Basis)\n\
basis_type		pw\n\
\n\
# #Parameters (4.Smearing)\n\
# smearing_method		gauss\n\
# smearing_sigma			0.0074\n\
# \n\
# #Parameters (5.Mixing)\n\
# mixing_type		pulay\n\
# mixing_beta		0.7\n\
".format(
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

    def create_KPT(self, path):
        with open(path+"/KPT", 'w') as f:
            f.write("K_POINTS\n\
0\n\
Gamma\n\
1 1 1 0 0 0\n")

    def create_files(self):
        with open('jobs', 'w') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                path = self.structures[i]
                os.mkdir(path)
                self.create_INPUT(path)
                self.create_STRU(path, self.cell[i], self.position[i], self.a)
                self.create_KPT(path)
                job.write("cd $temp_path\{0} && ".format(path) + self.abacus + ' > log\n')


    def cal(self):
        e_list = np.zeros(len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in range(len(self.structures)):
            outfile = './{0}/OUT.oftest/running_scf.log'.format(self.structures[i])
            try:
                get_total_e = subprocess.Popen(args="grep '!FINAL_ETOT_IS' %s|cut -d ' ' -f3" % outfile,
                                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read())  # maybe wrong
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
            e_list[i] = total_e
            # e_list[1][i] = kinetic_e
            # e_list[2][i] = hartree_e
        eslab_list = (e_list - self.natom * self.ground) * 1.60217663410e-16 /(2 * self.area * self.a ** 2 * 1e-20)
        with open('data', 'w') as f:
            f.write('Total energy(eV)\tsurface energy(mJ/m^2)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(len(self.structures)):
                f.write('%.12e\t%.12e\n' % (e_list[i], eslab_list[i]))



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
