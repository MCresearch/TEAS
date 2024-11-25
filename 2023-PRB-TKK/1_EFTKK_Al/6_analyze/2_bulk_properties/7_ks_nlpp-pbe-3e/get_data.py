import numpy as np
import os
import subprocess
from scipy.optimize import curve_fit


class Abacus:
    def __init__(self):
        self.abacus = "/home/dell/1_work/0_ABACUS_DEVELOP/build/abacus"
        self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.structures = ['fcc', 'bcc', 'hcp', 'sc']
        self.a = [4.06, 3.179285, 2.813975143, 2.65968504]
        self.covera = [1, 1, 1.6409235512, 1]
        self.cell = [np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]),
                     np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
                     np.array([[1, 0, 0], [-0.5, np.sqrt(3)/2, 0], [0, 0, 1]]),
                     np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])]
        self.position = [np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0]]),
                         np.array([[0, 0, 0], [2/3, 1/3, 0.5]]),
                         np.array([[0, 0, 0]])]
        self.natom = [len(each) for each in self.position]
        self.v = np.zeros(4)
        for i in range(4):
            self.v[i] = abs(np.dot(self.cell[i][0], np.cross(self.cell[i][1], self.cell[i][2])))
            self.v[i] *= self.covera[i] / self.natom[i]
        self.ecut = 60 # in Ry
        self.k = [20, 20, 12, 20]
        self.pseudo = 'Al.pbe-mt_fhi.UPF'
        self.id = 'Al'
        self.zatom = 13
        self.N = 11
        # self.create_files()
        # self.cal_e_v()

    def create_INPUT(self, path):
        # create .inpt file for PROFESS
        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			blpstest\n\
calculation     scf\n\
ntype			1\n\
nbands			15\n\
symmetry		1\n\
pseudo_dir      /home/dell/2_software/qe-6.6/pseudo/sg15_oncv_upf_2020-02-06\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{0}\n\
#dr2				1e-6\n\
scf_nmax		100\n\
\n\
\n\
#Parameters (3.Basis)\n\
basis_type		pw\n\
\n\
#Parameters (4.Smearing)\n\
smearing_method	gauss\n\
smearing_sigma	0.0074\n\
\n\
#Parameters (5.Mixing)\n\
mixing_type		pulay\n\
mixing_beta		0.7\n".format(
                self.ecut))

    def create_STRU(self, path, cell, position, a):
        # create .ion file for PROFESS
        with open(path+"/STRU", 'w') as f:
            f.write("ATOMIC_SPECIES\n\
{0} {1} {2}\n\
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
{0} {0} {0} 0 0 0\n".format(k))


    def create_files(self):
        coef = np.linspace(0.99, 1.01, self.N)
        with open('jobs', 'w') as job:
            job.write('temp_path=$PWD\n')
            for i in [0]:
            # for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i].copy()
                    cell[:, -1] *= self.covera[i]
                    a = self.a[i] * coef[j]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_INPUT(path)
                    self.create_STRU(path, cell, self.position[i], a)
                    self.create_KPT(path, self.k[i])
                    job.write("cd $temp_path\{0} && ".format(path) + self.abacus + '\n')

    def cal_e_v(self):
        # calculate e_v curve and return energy list
        e_list = np.zeros(self.N * len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in [0]:
        # for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}{1}/OUT.blpstest/running_scf.log'.format(self.structures[i], j)
                try:
                    get_total_e = subprocess.Popen(args="grep 'final etot' %s|cut -d ' ' -f5" % outfile,
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
        coef = np.linspace(0.99, 1.01, self.N)
        v0 = np.zeros(len(self.structures))
        e0 = np.zeros_like(v0)
        B = np.zeros_like(v0)
        for i in range(len(self.structures)):
            for j in range(self.N):
                v_list[self.N * i + j] = self.v[i] * (coef[j] * self.a[i]) ** 3
            e0[i], v0[i], B[i] = abacus.fit(v_list[i*11:(i+1)*11], e_list[i*11:(i+1)*11], v_list[i*11], v_list[i*11+10])
        with open(self.id + '.curve', 'w') as f:
            f.write('Volume per atom(Ang^3)\tEnergy per atom(eV)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(self.N * len(self.structures)):
                f.write('%.12e\t%.12e\n' % (v_list[i], e_list[i]))
        with open('data', 'w') as f:
            f.write('structure\tv0(Ang^3)\ta0(Ang)\tE0(eV)\tB(GPa)\n')
            for i in range(len(self.structures)):
                f.write('%s\t%.6f\t%.20f\t%.6f\t%.6f\n' % (self.structures[i], v0[i], (v0[i]/self.v[i])**(1/3), e0[i], B[i]))

    def dichotomy(self, f, a, b, eta=1e-7):
        # Solve the equation by dichotomy
        if f(a) * f(b) > 0:
            raise ValueError
        if abs(f(a)) <= eta:
            return a
        if abs(f(b)) <= eta:
            return b
        c = (a+b)/2
        while abs(f(c)) > eta:
            if f(c) * f(a) > 0:
                a = c
            else:
                b = c
            c = (a+b)/2
        return c

    def fit(self, volume, energy, v1, v2, v0=0):
        #volume in A^3, energy in eV
        def _energy(v, a, b, c, d, e, f):
            return a + b * v + c * v ** 2 + d * v ** 3 + e * v ** 4 + f * v ** 5

        def order1(v):
            return b + 2 * c * v + 3 * d * v ** 2 + 4 * e * v ** 3 + 5 * f * v ** 4
        try:
            a, b, c, d, e, f = curve_fit(_energy, volume, energy)[0]
        except:
            return 1e5, 1e5, 1e5
        try:
            v0 = self.dichotomy(order1, v1, v2)
        except:
            v0 = v0
        e0 = _energy(v0, a, b, c, d, e, f)
        B = (2 * c + 6 * d * v0 + 12 * e * v0 ** 2 + 20 * f * v0 ** 3) * \
            v0 * 160.2  # 160.2 for converting unit to GPa
        return e0, v0, B

if __name__ == '__main__':
    # input = '../get_ev/Al.curve'
    abacus = Abacus()
    # curve = np.loadtxt(input, skiprows=2).T
    # for i in range(4):
    #     e0, v0, B = abacus.fit(curve[0][i*11:(i+1)*11], curve[1][i*11:(i+1)*11], curve[0][i*11], curve[0][i*11+10])
    #     abacus.a[i] = (v0/abacus.v[i]) ** (1/3)
    abacus.create_files()
    abacus.cal_e_v()
