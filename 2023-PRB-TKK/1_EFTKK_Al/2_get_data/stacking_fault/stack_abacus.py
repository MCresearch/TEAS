import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit


class Abacus:
    def __init__(self):
        self.abacus = "abacus"
        self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"

        self.a = 3.97010529048443894240
        self.id = 'Al'
        self.structures = ['{0}{1}'.format(self.id, _) for _ in range(23)]
        self.cell1 = np.array([[1/2**0.5,0,0], [1/2**0.5*1/2,1/2**0.5*3**0.5/2,0], [0,0,3**0.5/3*20]])
        self.cell2 = np.array([[1/2**0.5,0,0], [1/2**0.5*1/2,1/2**0.5*3**0.5/2,0], [0,0,3**0.5/3*22]])
        self.cell3 = np.array([[1/2**0.5,0,0], [1/2**0.5*1/2,1/2**0.5*3**0.5/2,0], [0,0,3**0.5/3*21]])
        self.set_position()
        self.natom = np.array([len(each) for each in self.position])
        self.s = 3**0.5 / 2 * self.a ** 2

        self.ecut = 40 # in Ry
        self.pseudo = 'al.lda.lps'
        self.id = 'Al'
        self.zatom = 27
        self.N = 5
        # self.create_files()
        # self.cal_e_v()

        self.stack_file = open('stacking_energy.txt', 'w')
        self.etot_file = open('total_energy.txt', 'w')
        
        plt.figure()
        ax = plt.gca() 
        ax.tick_params(top=True,right = True)    

        self.create_files()
        etot_list, stack_list = self.cal_e_v()
        self.etot_file.writelines('%.12e\t' % _ for _ in etot_list)
        self.etot_file.write('\n')
        self.stack_file.writelines('%.6e\t' % _ for _ in stack_list)
        self.stack_file.write('\n')
        plt.plot(np.linspace(0,1,11), stack_list[1:12], linestyle='--',marker='^', label='KSDFT',markerfacecolor='w')
        plt.plot(np.linspace(1,2,11), stack_list[12:], linestyle='--',marker='^', markerfacecolor='w')

        self.etot_file.close()
        self.stack_file.close()
        plt.xlabel('d',fontsize=14)
        plt.ylabel('E$(mJ/m^2)$',fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(frameon=False,ncol=2,fontsize=12)
        plt.xlim(0,2)
        plt.ylim(0,)
        plt.tight_layout()
        # plt.show()
        plt.savefig('stacking.png', dpi=350)

    def set_position(self):
        A = np.array([0, 0, 0])
        B = np.array([1/3, 1/3, 0])
        C = np.array([2/3, 2/3, 0])
        self.position = [np.vstack([A,B,C,A,B,C,A,B,C,A,B,C,A,B,C,A,B,C,A,B,C])]
        self.position[0][:,-1] += np.linspace(0,20/21,21)
        position = np.vstack([A,B,C,A,B,C,A,B,C,A,C,B,A,C,B,A,C,B,A,C])
        position[:,-1] += np.linspace(0, 19/20, 20)
        for i in range(11):
            self.position.append(position.copy())
            position[5:14,0:2] += 1/30
        position = np.vstack([A,B,C,A,B,A,B,C,A,B,C,B,A,C,B,A,B,A,C,B,A,C])
        position[:,-1] += np.linspace(0, 21/22, 22)
        for i in range(11):
            self.position.append(position.copy())
            position[6:15,0:2] += 1/30

    def create_INPUT(self, path):
        # create .inpt file for PROFESS
        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			stack\n\
calculation     cell-relax\n\
ntype			1\n\
nbands			45\n\
symmetry		1\n\
pseudo_dir      /home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/real\n\
pseudo_type     blps\n\
pseudo_rcut     16\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{0}\n\
dr2				1e-6\n\
niter			100\n\
nstep           50\n\
force_thr_ev    0.000001\n\
fixed_axes      ab\n\
\n\
\n\
#Parameters (3.Basis)\n\
basis_type		pw\n\
\n\
#Parameters (4.Smearing)\n\
smearing		gauss\n\
sigma			0.0074\n\
\n\
#Parameters (5.Mixing)\n\
mixing_type		pulay\n\
mixing_beta		0.1\n".format(
                self.ecut))

    def create_STRU(self, path, cell, position):
        # create .ion file for PROFESS
        with open(path+"/STRU", 'w') as f:
            f.write("ATOMIC_SPECIES\n\
{0} {1} {2}\n\
\n\
LATTICE_CONSTANT\n\
{3}  // add lattice constant\n\
\n\
LATTICE_VECTORS\n".format(self.id, self.zatom, self.pseudo, self.a/self.bohr2a))
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('\nATOMIC_POSITIONS\nDirect\n\n{0}\n0.0\n{1}\n'.format(self.id, len(position)))
            for eachline in position:
                f.writelines('    ' + '%.12f' % each for each in eachline)
                f.write(' 0 0 1\n')

    def create_KPT(self, path):
        with open(path+"/KPT", 'w') as f:
            f.write("K_POINTS\n\
0\n\
Gamma\n\
10 10 1 0 0 0.5\n")

    def create_files(self):
        with open('jobs', 'w') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                path = self.structures[i]
                os.mkdir(path)
                self.create_INPUT(path)
                self.create_KPT(path)
                if self.natom[i] == 20:
                    self.create_STRU(path, self.cell1, self.position[i])
                elif self.natom[i] == 22:
                    self.create_STRU(path, self.cell2, self.position[i])
                elif self.natom[i] == 21:
                    self.create_STRU(path, self.cell3, self.position[i])
                job.write("cd $temp_path\{0} && ".format(path) + self.abacus + '\n')

    def cal_e_v(self):
        # calculate e_v curve and return energy list
        etot_list = np.zeros(len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in range(len(self.structures)):
            outfile = './{0}/OUT.stack/running_scf.log'.format(self.structures[i])
            try:
                get_total_e = subprocess.Popen(args="grep 'final etot' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read())  # maybe wrong
            except:
                total_e = 1e5
            etot_list[i] = total_e
        stack_list = etot_list.copy()
        stack_list[1:12] = (etot_list[1:12] - etot_list[1]) * 1.60217663410e-16 / self.s / 2 / 1e-20
        stack_list[12:] = (etot_list[12:] - etot_list[1] - 2*etot_list[0]/21) * 1.60217663410e-16 / self.s / 2 / 1e-20
        return etot_list, stack_list

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
