import numpy as np
import os
import subprocess
import xlrd
import matplotlib.pyplot as plt
# import xlwt
#from scipy.optimize import curve_fit


class PROFESS:
    def __init__(self):
        self.pPROFESS = "pPROFESS"
        # self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.data_file = '/home/dell/1_work/TKK_abacus/1_EFTKK_Al/6_analyze/2_bulk_properties/2_get_data/data.xls'
        data = xlrd.open_workbook(self.data_file)
        sheets = data.sheet_names()
        self.constants = [0] * len(sheets)
        for i in range(len(sheets)):
            worksheet = data.sheet_by_name(sheets[i])
            self.constants[i] = worksheet.cell_value(1,4)

        self.id = 'Al'
        self.structures = ['{0}{1}'.format(self.id, _) for _ in range(23)]
        self.cell1 = np.array([1/2**0.5, 1/2**0.5, 3**0.5/3*20, 90, 90, 60])
        self.cell2 = np.array([1/2**0.5, 1/2**0.5, 3**0.5/3*22, 90, 90, 60])
        self.cell3 = np.array([1/2**0.5, 1/2**0.5, 3**0.5/3*21, 90, 90, 60])
        self.set_position()
        self.natom = np.array([len(each) for each in self.position])
        self.area = 3**0.5 / 2

        self.ecut = 800
        self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/reci/'
        self.pseudo = 'al.lda.recpot'
        self.kernel_dir = "/home/dell/1_work/TKK_abacus/1_EFTKK_Al/5_results/"
        self.kernelfiles = ['1PROFESS_KERNEL.dat','3PROFESS_KERNEL.dat','4PROFESS_KERNEL.dat','5PROFESS_KERNEL.dat','8PROFESS_KERNEL.dat','11PROFESS_KERNEL.dat']

        self.rcut = [8,12, 16, 20,24, 32]
        self.work_dir = "/home/dell/1_work/TKK_abacus/1_EFTKK_Al/6_analyze/5_stacking_energy"

        self.stack_file = open('stacking_energy.txt', 'w')
        self.etot_file = open('total_energy.txt', 'w')
        self.stack_file.writelines('%s\t' % each[:-4] for each in self.kernelfiles)
        self.stack_file.write('\n')
        self.etot_file.writelines('%s\t' % each[:-4] for each in self.kernelfiles)
        self.etot_file.write('\n')
        plt.figure()
        ax = plt.gca() 
        ax.tick_params(top=True,right = True)
        self.colors = ['dodgerblue', 'k', 'r', 'orange','b','g']
        
        for i in range(len(self.kernelfiles)):
            self.kernelfile = self.kernelfiles[i]
            self.a = self.constants[i]
            self.s = (self.a / 2 **0.5) ** 2 * self.area
            os.chdir(self.work_dir)
            os.mkdir(self.kernelfile[:-4])
            os.chdir('./'+self.kernelfile[:-4])
            os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
            os.system('cp {0} .'.format(self.pseudo_dir+self.pseudo))
            self.create_files()
            etot_list, stack_list = self.calculate()
            self.etot_file.writelines('%.12e\t' % _ for _ in etot_list)
            self.etot_file.write('\n')
            self.stack_file.writelines('%.6e\t' % _ for _ in stack_list)
            self.stack_file.write('\n')
            plt.plot(np.linspace(0,1,11), stack_list[1:12], self.colors[i], linestyle='--',marker='^', label='rcut=%.1f' % self.rcut[i],markerfacecolor='w')
            plt.plot(np.linspace(1,2,11), stack_list[12:], self.colors[i], linestyle='--',marker='^', markerfacecolor='w')
        os.chdir(self.work_dir)
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

    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WT\nKERN {1}\nMINI CEL\nMETH CEL 3\nMETH ION CON\n'.format(
                self.ecut, self.kernelfile))

    def create_ion(self, ionfile, cell, position):
        # create .ion file for PROFESS
        with open(ionfile, 'w') as f:
            f.write('%BLOCK LATTICE_ABC\n')
            f.writelines('%.10f' % each + '    ' for each in cell)
            f.write('\n')
            f.write('%END BLOCK LATTICE_ABC\n')
            f.write('%BLOCK POSITIONS_FRAC\n')
            for eachline in position:
                f.write(self.id)
                f.writelines('    ' + '%.10f' % each for each in eachline)
                f.write('\n')
            f.write('%END BLOCK POSITIONS_FRAC\n')
            f.write('%BLOCK ION_OPTIMIZATION\n')
            for i in range(len(position)):
                f.write('0 0 1\n')
            f.write('%END BLOCK ION_OPTIMIZATION\n')
            f.write(
                '%BLOCK SPECIES_POT\n{0} {1}\n%END BLOCK SPECIES_POT\n'.format(self.id, self.pseudo))

    def create_files(self):
        with open('jobs', 'w') as job:
            for i in range(len(self.structures)):
                cell1 = self.cell1.copy()
                cell1[:3] *= self.a
                cell2 = self.cell2.copy()
                cell2[:3] *= self.a
                cell3 = self.cell3.copy()
                cell3[:3] *= self.a
                path = self.structures[i]
                os.mkdir(path)
                self.create_inpt(path+'/'+self.id+'.inpt')
                if self.natom[i] == 20:
                    self.create_ion(path+'/'+self.id+'.ion', cell1, self.position[i])
                elif self.natom[i] == 22:
                    self.create_ion(path+'/'+self.id+'.ion', cell2, self.position[i])
                elif self.natom[i] == 21:
                    self.create_ion(path+'/'+self.id+'.ion', cell3, self.position[i])
                job.write(self.pPROFESS + ' ' + path + '/' + self.id + '\n')

    def calculate(self):
        # perform structure optimization
        # e_list = np.zeros((2,self.N * len(self.structures)))
        etot_list = np.zeros(len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in range(len(self.structures)):
            outfile = './{0}/{1}.out'.format(self.structures[i], self.id)
            try:
                get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
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
