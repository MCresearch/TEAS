import numpy as np
import os
import subprocess
from scipy.optimize import curve_fit
import xlwt


class Profess:
    def __init__(self):
        self.pPROFESS = "pPROFESS"
        self.structures = ['bcc', 'fcc', 'sc', 'diamond']
        self.a = [3.3514983001409515, 4.216521201022214, 2.691873971428364, 5.604441917578286]
        # self.a = [3.431257319966733, 4.327318413514269, 2.737252293704602, 5.90960801642828]
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
        self.v = np.zeros(4)
        self.ecut = 900
        self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/reci/'
        self.pseudo = 'li.lda.recpot'
        # self.kernel_dir = "/home/dell/1_work/TKK_abacus/2_Li/1_kernels/"
        # self.kernelfiles = ['PROFESS_KERNEL.dat','1PROFESS_KERNEL.dat','3PROFESS_KERNEL.dat','4PROFESS_KERNEL.dat','5PROFESS_KERNEL.dat','8PROFESS_KERNEL.dat','11PROFESS_KERNEL.dat']
        # self.work_dir = "/home/dell/1_work/TKK_abacus/2_Li/2_bulk_properties/2_get_data"
        self.id = 'Li'
        self.N = 11
        self.workbook = xlwt.Workbook(encoding='utf-8')
        # for i in range(len(self.kernelfiles)):
            # self.kernelfile = self.kernelfiles[i]
            # self.a = self.constands[i]
        for j in range(4):
            self.v[j] = abs(np.dot(self.cell[j][0], np.cross(self.cell[j][1], self.cell[j][2])))
            self.v[j] *= self.covera[j] / self.natom[j]
        # os.chdir(self.work_dir)
        # os.mkdir(self.kernelfile[:-4])
        # os.chdir('./'+self.kernelfile[:-4])
        # os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
        os.system('cp {0} .'.format(self.pseudo_dir+self.pseudo))
        self.create_files()
        self.cal_e_v()
        # os.chdir(self.work_dir)
        self.workbook.save('data.xls')


    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WGC\nMAXI DEN 50\n'.format(
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
                '%BLOCK SPECIES_POT\n{0} {1}\n%END BLOCK SPECIES_POT\n'.format(self.id, self.pseudo))

    def create_files(self):
        coef = np.linspace(0.99, 1.01, self.N)
        with open('jobs', 'w') as job:
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i] * self.a[i] * coef[j]
                    cell[:, -1] *= self.covera[i]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_inpt(path+'/'+self.id+'.inpt')
                    self.create_ion(path+'/'+self.id+'.ion', cell, self.position[i])
                    job.write(self.pPROFESS + ' ' + path + '/' + self.id + '\n')

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
        coef = np.linspace(0.99, 1.01, self.N)
        v0 = np.zeros(len(self.structures))
        e0 = np.zeros_like(v0)
        B = np.zeros_like(v0)
        for i in range(len(self.structures)):
            for j in range(self.N):
                v_list[self.N * i + j] = self.v[i] * (coef[j] * self.a[i]) ** 3
            e0[i], v0[i], B[i] = self.fit(v_list[i*11:(i+1)*11], e_list[i*11:(i+1)*11], v_list[i*11], v_list[i*11+10])
        with open(self.id + '.curve', 'w') as f:
            f.write('Volume per atom(Ang^3)\tEnergy per atom(eV)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(self.N * len(self.structures)):
                f.write('%.12e\t%.12e\n' % (v_list[i], e_list[i]))
        worksheet = self.workbook.add_sheet("WT")
        worksheet.write(0, 0, 'structure')
        worksheet.write(0, 1, 'v0(Ang^3)')
        worksheet.write(0, 2, 'E0(eV)')
        worksheet.write(0, 3, 'B(GPa)')
        worksheet.write(0, 4, 'a0(Ang)')
        for i in range(len(self.structures)):
            worksheet.write(i+1, 0, self.structures[i])
            worksheet.write(i+1, 1, v0[i])
            worksheet.write(i+1, 2, e0[i])
            worksheet.write(i+1, 3, B[i])
            worksheet.write(i+1, 4, (v0[i]/self.v[i])**(1/3))

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
    profess = Profess()