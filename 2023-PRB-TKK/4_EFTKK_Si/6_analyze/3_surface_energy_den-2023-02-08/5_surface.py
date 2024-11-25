import numpy as np
import os
import subprocess
import xlrd
import xlwt
#from scipy.optimize import curve_fit


class PROFESS:
    def __init__(self):
        self.pPROFESS = "pPROFESS"
        # self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.data_file = '/home/dell/1_work/TKK_abacus/4_EFTKK_Si/6_analyze/2_bulk_properties/2_get_data/data.xls'
        data = xlrd.open_workbook(self.data_file)
        sheets = data.sheet_names()
        self.constants = [0] * len(sheets)
        self.ground_energies = [0] * len(sheets)
        for i in range(len(sheets)):
            worksheet = data.sheet_by_name(sheets[i])
            # self.constants[i] = worksheet.cell_value(1,4)
            self.constants[i] = 5.408352067736602
            self.ground_energies[i] = worksheet.cell_value(1,2)

        self.structures = ['cd100']
        self.cell = [np.array([[1/2**0.5, 0, 0], [0, 1/2**0.5, 0], [0, 0, 4]])]
        self.position = [np.array([[0, 0, 0], [0.5, 0 ,1/16], [0.5, 0.5, 2/16], [0, 0.5, 3/16], [0, 0, 4/16], [0.5, 0, 5/16], [0.5, 0.5, 6/16], [0, 0.5, 7/16], [0, 0, 8/16]])]
        self.natom = np.array([len(each) for each in self.position])
        self.area = np.array([np.linalg.norm(np.cross(cell[0], cell[1])) for cell in self.cell])

        self.ecut = 3600
        self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/reci/'
        self.pseudo = 'si.lda.recpot'
        self.kernel_dir = "/home/dell/1_work/TKK_abacus/4_EFTKK_Si/5_results/"
        self.prefix = [1, 4, 7, 10, 13]
        self.kernelfiles = [str(_) + 'PROFESS_KERNEL.dat' for _ in self.prefix]
        self.work_dir = "/home/dell/1_work/TKK_abacus/4_EFTKK_Si/6_analyze/3_surface_energy_den-2023-02-08/"
        self.id = 'Si'

        self.workbook = xlwt.Workbook(encoding='upf-8')
        # self.eslab_sheet = self.workbook.add_sheet("e_slab(mJ/m^2)")
        self.eslab_sheet = self.workbook.add_sheet("e_slab")
        self.etot_sheet = self.workbook.add_sheet('total energy')
        self.eslab_sheet.write(0, 0, 'kernel')
        self.etot_sheet.write(0, 0, 'kernel')
        for i in range(len(self.structures)):
            self.eslab_sheet.write(0, i+1, self.structures[i])
            self.etot_sheet.write(0, i+1, self.structures[i])

        for i in range(len(self.kernelfiles)):
            self.kernelfile = self.kernelfiles[i]
            self.ground_e = self.ground_energies[i]
            self.a = self.constants[i]
            os.chdir(self.work_dir)
            os.mkdir(self.kernelfile[:-4])
            os.chdir('./'+self.kernelfile[:-4])
            os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
            os.system('cp {0} .'.format(self.pseudo_dir+self.pseudo))
            self.create_files()
            etot_list, eslab_list = self.calculate()
            self.eslab_sheet.write(i+1, 0, self.kernelfile[:-4])
            self.etot_sheet.write(i+1, 0, self.kernelfile[:-4])
            for j in range(len(self.structures)):
                self.eslab_sheet.write(i+1, j+1, eslab_list[j])
                self.etot_sheet.write(i+1, j+1, etot_list[j])
        os.chdir(self.work_dir)
        self.workbook.save('surface_energy-2023-02-08.xls')


    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WT\nKERN {1}\nMAXI DEN 50\nPRIN DEN\n'.format(
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
                cell = self.cell[i] * self.a
                path = self.structures[i]
                os.mkdir(path)
                self.create_inpt(path+'/'+self.id+'.inpt')
                self.create_ion(path+'/'+self.id+'.ion', cell, self.position[i])
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
        eslab_list = (etot_list - self.natom * self.ground_e) * 1.60217663410e-16 /(2 * self.area * self.a ** 2 * 1e-20)
        return etot_list, eslab_list

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
