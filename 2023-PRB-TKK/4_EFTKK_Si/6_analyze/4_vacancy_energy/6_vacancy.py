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
            self.constants[i] = worksheet.cell_value(1,4)
            self.ground_energies[i] = worksheet.cell_value(1,2)

        # self.structures = ['111', '211', '222', '333', '444', '555', '666', '777']
        self.structures = ['333', '444', '555', '666', '777']
        # self.structures = ['111', '211', '222']
        self.unit_cell = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.unit_position = np.array([[0, 0, 0], [0.5, 0.5 ,0], [0.5, 0, 0.5], [0, 0.5, 0.5], [0.25, 0.25, 0.25], [0.75, 0.75, 0.25], [0.75, 0.25, 0.75], [0.25, 0.75, 0.75]])
        self.cell = []
        self.vacancy_position = []
        self.CD_position = []
        self.get_cells()
        self.natom = np.array([len(each) for each in self.vacancy_position])

        self.ecut = 1000
        self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/reci/'
        self.pseudo = 'si.lda.recpot'
        self.kernel_dir = "/home/dell/1_work/TKK_abacus/4_EFTKK_Si/5_results/"
        self.prefix = [1, 4, 7, 10, 13]
        self.kernelfiles = [str(_) + 'PROFESS_KERNEL.dat' for _ in self.prefix]
        self.work_dir = "/home/dell/1_work/TKK_abacus/4_EFTKK_Si/6_analyze/4_vacancy_energy"
        self.id = 'Si'

        self.workbook = xlwt.Workbook(encoding='upf-8')
        self.evacancy_sheet = self.workbook.add_sheet('vacancy energy')
        self.etot_sheet = self.workbook.add_sheet('total energy')
        self.evacancy_sheet.write(0, 0, 'kernel')
        self.etot_sheet.write(0, 0, 'kernel')
        for i in range(len(self.structures)):
            self.evacancy_sheet.write(0, i+1, self.structures[i])
            self.etot_sheet.write(0, i+1, self.structures[i])

        for i in range(len(self.kernelfiles)):
            self.kernelfile = self.kernelfiles[i]
            self.ground_e = self.ground_energies[i]
            self.a = self.constants[i]
            os.chdir(self.work_dir)
            # os.mkdir(self.kernelfile[:-4])
            os.chdir('./'+self.kernelfile[:-4])
            os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
            os.system('cp {0} .'.format(self.pseudo_dir+self.pseudo))
            self.create_files()
            etot_list, evacancy_list = self.calculate()
            self.evacancy_sheet.write(i+1, 0, self.kernelfile[:-4])
            self.etot_sheet.write(i+1, 0, self.kernelfile[:-4])
            for j in range(len(self.structures)):
                self.evacancy_sheet.write(i+1, j+1, evacancy_list[j])
                self.etot_sheet.write(i+1, j+1, etot_list[j])
        os.chdir(self.work_dir)
        self.workbook.save('vacancy_energy_CD_new.xls')

    def get_cells(self):
        for each in self.structures:
            size = np.array([int(i) for i in each])
            cell, position = self.extend(size, self.unit_cell, self.unit_position)
            self.cell.append(cell)
            self.vacancy_position.append(position[1:])
            self.CD_position.append(position[:])

    def extend(self, size, unit_cell, unit_position):
        cell = unit_cell * size.reshape(3,1)
        temp_position = unit_position / size
        position = temp_position.copy()
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    if [i,j,k] != [0,0,0]:
                        translation = np.array([i/size[0], j/size[1], k/size[2]])
                        position = np.vstack([position, temp_position+translation])
        return cell, position

    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WT\nKERN {1}\nMAXI DEN 50\n'.format(
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
                self.create_ion(path+'/'+self.id+'.ion', cell, self.vacancy_position[i])
                job.write(self.pPROFESS + ' ' + path + '/' + self.id + '\n')
                self.create_inpt(path+'/'+self.id+'_CD.inpt')
                self.create_ion(path+'/'+self.id+'_CD.ion', cell, self.CD_position[i])
                job.write(self.pPROFESS + ' ' + path + '/' + self.id + '_CD\n')

    def calculate(self):
        # perform structure optimization
        # e_list = np.zeros((2,self.N * len(self.structures)))
        etot_list = np.zeros(len(self.structures))
        eCD_list = np.zeros(len(self.structures))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for i in range(len(self.structures)):
            outfile = './{0}/{1}.out'.format(self.structures[i], self.id)
            CDfile = './{0}/{1}_CD.out'.format(self.structures[i], self.id)
            try:
                get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read())  # maybe wrong
                get_CD_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % CDfile,
                                            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                CD_e = float(get_CD_e.stdout.read())  # maybe wrong
            except:
                total_e = 1e5
                CD_e = 1e5
            etot_list[i] = total_e
            eCD_list[i] = CD_e
        evacancy_list = etot_list - self.natom/(self.natom+1) * eCD_list
        return etot_list, evacancy_list

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
