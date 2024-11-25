import numpy as np
import os
import subprocess

#from scipy.optimize import curve_fit


class PROFESS:
    def __init__(self, minNatom, maxNatom, ecut, stru_dir, kinetic, use_kernel, kernel_dir, kernel_prefix, xc):

        self.minNatom = minNatom
        self.maxNatom = maxNatom
        self.ecut = ecut
        if self.ecut == 0: self.ecut = 800
        self.kinetic = kinetic
        self.use_kernel = use_kernel
        self.kernel_dir = kernel_dir
        self.sub_dir = stru_dir.split('/')[-1]

        self.work_dir = os.getcwd()
        self.pPROFESS = "pPROFESS"
        # self.bohr2a = 0.529177249
        # self.ry2ev = "13.6056923"
        self.prefix = kernel_prefix
        self.kernelfile = '{0}PROFESS_KERNEL.dat'.format(kernel_prefix)

        self.xc = xc

        self.structures = []
        self.id = []
        self.natom = []
        self.a = []
        self.cell = []
        self.position = []
        self.find_stru(stru_dir)
        self.Nstru = len(self.a)
        print("Total {0} structures found:\n".format(self.Nstru))
        print(self.structures)
        self.covera = [1.] *  self.Nstru

        self.elements = ["Al", "Mg", "Li"]
        self.pseudo_dir = ""
        self.pseudo = dict()
        if self.xc == 0:
            self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/reci/'
            # self.pseudo_dir = '/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/PROFESS3.0/BLPSLibrary-master/LDA/real'
            self.pseudo = {"Al":'al.lda.recpot', "Mg":'mg.lda.recpot', "Li":'li.lda.recpot'}
        elif self.xc == 1: 
            self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/GGA/reci/'
            # self.pseudo_dir = '/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/PROFESS3.0/BLPSLibrary-master/GGA/real'
            self.pseudo = {"Al":'al.gga.recpot', "Mg":'mg.gga.recpot', "Li":'li.gga.recpot'}
        if use_kernel:
            os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
        for pp in self.pseudo.values():
            os.system('cp {0} .'.format(self.pseudo_dir+pp))
        self.create_files()
        self.calculate()

    def find_stru(self, stru_dir):
        all_stru = os.listdir(stru_dir)
        for each in all_stru:
            file = '{0}/{1}'.format(stru_dir, each)
            if self.read_stru(file):
                self.structures.append(each.split('.')[0])
    
    def read_stru(self, file):
        with open(file, 'r') as stru_file:
            title = stru_file.readline()
            title = title.replace('Li', '')
            title = title.replace('Mg', '')
            title = title.replace('Al', '')
            totnatom = sum([int(_) for _ in title.split()])
            if totnatom > self.maxNatom or totnatom <= self.minNatom:
                return False
            else:
                self.a.append(float(stru_file.readline()))
                self.cell.append(
                    np.array([[float(_) for _ in stru_file.readline().split()] for j in range(3)])
                )
                self.id.append(stru_file.readline().split())
                self.natom.append(np.array(stru_file.readline().split(), dtype=np.int64))
                stru_file.readline()
                self.position.append(
                    np.array([[float(_) for _ in stru_file.readline().split()[:3]] for j in range(totnatom)])
                )
                return True

    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('MINI DEN\nECUT	{0}\nKINE	{1}\nMAXI DEN 50\nPRIN DEN\n'.format(
            # f.write('MINI CEL\nMETH ION NON\nECUT	{0}\nKINE	{1}\nMAXI DEN 50\n'.format(
                self.ecut, self.kinetic))
            if self.use_kernel:
                f.write('KERN {0}\n'.format(self.kernelfile))
            if self.xc == 1:
                f.write('EXCH   PBE\n')

    def create_ion(self, ionfile, cell, position, id, natom):
        # create .ion file for PROFESS
        with open(ionfile, 'w') as f:
            f.write('%BLOCK LATTICE_CART\n')
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('%END BLOCK LATTICE_CART\n')
            f.write('%BLOCK POSITIONS_FRAC\n')

            elements = []
            for j in range(len(id)):
                elements += [id[j]] * natom[j]

            for i, eachline in enumerate(position):
                f.write(elements[i])
                f.writelines('    ' + '%.12f' % each for each in eachline)
                f.write('\n')
            f.write('%END BLOCK POSITIONS_FRAC\n')
            f.write('%BLOCK SPECIES_POT\n')
            for ion in id:
                f.write('{0} {1}\n'.format(ion, self.pseudo[ion]))
            f.write('%END BLOCK SPECIES_POT\n')

    def create_files(self):
        if not os.path.exists('log'): os.mkdir('log')
        with open('jobs', 'w') as job:
            for i in range(self.Nstru):
                cell = self.cell[i] * self.a[i]
                cell[:, -1] *= self.covera[i]
                path = self.structures[i]
                os.mkdir(path)
                self.create_inpt(path+'/'+path+'.inpt')
                self.create_ion(path+'/'+path+'.ion', cell, self.position[i], self.id[i], self.natom[i])
                job.write('{0} {2}/{3} > log/{4} && echo \"{4} done\"\n'.format(self.pPROFESS, self.sub_dir, path, path, path))

    def calculate(self):
        # perform structure optimization
        # e_list = np.zeros((2,self.N * self.Nstru))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        # for i in range(self.Nstru):
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
#         # for i in range(self.Nstru):
#         #     v = abs(np.dot(self.cell[i][0], np.cross(self.cell[i][1], self.cell[i][2])))
#         #     v *= self.a[i] ** 3 * self.covera[i] / self.natom[i]
#         #     for j in range(self.N):
#         #         v_list[self.N * i + j] = v * coef[j] ** 3
#         with open(self.id + '_target', 'w') as f:
#             f.write('1-5:fcc al 6-10:hcp 11-15:bcc 15-20:sc\n\
# Energy(eV)                   Hartree             natoms\n')
#             for i in range(self.N * self.Nstru):
#                 f.write('%.12e\t%.12e\t%d\n' % (e_list[0,i], e_list[1,i], self.natom[i//5]))

if __name__ == '__main__':
    stru_dir = '/home/dell/1_work/7_ABACUS_ML_OF/1_test/0_generate_data/2_ks-pbe/0_data_set/a_alloys/0_structures/'
    # stru_dir = '/home/dell/1_work/TKK_abacus/5_MgAl_alloy/4_28config/0_structures/Mg-Al'
    kernel_dir = "/home/dell/1_work/TKK_abacus/1_EFTKK_Al/5_results/"
    minNatom = int(input("Please enter the min number of atoms:"))
    maxNatom = int(input("Please enter the max number of atoms:"))
    ecut = float(input("Please enter the ECUT(in eV):"))
    kinetic = input("Please enter the name of KEDF:")
    use_kernel = int(input("Do you want to use kernel?"))
    xc = int(input("Please enter the type of XC (0:LDA, 1:PBE):"))
    kernel_prefix = 0
    work_dir = os.getcwd()
    if use_kernel: kernel_prefix = int(input("Please enter the prefix of kernel file (one of 1, 3, 4, 5, 8, 11):"))
    
    for sub_dir in ["0_LiMg", "1_MgAl", "2_LiAl", "3_LiMgAl", "4_single"]:
        total_dir = stru_dir + sub_dir
        if not os.path.exists(sub_dir): os.mkdir(sub_dir)
        os.chdir(sub_dir)
        profess = PROFESS(minNatom, maxNatom, ecut, total_dir, kinetic, use_kernel, kernel_dir, kernel_prefix, xc)
        os.chdir(work_dir)
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
