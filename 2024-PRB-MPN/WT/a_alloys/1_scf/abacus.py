import numpy as np
import os
import subprocess

#from scipy.optimize import curve_fit


class ABACUS:
    def __init__(self, minNatom, maxNatom, ecut, stru_dir, calculation, esolver_type, kinetic, use_kernel, kernel_dir, kernel_prefix, xc):

        self.ry2ev = 13.6056923
        self.bohr2a = 0.529177249
        self.minNatom = minNatom
        self.maxNatom = maxNatom
        self.ecut = ecut / self.ry2ev
        if self.ecut == 0: self.ecut = 800 / self.ry2ev
        self.calculation = calculation
        self.esolver_type = esolver_type
        self.kinetic = kinetic
        self.use_kernel = use_kernel
        self.kernel_dir = kernel_dir
        self.sub_dir = stru_dir.split('/')[-1]

        self.work_dir = os.getcwd()
        self.abacus = "/home/dell/1_work/7_ABACUS_ML_OF/0_abacus-develop/build_gnu_openmp_ttt_feg3_double/abacus"
        # self.abacus = "/home/dell/1_work/1_ABACUS_OFDFT/0_abacus-develop/build/abacus"
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

        # self.kspacing = 0.05
        self.symmetry_prec = 2e-4
        self.k = [1] * self.Nstru
        for i in range(self.Nstru):
            n = sum(self.natom[i])
            if n <= 4:
                self.k[i] = 20
            else:
                self.k[i] = int(20 / (n/4)**(1/3))
            if self.k[i] == 0: self.k[i] = 1
            if self.esolver_type == "ofdft" : self.k[i] = 1

        # self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/real/'
        self.elements = ["Al", "Mg", "Li"]
        self.mass = {"Al":26.981, "Mg":24.305, "Li":6.941}
        self.zvan = {"Al":3, "Mg":2, "Li":1}
        self.pseudo_dir = ""
        self.pseudo = dict()
        if self.xc == 0:
            self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/LDA/real/'
            self.pseudo = {"Al":'al.lda.lps', "Mg":'mg.lda.lps', "Li":'li.lda.lps.1'}
        elif self.xc == 1: 
            self.pseudo = {"Al":'al.gga.psp', "Mg":'mg.gga.psp', "Li":'li.gga.psp.1'}
            self.pseudo_dir = '/home/dell/2_software/g_pPROFESS/BLPSLibrary-master/GGA/real/'

        if use_kernel:
            os.system('cp {0} .'.format(self.kernel_dir+self.kernelfile))
        # for pp in self.pseudo.values():
        #     os.system('cp {0} .'.format(self.pseudo_dir+pp))
        self.create_files()
        # self.calculate()

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

    def create_INPUT(self, path, ntype, nbands):
        # create INPUT file for ABACUS
        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			oftest\n\
calculation     {0}\n\
esolver_type   {1}\n\
symmetry		1\n\
#symmetry_prec   {2}\n\
pseudo_dir      {3}\n\
pseudo_rcut     16\n\
nbands          {4}\n\
out_stru        1\n\
cal_force       1\n\
cal_stress      1\n\
out_chg         1\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{5}\n\
scf_nmax        100\n\
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
                self.calculation, self.esolver_type, self.symmetry_prec, self.pseudo_dir, nbands, self.ecut))

            if self.xc == 1:
                f.write("\n\
dft_functional  XC_GGA_X_PBE+XC_GGA_C_PBE\n")

            if self.calculation == "cell-relax":
                f.write("\n\
relax_nmax  100\n")

            if self.esolver_type == "ofdft":
                f.write("\n#OFDFT\n\
ntype           {0}\n\
of_kinetic      {1}\n\
of_method       tn\n\
of_conv         energy\n\
of_tole         2e-6\n".format(ntype, self.kinetic))

    def create_STRU(self, path, cell, position, a, id, natom):
        # create .ion file for PROFESS
        with open(path+"/STRU", 'w') as f:
            f.write("ATOMIC_SPECIES\n")
            for ion in id:
                f.write("{0} {1} {2} blps\n".format(ion, self.mass[ion], self.pseudo[ion]))
            f.write("\nLATTICE_CONSTANT\n\
{0}  // add lattice constant\n\
\n\
LATTICE_VECTORS\n".format(a/self.bohr2a))
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 

            recordatom = 0
            f.write('\nATOMIC_POSITIONS\nDirect\n\n')
            for j, ion in enumerate(id):
                f.write('{0}\n0.0\n{1}\n'.format(ion, natom[j]))
                for k in range(natom[j]):
                    f.writelines('    ' + '%.12f' % each for each in position[k+recordatom])
                    f.write(' 1 1 1\n')
                    # f.write(' 0 0 0\n')
                recordatom += natom[j]

    def create_KPT(self, path, k):
        with open(path+"/KPT", 'w') as f:
            f.write("K_POINTS\n\
0\n\
Gamma\n\
{0} {0} {0} 0 0 0\n".format(k))

    def create_files(self):
        if not os.path.exists('log'): os.mkdir('log')
        with open('../jobs', 'a') as job:
            job.write('temp_path=$PWD\n')
            for i in range(self.Nstru):
                cell = self.cell[i]
                cell[:, -1] *= self.covera[i]
                path = self.structures[i]

                nbands = 0
                for j, ion in enumerate(self.id[i]):
                    nbands += self.natom[i][j] * self.zvan[ion]
                nbands = int(nbands/2 * 1.25) + 1

                os.mkdir(path)
                self.create_INPUT(path, len(self.id[i]), nbands)
                self.create_STRU(path, cell, self.position[i], self.a[i], self.id[i], self.natom[i])
                self.create_KPT(path, self.k[i])
                job.write("cd $temp_path\{0}/{1} && {2} > ../log/{3} && echo \"{4} done\"\n".format(self.sub_dir, path, self.abacus, path, path))

#                 with open('{0}/submit.sh'.format(path), 'w') as shell:
#                     shell.write("#!/bin/bash\n\
# #SBATCH -J {0}\n\
# #SBATCH -p cn_nl\n\
# #SBATCH -N 1 \n\
# #SBATCH -o relax.out\n\
# #SBATCH -e relax.err\n\
# #SBATCH --no-requeue\n\
# #SBATCH -A mhchen_g1 \n\
# #SBATCH --qos=mhchencnnl\n\
# #SBATCH -n 16\n\
# \n\
# source /home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/env.sh\n\
# export OMP_NUM_THREADS=1\n\
# mpirun -n 16 abacus > ../log/{0} && echo \"{0} done\"\n".format(self.structures[i]))

    # def calculate(self):
    #     # perform structure optimization
    #     # e_list = np.zeros((2,self.N * self.Nstru))
    #     cal = subprocess.Popen(args='parallel < jobs',
    #                            shell=True, universal_newlines=False)
    #     cal.wait()
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
    kernel_dir = "/home/dell/1_work/TKK_abacus/1_EFTKK_Al/5_results/"
    calculation_dir = {0:"scf", 1:"cell-relax", 2:"ofdft"}

    minNatom = 0
    maxNatom = 100
    calculation = "scf"
    esolver_type = 'ofdft'
    ecut = 0
    xc = 1
    # minNatom = int(input("Please enter the min number of atoms:"))
    # maxNatom = int(input("Please enter the max number of atoms:"))
    # calculation = calculation_dir[int(input("Please choose the calculatin (0:scf, 1:cell-relax, 2:ofdft):"))]
    # ecut = float(input("Please enter the ECUT(in eV):"))
    # xc = int(input("Please enter the type of XC (0:LDA, 1:PBE):"))

    kinetic = ""
    use_kernel = 0
    kernel_prefix = 0
    work_dir = os.getcwd()
    if esolver_type == "ofdft":
        kinetic = "wt"
        use_kernel = 0
        # kinetic = input("Please enter the name of KEDF:")
        # use_kernel = int(input("Do you want to use kernel?"))
        # if use_kernel: kernel_prefix = int(input("Please enter the prefix of kernel file (one of 1, 3, 4, 5, 8, 11):"))
    for sub_dir in ["0_LiMg", "1_MgAl", "2_LiAl", "3_LiMgAl", "4_single"]:
        total_dir = stru_dir + sub_dir
        if not os.path.exists(sub_dir): os.mkdir(sub_dir)
        os.chdir(sub_dir)
        abacus = ABACUS(minNatom, maxNatom, ecut, total_dir, calculation, esolver_type, kinetic, use_kernel, kernel_dir, kernel_prefix, xc)
        os.chdir(work_dir)
    cal = subprocess.Popen(args='parallel < jobs',
                            shell=True, universal_newlines=False)
    cal.wait()
    # profess = PROFESS(minNatom, maxNatom, ecut, stru_dir, kinetic, use_kernel, kernel_dir, kernel_prefix)
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
