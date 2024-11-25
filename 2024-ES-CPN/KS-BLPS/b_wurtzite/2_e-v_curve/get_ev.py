import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit

# a_dir = {"AlP":5.457501663522375,
#          "AlAs":5.587728229478899, 
#          "AlSb": 6.095042106018323,
#          "GaP": 5.3201793898089935,
#          "GaAs": 5.457367361384721,
#          "GaSb": 5.943469048010328,
#          "InP": 5.943469048010328,
#          "InAs": 5.943469048010328,
#          "InSb": 6.313283469726961,
# }

files = {
    "AlP": "mp-8880-AlP.poscar",
    "AlAs": "mp-8881-AlAs.poscar",
    "AlSb": "mp-1018100-AlSb.poscar",
    "GaP": "mp-8882-GaP.poscar",
    "GaAs": "mp-8883-GaAs.poscar",
    "GaSb": "mp-1018059-GaSb.poscar",
    "InP": "mp-966800-InP.poscar",
    "InAs": "mp-1007652-InAs.poscar",
    "InSb": "mp-1007661-InSb.poscar",
}

class Abacus:
    def __init__(self, file_in, struc_in):
        print(struc_in)
        self.abacus = "/home/xianyuer/data/1_sunliang/1_work/0_ml_kedf/abacus_merge/abacus-develop/build_gnu_openmp/abacus"
        self.bohr2a = 0.529177249
        self.natom = dict()
        
        self.elements = ["Al", "Mg", "Li", "Si", "Ga", "In", "P", "As", "Sb"]
        self.mass = {"Al":26.981, "Mg":24.305, "Li":6.941, "Si":28.085, "Ga":69.723, "In":114.818, "P":30.974, "As":74.922, "Sb":121.760}
        self.zvan = {"Al":3, "Mg":2, "Li":1, "Si":4, "Ga":3, "In":3, "P":5, "As":5, "Sb":5}
        structure = struc_in
        self.structures = [structure]
        self.id = []
        self.a = []
        self.cell = []
        self.position = []
        self.covera = [1]
        self.totnatom = []
        self.ecut = 60 # in Ry
        # self.k = [12, 12, 16, 16, 12, 16]
        # self.smearing = [0.0002] + [0.0074] * 5
        self.pseudo = {'Li':'li.gga.psp.1', 'Mg':'mg.gga.psp', 'Al':'al.gga.psp', 'P':'p.gga.psp', 'Ga':'ga.gga.psp', 'As':'as.gga.psp', 'In':'in.gga.psp', 'Sb':'sb.gga.psp'}
        # self.natom = {'Ga':1, 'As':1}
        
        self.id = []
        self.read_stru(file_in)

        # self.id = ['Ga'] * self.natom['Ga'] + ['As'] * self.natom['As']
        # self.zatom = {'Li':3, 'Mg':12, 'Al':13, 'P':15, 'Ga':31, 'As':33, 'In':49, 'Sb':51}
        self.N = 11
        self.full_pw = 1
        self.full_pw_dim = 1
        self.kspacing = 0.05

        # self.create_files()
        # self.cal_e_v()
        
    def read_stru(self, file):
        with open(file, 'r') as stru_file:
            title = stru_file.readline()
            for element in self.elements:
                title = title.replace(element, '')
            totnatom_ = sum([int(_) for _ in title.split()])
            self.totnatom.append(totnatom_)
            # if totnatom > self.maxNatom or totnatom <= self.minNatom:
            #     return False
            # else:
            self.a.append(float(stru_file.readline()))
            self.cell.append(
                np.array([[float(_) for _ in stru_file.readline().split()] for j in range(3)])
            )
            self.id = stru_file.readline().split()
            print(self.id)
            natom_ = np.array(stru_file.readline().split(), dtype=np.int64)
            self.natom = dict(zip(self.id, natom_))
            print(self.natom)
            stru_file.readline()
            self.position.append(
                np.array([[float(_) for _ in stru_file.readline().split()[:3]] for j in range(totnatom_)])
            )
            return True

    def create_INPUT(self, path, nbands):
        # create INPUT file for ABACUS
        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			oftest\n\
calculation     scf\n\
symmetry		1\n\
#symmetry_prec   2e-4\n\
pseudo_dir      /home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/PROFESS3.0/BLPSLibrary-master/GGA/real\n\
# pseudo_dir      /home/xianyuer/data/1_sunliang/2_software/PP_ORB/BLPSLibrary-master/GGA/real\n\
pseudo_rcut     16\n\
nbands          {0}\n\
out_stru        1\n\
out_chg         1\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{1}\n\
scf_nmax        100\n\
\n\
#Parameters (3.Basis)\n\
basis_type		pw\n\
kspacing        {2}\n\
\n\
#Parameters (4.Smearing)\n\
smearing_method	gauss\n\
smearing_sigma	0.0074\n\
of_full_pw      1\n\
of_full_pw_dim  1\n\
\n\
#Parameters (5.Mixing)\n\
mixing_type		pulay\n\
mixing_beta		0.7\n".format(
                nbands, self.ecut, self.kspacing))

            f.write("\n\
dft_functional  XC_GGA_X_PBE+XC_GGA_C_PBE\n")

    def create_STRU(self, path, cell, position, a):
        # create .ion file for PROFESS
        with open(path+"/STRU", 'w') as f:
            f.write("ATOMIC_SPECIES\n\
{0} {1} {2} blps\n\
{3} {4} {5} blps\n\
\n\
LATTICE_CONSTANT\n\
{6}  // add lattice constant\n\
\n\
LATTICE_VECTORS\n".format(self.id[0], self.mass[self.id[0]], self.pseudo[self.id[0]], self.id[1], self.mass[self.id[1]], self.pseudo[self.id[1]], a/self.bohr2a))
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('\nATOMIC_POSITIONS\nDirect\n\n{0}\n0.0\n{1}\n'.format(self.id[0], self.natom[self.id[0]]))
            for i in range(self.natom[self.id[0]]):
                f.writelines('    ' + '%.12f' % each for each in position[i])
                f.write(' 1 1 1\n')
            f.write('{0}\n0.0\n{1}\n'.format(self.id[1], self.natom[self.id[1]]))
            for i in range(self.natom[self.id[0]],self.natom[self.id[0]]+self.natom[self.id[1]]):
                f.writelines('    ' + '%.12f' % each for each in position[i])
                f.write(' 1 1 1\n')

    def create_KPT(self, path):
        with open(path+"/KPT", 'w') as f:
            f.write("K_POINTS\n\
0\n\
Gamma\n\
1 1 1 0 0 0\n")

    def create_files(self):
        if not os.path.exists('log'): os.mkdir('log')
        coef = np.linspace(0.9, 1.1, self.N)
        with open('jobs', 'a') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                nbands = 0
                for j, ion in enumerate(self.id):
                    nbands += self.natom[ion] * self.zvan[ion]
                nbands = int(nbands/2 * 1.25) + 1
                for j in range(self.N):
                    cell = self.cell[i].copy()
                    cell[:, -1] *= self.covera[i]
                    a = self.a[i] * coef[j]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_INPUT(path, nbands)
                    self.create_STRU(path, cell, self.position[i], a)
                    self.create_KPT(path)
                    job.write("cd $temp_path\{0} && ".format(path) + self.abacus + ' > log\n')

                    with open('{0}/submit.sh'.format(path), 'w') as shell:
                        shell.write("#!/bin/bash\n\
#SBATCH -J {0}\n\
#SBATCH -p cn_nl\n\
#SBATCH -N 1 \n\
#SBATCH -o relax.out\n\
#SBATCH -e relax.err\n\
#SBATCH --no-requeue\n\
#SBATCH -A mhchen_g1 \n\
#SBATCH --qos=mhchencnnl\n\
#SBATCH -n 8\n\
\n\
source /lustre3/mhchen_pkuhpc/mhchen_coe/4_dengzichao/env.sh\n\
export OMP_NUM_THREADS=1\n\
mpirun -n 8 /home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/abacus-3.7/abacus-develop/build/abacus > ../log/{0} && echo \"{0} done\"\n".format(self.structures[i] + str(j)))



    def cal_e_v(self):
        # calculate e_v curve and return energy list
        e_list = np.zeros(self.N * len(self.structures))
        # cal = subprocess.Popen(args='parallel < jobs',
        #                        shell=True, universal_newlines=False)
        # cal.wait()
        for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}{1}/OUT.oftest/running_scf.log'.format(self.structures[i], j)
                try:
                    get_total_e = subprocess.Popen(args="grep '!FINAL_ETOT_IS' %s|cut -d ' ' -f3" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                    total_e = float(get_total_e.stdout.read())/self.totnatom[i]  # maybe wrong
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
            v *= self.a[i] ** 3 * self.covera[i] / self.totnatom[i]
            for j in range(self.N):
                v_list[self.N * i + j] = v * coef[j] ** 3
        with open('{0}.curve'.format(self.structures[0]), 'w') as f:
            f.write('Volume per atom(Ang^3)\tEnergy per atom(eV)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(self.N * len(self.structures)):
                f.write('%.12e\t%.12e\n' % (v_list[i], e_list[i]))

if __name__ == '__main__':
    stru = []
    stru_dir = "/home/xianyuer/data/1_sunliang/1_work/0_ml_kedf/1_test/0_generate_data/2_ks-pbe/0_data_set/b_wurtzite/0_structures/"
    for III in ["Al", "Ga", "In"]:
        for V in ["P", "As", "Sb"]:
            stru.append(III + V)
    # for each in stru:
    #     abacus = Abacus(stru_dir + files[each], each)
    #     abacus.create_files()
    # cal = subprocess.Popen(args='parallel < jobs',
    #                     shell=True, universal_newlines=False)
    # cal.wait()
    
    for each in stru:
        abacus = Abacus(stru_dir + files[each], each)
        abacus.cal_e_v()
    
    
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
