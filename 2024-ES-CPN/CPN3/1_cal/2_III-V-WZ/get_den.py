import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit

a_dir = {
    "AlP": 1.00499318959644856619,
    "AlAs": 1.00910657092604294505,
    "AlSb": 1.00296719422196956018,
    "GaP": 1.00141777494820782834,
    "GaAs": 1.00533848595329944331,
    "GaSb": 1.00465765524053662894,
    "InP": 1.00442538575014350677,
    "InAs": 1.00868172578763637404,
    "InSb": 1.00400591789164361778,
}

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
        self.abacus = "/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/abacus-mlkedf-r_min/abacus-develop/build/abacus"
        # self.abacus = "/home/xianyuer/data/1_sunliang/1_work/0_ml_kedf/abacus_merge/abacus-develop/build_gnu_openmp/abacus"
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
        
        self.id = []
        self.read_stru(file_in)
        self.a[0] = a_dir[structure]

        self.N = 1
        self.full_pw = 1
        self.full_pw_dim = 1

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

    def create_INPUT(self, path):
        def write_list(f, parameter, value):
            temp = parameter + "    "
            for each in value:
                temp += '{0}_'.format(each)
            f.write(temp[:-1] + "\n")

        work_dir = os.getcwd()
        folder = work_dir.split('/')[-4]

        temp = folder.split('-')[1:]
        setting = []
        setting_detail = []
        for each in temp:
            detail = []
            n_tail = 0
            for each_char in each:
                if each_char.isdigit():
                    n_tail -= 1
                    detail.append(int(each_char))
            setting.append(each[:n_tail])
            setting_detail.append(detail)
            
        temp = folder.split('-')[0].split('g')[1]
        chi_xi_dir = {4:0.6, 3:0.6, 2:0.6, 1.5:0.8, 1:1, 0.75:1.5, 0.5:3, 0.35:6, 0.25:12}
        kernel_scaling = [float(_) for _ in temp.split('_')]
        chi_p = [0.2]
        chi_q = [0.1]
        chi_xi = []
        chi_pnl = []
        chi_qnl = []
        for each in kernel_scaling:
            chi_xi.append(chi_xi_dir[each])
            chi_pnl.append(0.2)
            chi_qnl.append(0.1)
        nkernel = len(kernel_scaling)
        
        # nlayer = folder.split('_')[1]
        # nnode = folder.split('_')[2].split('-')[0]
        nlayer = 3
        nnodes = [100, 100, 100, 100] * 5
        nnode = nnodes[int(folder[0])]

        with open(path+"/INPUT", 'w') as f:
            f.write("INPUT_PARAMETERS\n\
#Parameters (1.General)\n\
suffix			oftest\n\
calculation     scf\n\
esolver_type    ofdft\n\
ntype			2\n\
symmetry		1\n\
out_chg         1\n\
# pseudo_dir      /home/xianyuer/data/1_sunliang/2_software/PP_ORB/BLPSLibrary-master/GGA/real\n\
pseudo_dir      /home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/PROFESS3.0/BLPSLibrary-master/GGA/real/\n\
pseudo_rcut     16\n\
\n\
#Parameters (2.Iteration)\n\
ecutwfc			{0}\n\
scf_nmax        200\n\
dft_functional  XC_GGA_X_PBE+XC_GGA_C_PBE\n\
\n\
\n\
#OFDFT\n\
#of_kinetic     tf\n\
#of_kinetic     vw\n\
of_kinetic      ml\n\
of_method       tn\n\
of_conv         energy\n\
of_tole         2e-6\n\
of_full_pw      {1}\n\
of_full_pw_dim  {2}\n\
of_ml_device    gpu\n\
of_ml_feg       3\n\
\n\
#Parameters (3.Basis)\n\
basis_type		pw\n\
of_ml_nkernel   {3}\n".format(
                self.ecut, self.full_pw, self.full_pw_dim, nkernel))

            write_list(f, "of_ml_chi_xi", chi_xi)
            write_list(f, "of_ml_chi_p", chi_p)
            write_list(f, "of_ml_chi_q", chi_q)
            write_list(f, "of_ml_chi_pnl", chi_pnl)
            write_list(f, "of_ml_chi_qnl", chi_qnl)
            write_list(f, "of_ml_kernel", [1] * nkernel)
            write_list(f, "of_ml_yukawa_alpha", [1] * nkernel)
            
            write_list(f, "of_ml_kernel_scaling", kernel_scaling)

            for i, each in enumerate(setting):
                write_list(f, 'of_ml_{0}'.format(each), setting_detail[i])

            f.write('of_ml_nnode    {0}\n'.format(nnode))
            f.write('of_ml_nlayer   {0}\n'.format(nlayer))

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
        work_dir = os.getcwd()
        coef = np.array([1.0])
        with open('jobs', 'a') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i].copy()
                    cell[:, -1] *= self.covera[i]
                    a = self.a[i] * coef[j]
                    path = self.structures[i]
                    # path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_INPUT(path)
                    self.create_STRU(path, cell, self.position[i], a)
                    self.create_KPT(path)
                    os.chdir(path)
                    cal = subprocess.Popen(args="cp ../../../../model/net60000.pt ./net.pt && " + self.abacus + ' > log\n', shell=True, universal_newlines=False)
                    cal.wait()
                    os.chdir(work_dir)
                    job.write("cd $temp_path\{0} && cp ../../../../model/net60000.pt ./net.pt && ".format(path) + self.abacus + ' > log\n')



    def cal_e_v(self):
        # calculate e_v curve and return energy list
        e_list = np.zeros(self.N * len(self.structures))
        # cal = subprocess.Popen(args='parallel < jobs',
        #                        shell=True, universal_newlines=False)
        # cal.wait()
        for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}/OUT.oftest/running_scf.log'.format(self.structures[i])
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
        coef = np.array([1.0])
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
    stru_dir = "/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/1_work/9_mlkedf/0_generate_data/2_ks-pbe/0_data_set/b_wurtzite/0_structures/"
    # stru_dir = "/home/xianyuer/data/1_sunliang/1_work/0_ml_kedf/1_test/0_generate_data/2_ks-pbe/0_data_set/b_wurtzite/0_structures/"
    for III in ["Al", "Ga", "In"]:
        for V in ["P", "As", "Sb"]:
            stru.append(III + V)
    for each in stru:
        abacus = Abacus(stru_dir + files[each], each)
        abacus.create_files()
    # cal = subprocess.Popen(args='parallel < jobs',
    #                     shell=True, universal_newlines=False)
    # cal.wait()
    
    # for each in stru:
    #     abacus = Abacus(stru_dir + files[each], each)
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
