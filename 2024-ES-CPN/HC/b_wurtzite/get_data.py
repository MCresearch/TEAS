import numpy as np
import os
import subprocess
from scipy.optimize import curve_fit

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

lambda_dir = {"AlP":0.012,
         "AlAs":0.0125,
         "AlSb": 0.012,
         "GaP": 0.01,
         "GaAs": 0.013,
         "GaSb": 0.01,
         "InP": 0.012,
         "InAs": 0.0142,
         "InSb": 0.012,
}

beta_dir = {"AlP":0.845,
         "AlAs":0.825,
         "AlSb": 0.750,
         "GaP": 0.791,
         "GaAs": 0.783,
         "GaSb": 0.720,
         "InP": 0.885,
         "InAs": 0.875,
         "InSb": 0.810,
}

class Abacus:
    def __init__(self, file_in, struc_in):
        print(struc_in)
        self.pPROFESS = "PROFESS"
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
        self.ecut = 3200 # in Ry
        # self.k = [12, 12, 16, 16, 12, 16]
        # self.smearing = [0.0002] + [0.0074] * 5
        self.pseudo = {'Li':'li.gga.recpot.1', 'Mg':'mg.gga.recpot', 'Al':'al.gga.recpot', 'P':'p.gga.recpot', 'Ga':'ga.gga.recpot', 'As':'as.gga.recpot', 'In':'in.gga.recpot', 'Sb':'sb.gga.recpot'}
        
        self.id = []
        self.read_stru(file_in)

        self.v = np.zeros(1)
        for i in range(1):
            self.v[i] = abs(np.dot(self.cell[i][0], np.cross(self.cell[i][1], self.cell[i][2])))
            self.v[i] *= self.covera[i] / self.totnatom[i]

        self.N = 11
        self.full_pw = 1
        self.full_pw_dim = 1
        self.lam = lambda_dir[structure]
        self.beta = beta_dir[structure]

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

    def create_inpt(self, inptfile, beta, hcla):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	HC\nNNDI 1\nPARA BETA {1}\nPARA HCLA {2}\nEXCH   PBE\n'.format(
                self.ecut, beta, hcla))

    def create_ion(self, ionfile, cell, position):
        # create .ion file for PROFESS
        with open(ionfile, 'w') as f:
            f.write('%BLOCK LATTICE_CART\n')
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('%END BLOCK LATTICE_CART\n')
            f.write('%BLOCK POSITIONS_FRAC\n')
            posi_index = 0
            for i in range(len(self.id)):
                for j in range(self.natom[self.id[i]]):
                    f.write(self.id[i])
                    f.writelines('    ' + '%.12f' % each for each in position[posi_index])
                    f.write('\n')
                    posi_index += 1
            f.write('%END BLOCK POSITIONS_FRAC\n')
            f.write('%BLOCK SPECIES_POT\n')
            for i in range(len(self.id)):
                f.write('{0} {1}\n'.format(self.id[i], self.pseudo[self.id[i]]))
            f.write('%END BLOCK SPECIES_POT\n')

    def create_files(self):
        coef = np.linspace(0.99, 1.01, self.N)
        with open('jobs', 'a') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i] * self.a[i] * coef[j]
                    cell[:, -1] *= self.covera[i]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_inpt(path+'/'+self.structures[i]+'.inpt', self.beta, self.lam)
                    self.create_ion(path+'/'+self.structures[i]+'.ion', cell, self.position[i])
                    job.write(self.pPROFESS + ' ' + path + '/' + self.structures[i] + ' > '+ path + '/log\n')

    def cal_e_v(self):
        # calculate e_v curve and return energy list
        e_list = np.zeros(self.N * len(self.structures))
        # cal = subprocess.Popen(args='parallel < jobs',
        #                        shell=True, universal_newlines=False)
        # cal.wait()
        for i in range(len(self.structures)):
            for j in range(self.N):
                outfile = './{0}{1}/{2}.out'.format(self.structures[i], j, self.structures[i])
                try:
                    get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
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
        coef = np.linspace(0.99, 1.01, self.N)
        v0 = np.zeros(len(self.structures))
        e0 = np.zeros_like(v0)
        B = np.zeros_like(v0)
        for i in range(len(self.structures)):
            for j in range(self.N):
                v_list[self.N * i + j] = self.v[i] * (coef[j] * self.a[i]) ** 3
            e0[i], v0[i], B[i] = abacus.fit(v_list[i*11:(i+1)*11], e_list[i*11:(i+1)*11], v_list[i*11], v_list[i*11+10])
        with open('{0}.curve'.format(self.structures[0]), 'w') as f:
            f.write('Volume per atom(Ang^3)\tEnergy per atom(eV)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(self.N * len(self.structures)):
                f.write('%.12e\t%.12e\n' % (v_list[i], e_list[i]))
        with open('data', 'a') as f:
            if os.path.getsize('data') == 0:
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
    stru = []
    # stru_dir = "/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/1_work/9_mlkedf/0_generate_data/2_ks-pbe/0_data_set/b_wurtzite/0_structures/"
    stru_dir = "/home/xianyuer/data/1_sunliang/1_work/0_ml_kedf/1_test/0_generate_data/2_ks-pbe/0_data_set/b_wurtzite/0_structures/"
    for III in ["Al", "Ga", "In"]:
        for V in ["P", "As", "Sb"]:
            stru.append(III + V)
    for each in stru:
        input = '../0_e-v_curve/{}.curve'.format(each)
        curve = np.loadtxt(input, skiprows=2).T
        abacus = Abacus(stru_dir + files[each], each)
        for i in range(1):
            e0, v0, B = abacus.fit(curve[0][i*11:(i+1)*11], curve[1][i*11:(i+1)*11], curve[0][i*11], curve[0][i*11+10])
            abacus.a[i] = (v0/abacus.v[i]) ** (1/3)
            print(abacus.a[i])
        abacus.create_files()
    cal = subprocess.Popen(args='parallel < jobs',
                        shell=True, universal_newlines=False)
    cal.wait()
    
    for each in stru:
        input = '../0_e-v_curve/{}.curve'.format(each)
        curve = np.loadtxt(input, skiprows=2).T
        abacus = Abacus(stru_dir + files[each], each)
        for i in range(1):
            e0, v0, B = abacus.fit(curve[0][i*11:(i+1)*11], curve[1][i*11:(i+1)*11], curve[0][i*11], curve[0][i*11+10])
            abacus.a[i] = (v0/abacus.v[i]) ** (1/3)
            print(abacus.a[i])
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
