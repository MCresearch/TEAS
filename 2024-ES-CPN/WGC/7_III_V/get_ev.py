import numpy as np
import os
import subprocess
#from scipy.optimize import curve_fit

a_dir = {"AlP":5.457501663522375,
         "AlAs":5.587728229478899, 
         "AlSb": 6.095042106018323,
         "GaP": 5.3201793898089935,
         "GaAs": 5.457367361384721,
         "GaSb": 5.943469048010328,
         "InP": 5.943469048010328,
         "InAs": 5.943469048010328,
         "InSb": 6.313283469726961,
}

gamma_dir = {"AlP":4.2,
         "AlAs":4.2,
         "AlSb": 4.2,
         "GaP": 4.2,
         "GaAs": 4.2,
         "GaSb": 4.2,
         "InP": 4.2,
         "InAs": 4.2,
         "InSb": 4.2,
}

class Abacus:
    def __init__(self, natom_in):
        self.pPROFESS = "PROFESS"
        self.bohr2a = 0.529177249
        self.natom = natom_in
        
        structure = ""
        for each in natom_in.keys():
            structure += each
        self.structures = [structure]
        self.a = [a_dir[structure]]
        # self.ry2ev = "13.6056923"
        # self.structures = ['ZB']
        # self.a = [5.457367361384721]
        self.covera = [1]
        self.cell = [np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])]
        self.position = [np.array([[0, 0, 0], [1/4, 1/4, 1/4]])]
        self.totnatom = [len(each) for each in self.position]
        self.ecut = 3200 # in Ry
        # self.k = [12, 12, 16, 16, 12, 16]
        # self.smearing = [0.0002] + [0.0074] * 5
        self.pseudo = {'Li':'li.gga.recpot.1', 'Mg':'mg.gga.recpot', 'Al':'al.gga.recpot', 'P':'p.gga.recpot', 'Ga':'ga.gga.recpot', 'As':'as.gga.recpot', 'In':'in.gga.recpot', 'Sb':'sb.gga.recpot'}
        # self.natom = {'Ga':1, 'As':1}
        
        self.id = []
        for each in natom_in.keys():
            for i in range(natom_in[each]):
                self.id.append(each)
        # self.id = ['Ga'] * self.natom['Ga'] + ['As'] * self.natom['As']
        self.zatom = {'Li':3, 'Mg':12, 'Al':13, 'P':15, 'Ga':31, 'As':33, 'In':49, 'Sb':51}
        self.N = 11
        self.full_pw = 1
        self.full_pw_dim = 1
        self.gamma = gamma_dir[structure]
        # self.create_files()
        # self.cal_e_v()

    def create_inpt(self, inptfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WGC\nPARA GAMM {1}\nEXCH   PBE\n'.format(
                self.ecut, self.gamma))
            
    def create_ion(self, ionfile, cell, position):
        # create .ion file for PROFESS
        with open(ionfile, 'w') as f:
            f.write('%BLOCK LATTICE_CART\n')
            for eachline in cell:
                f.writelines('%.12f' % each + '    ' for each in eachline)
                f.write('\n') 
            f.write('%END BLOCK LATTICE_CART\n')
            f.write('%BLOCK POSITIONS_FRAC\n')
            for i in range(len(self.id)):
                f.write(self.id[i])
                f.writelines('    ' + '%.12f' % each for each in position[i])
                f.write('\n')
            f.write('%END BLOCK POSITIONS_FRAC\n')
            f.write('%BLOCK SPECIES_POT\n')
            for i in range(len(self.id)):
                f.write('{0} {1}\n'.format(self.id[i], self.pseudo[self.id[i]]))
            f.write('%END BLOCK SPECIES_POT\n')

    def create_files(self):
        coef = np.linspace(0.9, 1.1, self.N)
        with open('jobs', 'a') as job:
            job.write('temp_path=$PWD\n')
            for i in range(len(self.structures)):
                for j in range(self.N):
                    cell = self.cell[i] * self.a[i] * coef[j]
                    cell[:, -1] *= self.covera[i]
                    path = self.structures[i] + str(j)
                    os.mkdir(path)
                    self.create_inpt(path+'/'+self.structures[i]+'.inpt')
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
        coef = np.linspace(0.9, 1.1, self.N)
        for i in range(len(self.structures)):
            v = abs(np.dot(self.cell[i][0], np.cross(self.cell[i][1], self.cell[i][2])))
            v *= self.a[i] ** 3 * self.covera[i] / self.totnatom[i]
            for j in range(self.N):
                v_list[self.N * i + j] = v * coef[j] ** 3
        with open(self.structures[0] + '.curve', 'w') as f:
            f.write('Volume per atom(Ang^3)\tEnergy per atom(eV)\n')
            f.writelines(each + '\t' for each in self.structures)
            f.write('\n')
            for i in range(self.N * len(self.structures)):
                f.write('%.12e\t%.12e\n' % (v_list[i], e_list[i]))

if __name__ == '__main__':
    stru = []
    for III in ["Al", "Ga", "In"]:
        for IV in ["P", "As", "Sb"]:
            stru.append({III:1, IV:1})
    for each in stru:
        abacus = Abacus(each)
        abacus.create_files()
    cal = subprocess.Popen(args='parallel < jobs',
                        shell=True, universal_newlines=False)
    cal.wait()
    
    for each in stru:
        abacus = Abacus(each)
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