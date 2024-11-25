import numpy as np
import os
import subprocess
from scipy.optimize import curve_fit


class Profess:
    def __init__(self, settings):
        self.pPROFESS = "/home/mhchen_pkuhpc/mhchen_coe/lustre2/1_sunliang/2_software/PROFESS3.0_modified/pPROFESS"
        # self.pPROFESS = 'pPROFESS'
        self.create_files(settings)

    def create_inpt(self, inptfile, ecut, kernelfile):
        # create .inpt file for PROFESS
        with open(inptfile, 'w') as f:
            f.write('ECUT	{0}\nKINE	WT\nKERN {1}\nMAXI DEN 50\n'.format(
                ecut, kernelfile))

    def create_ion(self, ionfile, a, pseudo, c=0, structrue='cd'):
        # create .ion file for PROFESS
        if c == 0:
            c = a
        if structrue == 'cd':
            with open(ionfile, 'w') as f:
                f.write(
                    '%BLOCK LATTICE_CART\n0 {0} {0}\n{0} 0 {0}\n{0} {0} 0\n%END BLOCK LATTICE_CART\n'.format(a/2))
                f.write('%BLOCK POSITIONS_FRAC\n')
                f.write('Si 0 0 0\n')
                f.write('Si 0.25 0.25 0.25\n')
                f.write('%END BLOCK POSITIONS_FRAC\n')
                f.write(
                    '%BLOCK SPECIES_POT\nSi {0}\n%END BLOCK SPECIES_POT\n'.format(pseudo))
        elif structrue == 'beta-tin':
            with open(ionfile, 'w') as f:
                f.write(
                    '%BLOCK LATTICE_CART\n-{0} {0} {1}\n{0} -{0} {1}\n{0} {0} -{1}\n%END BLOCK LATTICE_CART\n'.format(a/2, c/2))
                f.write('%BLOCK POSITIONS_FRAC\n')
                f.write('Si 0 0 0\n')
                f.write('Si 0.75 0.25 0.5\n')
                f.write('%END BLOCK POSITIONS_FRAC\n')
                f.write(
                    '%BLOCK SPECIES_POT\nSi {0}\n%END BLOCK SPECIES_POT\n'.format(pseudo))

    def create_files(self, settings):
        with open('jobs', 'w') as job:
            job.write(self.pPROFESS + ' {0}/{1}\n'.format(settings.structures[1], settings.name))
            for i, a in enumerate(settings.aC_list):
                path = 'cd-' + settings.name + str(i)
                os.mkdir(path)
                self.create_inpt(path+'/'+settings.name+'.inpt', settings.ecut, settings.kernelfile)
                self.create_ion(path+'/'+settings.name +
                                '.ion', a, settings.pseudo)
                job.write(self.pPROFESS +
                          ' '+path+'/'+settings.name+'\n')
            for i, a in enumerate(settings.aB_list):
                path = 'beta-tin-' + settings.name + str(i)
                os.mkdir(path)
                self.create_inpt(path+'/'+settings.name+'.inpt', settings.ecut, settings.kernelfile)
                self.create_ion(path+'/'+settings.name +
                                '.ion', a, settings.pseudo, settings.cB_list[i], 'beta-tin')
                job.write(self.pPROFESS +
                          ' '+path+'/'+settings.name+'\n')

    def cal_e_v(self, n, name, structures=[]):
        # calculate e_v curve and return energy list
        N = 2*n + len(structures)
        e_list = np.zeros((2, N))
        cal = subprocess.Popen(args='parallel < jobs',
                               shell=True, universal_newlines=False)
        cal.wait()
        for each in (structures[0], structures[2]):
            cal = subprocess.Popen(args='mpirun -n 23 {0} {1}/{2}'.format(self.pPROFESS, each, name), shell=True, universal_newlines=False)
            cal.wait()
        for i in range(N):
            if i < n:
                outfile = './cd-{0}{1}/{0}.out'.format(name, i)
            elif i < 2*n and i >= n:
                outfile = './beta-tin-{0}{1}/{0}.out'.format(name, i - n)
            elif i >= 2*n:
                outfile = './{0}/{1}.out'.format(structures[i-2*n], name)
            try:
                get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read())  # maybe wrong
                # get_kinetic_e = subprocess.Popen(args="grep 'TOTAL KINETIC ENERGY' %s|tr -s ' '|cut -d ' ' -f6" % outfile,
                #                                  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                # kinetic_e = float(get_kinetic_e.stdout.read())
                get_hartree_e = subprocess.Popen(args="grep 'Coulombic Energy' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                 stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                hartree_e = float(get_hartree_e.stdout.read())
            except:
                total_e = 1e5
                # kinetic_e = 1e5
                hartree_e = 1e5
            e_list[0,i] = total_e
            # e_list[1,i] = kinetic_e
            e_list[1,i] = hartree_e
        return e_list

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
