import numpy as np
# import os
import subprocess
# import matplotlib.pyplot as plt

id = 'Al'
ecut = 900
pseudo = 'al.lda.recpot'
PROFESS = '/home/dell/2_software/g_pPROFESS/pPROFESS'


# def create_inpt(inptfile, ecut):
#     # create .inpt file for PROFESS
#     with open(inptfile, 'w') as f:
#         f.write('ECUT	{0}\nKINE	WT\nKERN PROFESS_KERNEL.dat\nMINI CEL\nMETH CEL 3\nMETH ION NON\n'.format(ecut))

# def create_ion(ionfile, id, cell, position, pseudo):
#     # create .ion file for PROFESS
#     with open(ionfile, 'w') as f:
#         f.write('%BLOCK LATTICE_ABC\n')
#         for eachline in cell:
#             f.writelines('%.10f' % each + '    ' for each in eachline)
#             f.write('\n')
#         f.write('%END BLOCK LATTICE_ABC\n')
#         f.write('%BLOCK POSITIONS_FRAC\n')
#         for eachline in position:
#             f.write(id)
#             f.writelines('    ' + '%.10f' % each for each in eachline)
#             f.write('\n')
#         f.write('%END BLOCK POSITIONS_FRAC\n')
#         f.write(
#             '%BLOCK SPECIES_POT\n{0} {1}\n%END BLOCK SPECIES_POT\n'.format(id, pseudo))

paths = ['fcc111', 'fcc100', 'fcc110', 'perfectfcc']
nums = [7,7,9]
s = [(3.968 / 2**0.5) ** 2 * 3 ** 0.5 / 2 * 1e-20, (3.968 / 2**0.5) ** 2 * 1e-20,  3.968 * 3.968/2**0.5 * 1e-20]

def create_files():
    with open('jobs', 'w') as job:
        for path in paths:
            job.write(PROFESS + ' '+path+'/al\n')


def cal_e_v():
    # calculate e_v curve and return energy list
    n = 3
    create_files()
    cal = subprocess.Popen(args='parallel < jobs',
                           shell=True, universal_newlines=False)
    cal.wait()
    get_efcc = subprocess.Popen(args="grep 'TOTAL ENERGY' ./perfectfcc/al.out|tr -s ' ' |cut -d ' ' -f5",  stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    efcc = float(get_efcc.stdout.read())
    
    with open('surface', 'w') as f:
        for i in range(n):
            outfile = './{0}/al.out'.format(paths[i])
            try:
                get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read())  # maybe wrong
                print(total_e)
                surface_e = (total_e - nums[i] * efcc) * 1.602e-16 / s[i] / 2
                print(surface_e)
            except:
                surface_e = 1e5
            f.write(str(surface_e)+'\n')


cal_e_v()
