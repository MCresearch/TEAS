import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt

id = 'Al'
ecut = 800
pseudo = 'al_HC.lda.recpot'
PROFESS = '/home/dell/2_software/g_pPROFESS/pPROFESS'


def create_inpt(inptfile, ecut):
    # create .inpt file for PROFESS
    with open(inptfile, 'w') as f:
        f.write('ECUT	{0}\nKINE	WT\nMINI CEL\nMETH CEL 3\nMETH ION CON\n'.format(ecut))

def create_ion(ionfile, id, cell, position, pseudo):
    # create .ion file for PROFESS
    with open(ionfile, 'w') as f:
        f.write('%BLOCK LATTICE_ABC\n')
        for eachline in cell:
            f.writelines('%.10f' % each + '    ' for each in eachline)
            f.write('\n')
        f.write('%END BLOCK LATTICE_ABC\n')
        f.write('%BLOCK POSITIONS_FRAC\n')
        for eachline in position:
            f.write(id)
            f.writelines('    ' + '%.10f' % each for each in eachline)
            f.write('\n')
        f.write('%END BLOCK POSITIONS_FRAC\n')
        f.write('%BLOCK ION_OPTIMIZATION\n')
        for i in range(len(position)):
            f.write('0 0 1\n')
        f.write('%END BLOCK ION_OPTIMIZATION\n')
        f.write(
            '%BLOCK SPECIES_POT\n{0} {1}\n%END BLOCK SPECIES_POT\n'.format(id, pseudo))

def create_files(positions):
    with open('jobs', 'w') as job:
        for i, position in enumerate(positions):
            path = id + str(i)
            os.mkdir(path)
            create_inpt(path+'/'+id+'.inpt',ecut)
            if len(position) == 20:
                create_ion(path+'/'+id+'.ion', id, cell1, position, pseudo)
            elif len(position) == 22:
                create_ion(path+'/'+id+'.ion', id, cell2, position, pseudo)
            elif len(position) == 21:
                create_ion(path+'/'+id+'.ion', id, cell3, position, pseudo)
            job.write(PROFESS +' '+path+'/'+id+'\n')

def cal_e_v(positions):
    # calculate e_v curve and return energy list
    n = len(positions)
    create_files(positions)
    elist = np.zeros(n)
    cal = subprocess.Popen(args='parallel < jobs',
                            shell=True, universal_newlines=False)
    cal.wait()
    with open('stacking', 'w') as f:
        for i in range(n):
            outfile = './{0}{1}/{0}.out'.format(id, i)
            try:
                get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % outfile,
                                                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
                total_e = float(get_total_e.stdout.read())*1.602e-16/1e-20/(2.8064544120**2*3**0.5/2)  # maybe wrong
                #if i < 11:
                #    total_e *= 8
                #else:
                #    total_e *= 12
            except:
                total_e = 1e5
            elist[i] = total_e
            f.write(str(total_e)+'\n')
    plt.figure()
    plt.plot(np.linspace(0,1,11), (elist[1:12]-elist[1])/2)
    plt.plot(np.linspace(1,2,11), (elist[12:]-2*elist[0]/21-elist[1])/2)
    plt.savefig('stacking')

cell1 = np.array([[2.8057997077482204, 2.8057997077482204, 45.81851736288869, 90, 90, 60]])
cell2 = np.array([[2.8057997077482204, 2.8057997077482204, 45.81851736288869/20*22, 90, 90, 60]])
cell3 = np.array([[2.8057997077482204, 2.8057997077482204, 45.81851736288869/20*21, 90, 90, 60]])

A = np.array([0, 0, 0])
B = np.array([1/3, 1/3, 0])
C = np.array([2/3, 2/3, 0])
positions = [np.vstack([A,B,C,A,B,C,A,B,C,A,B,C,A,B,C,A,B,C,A,B,C])]
positions[0][:,-1] += np.linspace(0,20/21,21)
position = np.vstack([A,B,C,A,B,C,A,B,C,A,C,B,A,C,B,A,C,B,A,C])
position[:,-1] += np.linspace(0, 19/20, 20)
for i in range(11):
    positions.append(position.copy())
    position[5:14,0:2] += 1/30
position = np.vstack([A,B,C,A,B,A,B,C,A,B,C,B,A,C,B,A,B,A,C,B,A,C])
position[:,-1] += np.linspace(0, 21/22, 22)
for i in range(11):
    positions.append(position.copy())
    position[6:15,0:2] += 1/30
#create_files(positions)
cal_e_v(positions)

