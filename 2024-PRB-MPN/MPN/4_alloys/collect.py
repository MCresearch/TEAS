import numpy as np
import os
import subprocess

Bohr2A = 0.529177249

def fetch_stru_profess(stru_file):
    if not os.path.exists(stru_file):
        print('WARNING: {0} doesnot exist.'.format(stru_file))
        return 100, 100, dict()
    with open(stru_file) as f:
        f.readline()
        tempCell = np.array([[float(_) for _ in f.readline().split()] for j in range(3)])
        volume = abs(np.dot(tempCell[0], np.cross(tempCell[1], tempCell[2])))
        f.readline()
        f.readline()

        localAtoms = dict()
        while True:
            temp = f.readline()
            if temp[0] == "%":
                break
            else:
                id = temp.split()[0]
                if id not in localAtoms.keys():
                    localAtoms[id] = 1
                else:
                    localAtoms[id] += 1

        totalAtom = sum(localAtoms.values())
        print("fetch {0} done.\n".format(stru_file))
        return totalAtom, volume/totalAtom, localAtoms

def fetch_out_profess(out_file):
    try:
        get_total_e = subprocess.Popen(args="grep 'TOTAL ENERGY' %s|tr -s ' '|cut -d ' ' -f5" % out_file,
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        total_e = float(get_total_e.stdout.read())
        print("fetch {0} done.\n".format(out_file))
        return total_e
    except:
        print('WARNING: check {0}, total energy is not found!'.format(out_file))
        return 1e5

# def fetch_stru_abacus(stru_file):
#     if not os.path.exists(stru_file):
#         print('WARNING: {0} doesnot exist.'.format(stru_file))
#         return 100, 100, dict()
#     with open(stru_file) as f:
#         localAtoms = dict()
#         f.readline()
#         while True:
#             temp = f.readline()
#             if len(temp) == 1:
#                 break
#             else:
#                 localAtoms[temp.split()[0]] = 0
        
#         f.readline()
#         a = float(f.readline()) * Bohr2A
#         f.readline()
#         f.readline()

#         tempCell = np.array([[float(_) for _ in f.readline().split()[:3]] for j in range(3)])
#         volume = abs(np.dot(tempCell[0], np.cross(tempCell[1], tempCell[2])))
#         for _ in range(4): f.readline()
#         while True:
#             temp = f.readline()
#             if not temp: break
#             temp = temp.split()
#             if len(temp) > 0:
#                 if temp[0] in localAtoms.keys():
#                     element = temp[0]
#                     f.readline()
#                     localAtoms[element] = float(f.readline().split()[0])

#         totalAtom = sum(localAtoms.values())
#         print("fetch {0} done.\n".format(stru_file))
#         return totalAtom, volume/totalAtom, localAtoms

def fetch_stru_abacus(stru_file):
    try:
        get_volume = subprocess.Popen(args="grep Volume %s |tail -n1| tr -s ' '|cut -d ' ' -f5" % stru_file,
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        total_volume =float(get_volume.stdout.read())
    except:
        print('WARNING: check {0}, total volume is not found!'.format(stru_file))
        total_volume = 100
    if not os.path.exists(stru_file):
        print('WARNING: {0} doesnot exist.'.format(stru_file))
        return 100, 100, dict()
    with open(stru_file) as f:
        ntype = 0
        found_type = 0
        localAtoms = dict()
        temp = f.readline()
        while True:
            if "ntype" in temp:
                ntype = int(temp.split()[-1])
            if "atom label =" in temp:
                element = temp.split()[-1]
                f.readline()
                f.readline()
                f.readline()
                localAtoms[element] = int(f.readline().split()[-1])
                found_type += 1
                if found_type == ntype:
                    break
            temp = f.readline()
            if len(temp) == 0: break
    # print(ntype)
    totalAtom = sum(localAtoms.values())
    if totalAtom == 0: print(localAtoms)
    return totalAtom, total_volume/totalAtom, localAtoms

def fetch_out_abacus(out_file):
    try:
        get_total_e = subprocess.Popen(args="grep '!FINAL_ETOT_IS' %s|tr -s ' '|cut -d ' ' -f3" % out_file,
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        total_e = float(get_total_e.stdout.read())
        print("fetch {0} done.\n".format(out_file))
        return total_e
    except:
        print('WARNING: check {0}, total energy is not found!'.format(out_file))
        return 1e5

def cal_form(atoms, totnatom, energyPerAtom, reference):
    id = atoms.keys()
    formEnergy = energyPerAtom
    for each in id:
        formEnergy -= atoms[each]/totnatom * reference[each]
    return formEnergy * 1000

def main():
    # software = int(input("Please choose the software (0:ABACUS, 1:PROFESS):"))
    software = 0

    reference_file = {"./4_single/fcc_Al":"Al", "./4_single/hcp_Mg":"Mg", "./4_single/bcc_Li":"Li"}
    reference = dict()

    # structures = os.listdir('.')
    structures = []
    if software == 1:
        for each in os.walk('.'):
            found_inpt = False
            found_ion = False
            for each_file in each[2]:
                if "inpt" in each_file:
                    found_inpt = True
                if "ion" in each_file:
                    found_ion = True
            if found_inpt and found_ion:
                structures.append(each[0])

    elif software == 0:
        for each in os.walk('.'):
            if 'INPUT' in each[2] and 'STRU' in each[2] and 'KPT' in each[2]:
                structures.append(each[0])

    structures = sorted(structures)
    
    Nstru = len(structures)
    atoms = [0] * Nstru
    totnatom = np.zeros(Nstru)
    volumePerAtom = np.zeros(Nstru) # A^3
    energyPerAtom = np.zeros(Nstru) # eV
    
    # if software == 0: reference["Al"] = -57.949012

    formationEnergy = np.zeros(Nstru) # meV

    if software == 1:
        for i in range(Nstru):
            stru_file = '{0}/{1}.ion'.format(structures[i], structures[i].split('/')[-1])
            totnatom[i], volumePerAtom[i], atoms[i] = fetch_stru_profess(stru_file)
            out_file = '{0}/{1}.out'.format(structures[i], structures[i].split('/')[-1])
            energyPerAtom[i] = fetch_out_profess(out_file)/totnatom[i]
            if structures[i] in reference_file.keys():
                reference[reference_file[structures[i]]] = energyPerAtom[i]

    elif software == 0:
        for i in range(Nstru):
            stru_file = '{0}/OUT.oftest/running_scf.log'.format(structures[i])
            # stru_file = '{0}/OUT.oftest/STRU_ION_D'.format(structures[i])
            totnatom[i], volumePerAtom[i], atoms[i] = fetch_stru_abacus(stru_file)
            out_file = '{0}/OUT.oftest/running_*.log'.format(structures[i])
            energyPerAtom[i] = fetch_out_abacus(out_file)/totnatom[i]
            if structures[i] in reference_file.keys():
                reference[reference_file[structures[i]]] = energyPerAtom[i]

    for i in range(Nstru):
        formationEnergy[i] = cal_form(atoms[i], totnatom[i], energyPerAtom[i], reference)
    
    with open("data", 'w') as f:
        f.write('structure\tnatom\tcomponent\tvolume(A^3/atom)\tenergy(eV/atom)\tformation(meV)\n')
        for i in range(Nstru):
            component = ''
            for key in atoms[i].keys():
                component = component + '{0}{1}'.format(key, atoms[i][key])
            f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(
                structures[i],
                totnatom[i],
                component,
                volumePerAtom[i],
                energyPerAtom[i],
                formationEnergy[i]
            ))

main()