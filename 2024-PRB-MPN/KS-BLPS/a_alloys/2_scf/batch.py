import numpy as np
import os

maxWork = 10

def collect_stru(dir):
    strus = []
    for each in os.scandir(dir):
        if each.is_dir() and each.name != 'log':
            strus.append(each.name)
    return strus

def collect_log(log_dir):
    log = os.listdir(log_dir)
    return log

def collect_works(work_dir):
    works = []
    for config in os.scandir(work_dir):
        if config.is_dir():
            works.append(config.name)
    return works

def sbatch(batch_dir):
    work_dir = os.getcwd()
    os.chdir(batch_dir)
    os.system('sbatch ./submit.sh')
    os.chdir(work_dir)
    print("sbatch {0} done".format(batch_dir))

def main():
    work_dir = os.getcwd()
    works = collect_works(work_dir)
    totWork = len(works)
    strus = [0] * totWork
    logs = [0] * totWork
    for i in range(totWork):
        strus[i] = collect_stru('./{0}/'.format(works[i]))
        logs[i] = collect_log('./{0}/log/'.format(works[i]))
    print('Collect done:')
    for i in range(totWork):
        print('{0}:\t{1}\t{2}/{3}'.format(i, works[i], len(logs[i]), len(strus[i])))
    
    n = 0
    while n < maxWork:
        indx = int(input("Select the order of the sturcture you want to sbatch ({0} left):".format(maxWork - n)))
        while len(strus[indx]) == len(logs[indx]):
            indx = int(input('All works of {0} are done, select another one please:'.format(works[indx])))

        work = works[indx]
        for each in strus[indx]:
            if each not in logs[indx]:
                sbatch("{0}/{1}/".format(work, each))
                n+=1
                if n >= maxWork:
                    break

main()
