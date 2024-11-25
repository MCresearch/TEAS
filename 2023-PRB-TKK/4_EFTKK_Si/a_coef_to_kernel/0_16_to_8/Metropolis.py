import numpy as np
import random
import os


class Metropolis:
    def __init__(self, settings, create_kernel, profess):
        self.n = len(settings.aC_list)
        # self.target, self.V0 = self.get_target(settings, profess)
        self.target = self.get_target(settings)
        self.lowest, temp, diff = self.get_residual(
            settings.coef, settings, create_kernel, profess)
        self.dump_files(settings.coef, e_v=temp, prefix='origin')
        logfile = open('log', 'w')
        logfile.write(
            'T    E_tot    E_har    J    G    residual    lowest\n')
        logfile.write('target    ')
        logfile.writelines('%.6e    ' % sum(each) for each in self.target)
        logfile.write('0    0\n')
        logfile.write('original    ')
        logfile.writelines('%.6f    ' % each for each in diff)
        logfile.write('%.6f    ' % self.lowest)
        logfile.write('%.6f\n' % self.lowest)
        logfile.close()
        os.rename(settings.kernelfile, 'origin_' + settings.kernelfile)
        self.true_lowest = self.lowest
        self.opt_coef = settings.coef.copy()
        self.metropolis(settings, create_kernel, profess)
        opt_e_v = self.get_residual(
            self.opt_coef, settings, create_kernel, profess)[1]
        self.dump_files(self.opt_coef, opt_e_v, prefix='opt')
        os.rename(settings.kernelfile, 'opt_' + settings.kernelfile)
        e_v = self.get_residual(settings.coef, settings,
                                create_kernel, profess)[1]
        self.dump_files(settings.coef, e_v)

    def get_target(self, settings):
        targetfile = np.loadtxt(settings.targetfile, skiprows=2).T
        # [total energy, kinetic energy, hartree energy, natoms]
        self.natoms = targetfile[2]
        target = targetfile[:2]
        #target = targetfile[:3] * settings.Ha2eV
        # target_temp = np.hstack([targetC, targetB])*settings.Ha2eV/2
        # target = [each for each in target_temp]
        # V0C = settings.vC_list[5]
        # V0C, target_BC = profess.fit(
        #     settings.vC_list, target[0][:self.n], settings.vC_list[0], settings.vC_list[-1], V0C)[1:]
        # V0B = settings.vB_list[5]
        # V0B, target_BB = profess.fit(
        #     settings.vB_list, target[0][self.n:], settings.vB_list[0], settings.vB_list[-1], V0B)[1:]
        # target.append(np.array([target_BC, target_BB]))
        return target

    def dump_files(self, coef, e_v=np.zeros(0), prefix=''):
        with open(prefix+'coef.dat', 'w') as f:
            f.writelines('%.12e\n' % each for each in coef)
        if e_v.size:
            with open(prefix+'e_v.dat', 'w') as f:
                f.writelines('%.6e\n' % each for each in e_v)

    def accept(self, residual, T):
        if residual <= self.lowest:
            return True
        else:
            x = random.random()
            probability = np.exp(-(residual-self.lowest)/T)
            return True if x < probability else False

    def move(self, k, settings, create_kernel, profess):
        coef_temp = settings.coef.copy()
        coef_temp[k] *= (1 + settings.rate[k] * random.uniform(-1, 1))
        coef_temp[k] += settings.bias * random.uniform(-1, 1)
        coef_temp[-1] = (1.6 - np.dot(create_kernel.jG[:-1, 0],
                                      coef_temp[:-1]))/create_kernel.jG[-1, 0]
        residual, temp, diff = self.get_residual(
            coef_temp, settings, create_kernel, profess)
        return coef_temp, residual, temp, diff

    def get_residual(self, coef, settings, create_kernel, profess):
        create_kernel.get_recip_kernel(coef, settings)
        create_kernel.output_kernel(settings)
        temp = profess.cal_e_v(self.n, settings.name, settings.structures)
        # BC = profess.fit(
        #     settings.vC_list, temp[0][:self.n], settings.vC_list[0], settings.vC_list[-1], self.V0[0])[-1]
        # BB = profess.fit(
        #     settings.vB_list, temp[0][self.n:], settings.vB_list[0], settings.vB_list[-1], self.V0[1])[-1]
        # temp.append(np.array([BC, BB]))
        diff = np.sum(abs(temp - self.target) / self.natoms, axis=1)
        diff = np.append(diff, abs(np.dot(coef, create_kernel.jle_inte)))
        diff = np.append(diff, abs(np.dot(coef, create_kernel.jG_inte)))
        residual = np.dot(diff, settings.weight)
        return residual, temp[0], diff

    def metropolis(self, settings, create_kernel, profess):
        logfile = open('log', 'a')
        for i in range(settings.nT):
            t = settings.T * settings.cooling_rate ** i
            nums = np.zeros(settings.n_bessel-1)
            for j in range(settings.nsteps):
                if j % settings.ncheck == 0 and j != 0:
                    for k in range(7):
                        if nums[k]/settings.ncheck >= settings.accept_high:
                            settings.rate[k] *= 1.2
                        elif nums[k]/settings.ncheck <= settings.accept_low:
                            settings.rate[k] *= 0.6
                        nums[k] = 0
                for k in range(7):
                    coef_temp, residual, temp, diff = self.move(
                        k, settings, create_kernel, profess)
                    if residual <= self.true_lowest:
                        self.true_lowest = residual
                        self.opt_coef = coef_temp.copy()
                    # if residual <= 1:
                    #     dump_files(settings.coef=coef_temp, e_v=temp,
                    #                prefix='T%.3fn%d' % (t, j))
                    #     os.rename(kernelfile, 'T%.3fn%d' % (t, j) + kernelfile)
                    if self.accept(residual, t):
                        settings.coef = coef_temp.copy()
                        nums[k] += 1
                        self.lowest = residual
                    logfile.write('%.6e    ' % t)
                    logfile.writelines('%.6e    ' % each for each in diff)
                    logfile.write('%.6e    ' % residual)
                    logfile.write('%.6e\n' % self.lowest)
                    logfile.flush()
        logfile.close()
