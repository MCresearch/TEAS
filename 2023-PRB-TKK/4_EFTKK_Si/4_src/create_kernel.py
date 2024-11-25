import numpy as np
from ctypes import *


class CreateKernel:
    def __init__(self, settings):
        self.jle = self.get_real_bessel(settings)
        self.jG = self.get_recip_bessel(settings)
        self.jle_inte = np.array(
            [self.cotes(settings.r, each * settings.r * settings.r) for each in self.jle])
        self.jG_inte = np.array([self.cotes(settings.g[1:], each[1:] * settings.g[1:]
                                            * settings.g[1:] * settings.g[1:] * settings.g[1:]) for each in self.jG])

    def get_real_bessel(self, settings):
        bessel = cdll.LoadLibrary("./create_Bessel.so")
        nmax = 1610  # defined in create_Bessel.cpp
        assert(settings.nr <= nmax)
        jle = ((c_double*nmax)*settings.n_bessel)()
        bessel.allocate_Bessel.argtypes = [
            c_double, c_double, c_int, c_int, (c_double*nmax)*settings.n_bessel]
        bessel.allocate_Bessel(settings.dr, settings.rcut,
                               settings.n_bessel, settings.l, jle)
        return np.array(jle)[:, :settings.nr]

    def get_real_kernel(self, settings):
        kernel = np.zeros_like(self.jle[0])
        for i in range(settings.n_bessel):
            kernel += self.jle[i] * settings.coef[i]
        return kernel

    def get_recip_bessel(self, settings):
        # g = np.linspace(0, ng * dg, ng + 1)
        # r = np.linspace(dr, nr * dr, nr)
        jG = np.zeros((settings.n_bessel, len(settings.g)))
        for i in range(settings.n_bessel):
            jG[i][1:] = self.FFT_R2G(
                settings.r, settings.dr, settings.g[1:], self.jle[i])
            jG[i, 0] = self.FFT_R2G0(settings.r, settings.dr, self.jle[i])
        return jG

    def get_recip_kernel(self, coef, settings):
        self.gkernel = np.zeros_like(self.jG[0])
        for i in range(settings.n_bessel):
            self.gkernel += self.jG[i] * coef[i]
        self.gkernel -= 1.6

    def output_kernel(self, settings):
        n = len(settings.g)
        with open(settings.kernelfile, 'w') as f:
            f.write('{0}\n'.format(settings.kinetic))
            f.write('1.182014212098715E+00\n-1.000000000000000E+00\n')
            f.write('{0}\n'.format(n))
            for i in range(n):
                f.writelines('%.5e\t%.12e\n' %
                             (settings.g[i], self.gkernel[i]))

    def FFT_R2G(self, r, dr, g, f_R):
        return np.array([sum(f_R * np.sin(r * _g) * r) / _g * 4 * np.pi * dr for _g in g])

    def FFT_R2G0(self, r, dr, f_R):
        return sum(f_R * r * r) * 4 * np.pi * dr

    def cotes(self, r, f):
        assert len(r) % 4 == 1
        n = len(r) // 4
        P = 0
        for i in range(n):
            P += (r[4*i+4]-r[4*i])*(7*f[4*i]+32*f[4*i+1] +
                                    12*f[4*i+2]+32*f[4*i+3]+7*f[4*i+4])/90
        return P
