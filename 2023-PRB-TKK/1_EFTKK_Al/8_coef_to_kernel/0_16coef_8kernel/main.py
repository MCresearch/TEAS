import numpy as np

from create_kernel import CreateKernel
from profess import Profess
from settings import Settings
from Metropolis import Metropolis

settings = Settings()
create_kernel = CreateKernel(settings)
settings.coef[-1] = (1.6 - np.dot(create_kernel.jG[:-1, 0],
                                  settings.coef[:-1]))/create_kernel.jG[-1, 0]
profess = Profess(settings)
metropolis = Metropolis(settings, create_kernel, profess)
