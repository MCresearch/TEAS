import numpy as np
import matplotlib.pyplot as plt

path = '/home/dell/1_work/TKK_abacus/1_EFTKK_Al/5_results/'
files = ['1PROFESS_KERNEL.dat','3PROFESS_KERNEL.dat','4PROFESS_KERNEL.dat','5PROFESS_KERNEL.dat']
# files = ['1PROFESS_KERNEL.dat','3PROFESS_KERNEL.dat','4PROFESS_KERNEL.dat','5PROFESS_KERNEL.dat','8PROFESS_KERNEL.dat','11PROFESS_KERNEL.dat']
rcut = [8, 12,16,20]
colors = ['dimgray', 'royalblue', 'r', 'sandybrown']
# rcut = [8, 12,16,20, 24, 32]

target_file = '/home/dell/2_software/tools/TKK/2013-09-23-Jlq_Kernel_1.6Fix/TKK/rcut32/DUMP_KG.dat'
wt = np.loadtxt(target_file).T

plt.figure()
ax = plt.gca() 
ax.tick_params(top=True,right = True)
plt.plot(wt[0], wt[2]-1.6, '-',c='k', label='WT')
for i in range(len(files)):
    kernel = np.loadtxt(path+files[i], skiprows=4).T
    plt.plot(kernel[0], kernel[1], '--', c=colors[i], label='$\lambda_c$=%.1f' % rcut[i]) 
plt.xlabel("$\eta=q / (2k_F)$",fontsize=18)
plt.ylabel("$W(\eta)$",fontsize=18)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(0,6)
plt.ylim(-2.8, 0.6)
plt.text(-1, 0.6, "(c)", fontsize=14)
plt.tight_layout()
plt.legend(frameon=False,ncol=1,fontsize=12)
plt.savefig('kernel.png', dpi=350)