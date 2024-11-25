import numpy as np
import matplotlib.pyplot as plt

path = '/home/dell/1_work/TKK_abacus/1_EFTKK_Al/5_results/'
files = ['1PROFESS_KERNEL.dat']
# files = ['1PROFESS_KERNEL.dat','3PROFESS_KERNEL.dat','4PROFESS_KERNEL.dat','5PROFESS_KERNEL.dat','8PROFESS_KERNEL.dat','11PROFESS_KERNEL.dat']
rcut = [16,20]
colors = ['dimgray', 'royalblue', 'r', 'sandybrown']
# rcut = [8, 12,16,20, 24, 32]

# target_file = './4origin_PROFESS_KERNEL.dat'
target_file = './TWT'
twt = np.loadtxt(target_file, skiprows=4).T

plt.figure()
ax = plt.gca() 
ax.tick_params(top=False,right = False, bottom=False, left=False)
plt.axis('off')
# plt.plot(twt[0], twt[1], '-',c='orange', label='twt', linewidth=5)
for i in range(len(files)):
    kernel = np.loadtxt(path+files[i], skiprows=4).T
    plt.plot(kernel[0], kernel[1], '-',c='r', linewidth=5) 
# plt.xlabel("$\eta=q / (2k_F)$",fontsize=18)
# plt.ylabel("$W(\eta)$",fontsize=18)
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
plt.xlim(0,6)
plt.ylim(-2.8, 0.6)
# plt.text(-1, 0.6, "(c)", fontsize=14)
plt.tight_layout()
# plt.legend(frameon=False,ncol=1,fontsize=12)
plt.savefig('tkk.png', dpi=350)
# plt.savefig('twt.png', dpi=350)