import numpy as np
import matplotlib.pyplot as plt

path = '/home/dell/1_work/TKK_abacus/4_EFTKK_Si/5_results/'
prefix = [1, 4, 7]
# prefix = [1, 4, 7, 10]
files = [str(_) + 'opte_v.dat' for _ in prefix]
rcut = [8, 12,16,20]

target_file = '/home/dell/1_work/TKK_abacus/4_EFTKK_Si/3_get_target/Si_target'
target = np.loadtxt(target_file, skiprows=2).T
natom = target[-1]
target_ev = target[0]
structures = ['CD']*11 + [r'$\beta$-tin']*11 + ['cd(100)', r'CD($1\times1\times1$)',r'CD($2\times1\times1$)']

x = np.arange(len(target_ev))
markers = ['h', 'p', 's', '^']
# markers = ['1', '2', '3', '4']
colors = ['sandybrown', 'royalblue', 'r']
# colors = ['dimgray', 'royalblue', 'r', 'sandybrown']

plt.figure()
ax = plt.gca() 
ax.tick_params(top=True,right = True)
for i in range(len(files)):
    result = np.loadtxt(path+files[i])
    # plt.scatter(x, result/natom, s=30, marker=markers[i], c=colors[i], label='$\lambda_c$=%.1f' % rcut[i])
    plt.scatter(x, result/natom, s=40, marker=markers[i], c='none', edgecolors=colors[i], label='$\lambda_c$=%.1f' % rcut[i])
plt.scatter(x, target_ev/natom, c='k',marker='x',s=40, label='KS-BLPS')
# plt.xlabel('structures',fontsize=18)
plt.ylabel('Energy per atom/eV',fontsize=18)
plt.xticks([5, 16, 22, 23.5], ['CD', r'$\beta$-tin', 'surface', 'vacancy'], fontsize=14, rotation=45)
plt.yticks(fontsize=14)
plt.legend(frameon=False,ncol=1,fontsize=14)
ymin = -111
ymax = -108.5
plt.ylim(ymin, ymax)
plt.xlim(-0.5, len(target_ev)-0.5)
plt.vlines([10.5, 21.5, 22.5], ymin, ymax, linestyles=':', colors=['gray'] * 3)
plt.text(3, -108.65, r"$\rm{a_0}=5.41\ \rm{\AA}$", fontsize=10)
plt.text(14, -108.65, r"$\rm{a_0}=4.75\ \rm{\AA}$", fontsize=10)
plt.text(14, -108.8, r"$\rm{c_0}=2.59\ \rm{\AA}$", fontsize=10)
plt.text(-6.7/26*23, -108.5, "(b)", fontsize=16)
plt.tight_layout()
# plt.show()
plt.savefig('fittings_bigfont.png', dpi=350)