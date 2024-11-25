import numpy as np
import matplotlib.pyplot as plt

path = '/home/dell/1_work/TKK_abacus/1_EFTKK_Al/5_results/'
files = ['1opte_v.dat','3opte_v.dat','4opte_v.dat']
# files = ['1opte_v.dat','3opte_v.dat','4opte_v.dat','5opte_v.dat']
# files = ['1opte_v.dat','3opte_v.dat','4opte_v.dat','5opte_v.dat','8opte_v.dat','11opte_v.dat']
rcut = [8, 12,16,20]
# rcut = [8, 12,16,20, 24, 32]

target_file = '/home/dell/1_work/TKK_abacus/1_EFTKK_Al/3_get_target/Al_target'
target = np.loadtxt(target_file, skiprows=2).T
natom = target[-1]
target_ev = target[0] / natom
target_ev[20], target_ev[21] = target_ev[21], target_ev[20]
structures = ['fcc']*5 + ['hcp']*5 + ['bcc']*5 + ['sc']*5 + ['fcc(100)','fcc(110)','fcc(111)',r'fcc($1\times1\times1$)',r'fcc($2\times1\times1$)',r'fcc($2\times2\times1$)']

x = np.arange(len(target_ev))
markers = ['h', 'p', 's', '^']
# markers = ['1', '2', '3', '4']
colors = ['sandybrown', 'royalblue', 'r']
# colors = ['dimgray', 'royalblue', 'r', 'sandybrown']

plt.figure()
ax = plt.gca() 
ax.tick_params(top=True,right = True)
for i in range(len(files)):
    result = np.loadtxt(path+files[i]) / natom
    result[20], result[21] = result[21], result[20]
    # plt.scatter(x, result, s=30, c=colors[i], marker=markers[i], label='$\lambda_c$=%.1f' % rcut[i])
    plt.scatter(x, result, s=40, c='none', edgecolors=colors[i], marker=markers[i], label='$\lambda_c$=%.1f' % rcut[i])
plt.scatter(x, target_ev, c='k',marker='x',s=40, label='KS-BLPS')
# plt.xlabel('structures',fontsize=18)
plt.ylabel('Energy per atom/eV',fontsize=18)
# plt.xticks(x, structures, rotation='vertical',fontsize=12)
plt.xticks([2, 7, 12, 17, 21, 24], ['fcc', 'hcp', 'bcc', 'sc', 'surfaces', 'vacancy'], fontsize=14, rotation=45)
plt.yticks(fontsize=14)
ymin = -58.2
ymax = -56.4
plt.ylim(ymin, ymax)
plt.xlim(-0.5, len(target_ev)-0.5)
plt.vlines([4.5, 9.5, 14.5, 19.5, 22.5], ymin, ymax, linestyles=':', colors=['gray'] * 5)
plt.legend(frameon=False,fontsize=14)

plt.text(0, -56.5, r"$\rm{a_0}=3.97\ \rm{\AA}$", fontsize=10)
plt.text(5, -56.5, r"$\rm{a_0}=2.81\ \rm{\AA}$", fontsize=10)
plt.text(5, -56.6, r"$\rm{c_0}=4.61\ \rm{\AA}$", fontsize=10)
plt.text(10, -56.5, r"$\rm{a_0}=3.18\ \rm{\AA}$", fontsize=10)
plt.text(15, -56.5, r"$\rm{a_0}=2.66\ \rm{\AA}$", fontsize=10)
plt.text(-6, -56.4, "(a)", fontsize=16)
plt.tight_layout()
# plt.show()
plt.savefig('fittings_bigfont.png', dpi=350)