import os
import numpy as np
from matplotlib import pyplot as plt

os.spawnl(os.P_WAIT, "../plot_wrm", "vangenuchten.xml")
sat = np.loadtxt('van_Genuchten_sat.txt')
pc = np.loadtxt('van_Genuchten_pc.txt')
krel = np.loadtxt('van_Genuchten_krel.txt')

plt.figure()

plt.subplot(121)
plt.plot(sat, pc, 'b-x')
plt.xlim(0,1)
plt.xlabel("saturation")
plt.ylabel("p_c")

plt.subplot(122)
plt.plot(sat, krel, 'b-x')
plt.xlim(0,1)
plt.xlabel("saturation")
plt.ylabel("k_rel")

plt.show()
