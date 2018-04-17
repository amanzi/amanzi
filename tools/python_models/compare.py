import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import axes_grid1
import sys,os
sys.path.append(os.path.join(os.environ['HOME'],'research','python'))
import colors

import permafrost_model_explicit
import permafrost_model
import wc_T

# model explicit
pm = permafrost_model_explicit.PermafrostModel()
wc = wc_T.WC_T(pm)

# model implicit
pm2 = permafrost_model.PermafrostModel()
wc2 = wc_T.WC_T(pm2)

# water content for each case
p0s = np.array([97000,95000,90000])
WC0s = np.array([wc.wc(273.2,p0) for p0 in p0s])

def get_reasonable_Ts():
    Ts1 = np.arange(273.15-10, 273.15-5, 1.)
    Ts2 = np.arange(273.15-5, 273.15-1, .1)
    Ts3 = np.arange(273.15-1, 273.15-.01, .001)
    Ts3a = np.arange(273.15-.01, 273.15+.01, .0001)
    Ts3b = np.arange(273.15+.01, 273.15+1, .001)
    Ts4 = np.arange(273.15+1, 273.15+5, .1)
    Ts5 = np.arange(273.15+5, 273.15+10, 1.)
    Ts = np.concatenate((Ts1, Ts2, Ts3, Ts3a, Ts3b, Ts4, Ts5))
    return Ts
Ts = get_reasonable_Ts()
print len(Ts)

fig = plt.figure(figsize=(10,8))
axs = []

ax1 = []
ax1.append(fig.add_subplot(221))
div1 = axes_grid1.make_axes_locatable(ax1[0])
ax1.append(div1.append_axes("right",size=3.,pad=0,sharey=ax1[0]))
axs.append(ax1)

ax2 = []
ax2.append(fig.add_subplot(222))
div2 = axes_grid1.make_axes_locatable(ax2[0])
ax2.append(div2.append_axes("right",size=3.,pad=0,sharey=ax2[0]))
axs.append(ax2)

ax3 = []
ax3.append(fig.add_subplot(223))
div3 = axes_grid1.make_axes_locatable(ax3[0])
ax3.append(div3.append_axes("right",size=3.,pad=0,sharey=ax3[0]))
axs.append(ax3)

ax4 = []
ax4.append(fig.add_subplot(224))
div4 = axes_grid1.make_axes_locatable(ax4[0])
ax4.append(div4.append_axes("right",size=3.,pad=0,sharey=ax4[0]))
axs.append(ax4)

def plot(c, Ts, pTs, s, pTs2, s2):
    axs[0][0].plot(Ts, pTs[:,0], '-', color=c)
    axs[0][0].plot(Ts, pTs2[:,0], '-', color=c)
    axs[0][1].plot(Ts, pTs[:,0], '-', color=c)
    axs[0][1].plot(Ts, pTs2[:,0], '-', color=c)
    axs[0][0].plot([273.15,273.15],[-3.e7, .5e7],'k')
    axs[0][1].plot([273.15,273.15],[-3.e7, .5e7],'k')

    axs[1][0].plot(Ts, s[:,0], '-', color=c)
    axs[1][0].plot(Ts, s2[:,0], '-', color=c)
    axs[1][1].plot(Ts, s[:,0], '-', color=c)
    axs[1][1].plot(Ts, s2[:,0], '-', color=c)
    axs[1][0].plot([273.15,273.15],[0,1],'k')
    axs[1][1].plot([273.15,273.15],[0,1],'k')

    axs[2][0].plot(Ts, s[:,1], '-', color=c)
    axs[2][0].plot(Ts, s2[:,1], '-', color=c)
    axs[2][1].plot(Ts, s[:,1], '-', color=c)
    axs[2][1].plot(Ts, s2[:,1], '-', color=c)
    axs[2][0].plot([273.15,273.15],[0,1],'k')
    axs[2][1].plot([273.15,273.15],[0,1],'k')

    axs[3][0].plot(Ts, s[:,2], '-', color=c)
    axs[3][0].plot(Ts, s2[:,2], '-', color=c)
    axs[3][1].plot(Ts, s[:,2], '-', color=c)
    axs[3][1].plot(Ts, s2[:,2], '-', color=c)
    axs[3][0].plot([273.15,273.15],[0,1],'k')
    axs[3][1].plot([273.15,273.15],[0,1],'k')


pTs = np.array([[101325.,T] for T in Ts])
s = np.array([pm.saturations_Tp(T,p) for p,T in pTs])
s2 = np.array([pm2.saturations_Tp(T,p) for p,T in pTs])
plot('r', Ts, pTs, s, pTs, s2)

colors = ['goldenrod','g','b']
for i,WC0 in enumerate(WC0s):
    pTs = np.array([[wc.pressure(T,WC0),T] for T in Ts])
    s = np.array([pm.saturations_Tp(T,p) for p,T in pTs])

    pTs2 = np.array([[wc2.pressure(T,WC0),T] for T in Ts])
    s2 = np.array([pm2.saturations_Tp(T,p) for p,T in pTs2])
    plot(colors[i],Ts,pTs,s,pTs2,s2)


axs[0][0].set_ylabel("pressure")
axs[0][0].set_xlabel("temperature")
axs[0][0].set_xticks([265.,270, 275, 280])
axs[0][1].set_xlim(273.14,273.16)
axs[0][1].set_xticks([273.14,273.16])

axs[1][0].set_ylabel("gas saturation")
axs[1][0].set_xlabel("temperature")
axs[1][0].set_ylim(-.01,1.01)
axs[1][0].set_xticks([265.,270, 275, 280])
axs[1][1].set_ylim(-.01,1.01)
axs[1][1].set_xlim(273.14,273.16)
axs[1][1].set_xticks([273.14,273.16])

axs[2][0].set_ylabel("liquid saturation")
axs[2][0].set_xlabel("temperature")
axs[2][0].set_ylim(-.01,1.01)
axs[2][0].set_xticks([265.,270, 275, 280])
axs[2][1].set_ylim(-.01,1.01)
axs[2][1].set_xlim(273.14,273.16)
axs[2][1].set_xticks([273.14,273.16])

axs[3][0].set_ylabel("ice saturation")
axs[3][0].set_xlabel("temperature")
axs[3][0].set_ylim(-.01,1.01)
axs[3][0].set_xticks([265.,270, 275, 280])
axs[3][1].set_ylim(-.01,1.01)
axs[3][1].set_xlim(273.14,273.16)
axs[3][1].set_xticks([273.14,273.16])

plt.show()
