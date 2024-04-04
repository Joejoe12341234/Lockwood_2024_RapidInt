from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math

##### Leadtime plot

fig=plt.figure(figsize=(10, 14))
import scipy.io
mat = scipy.io.loadmat('/LEADTIME2.mat')
SSP = mat['SSP']
NCEP = mat['NCEP']
X = mat['X']
ax=fig.add_subplot(312)
ax.set_xlim([0,200])
line2, =  plt.plot(np.transpose(X),NCEP,'b',linewidth=3.0, label='NCEP')
line1, =  plt.plot(np.transpose(X),np.nanmean(SSP,0),'r',linewidth=3.0, label='SSP2-4.5')
for v in np.arange(6):
    plt.plot(np.transpose(X),SSP[v,:],'r',linewidth=0.3)
ax.set_ylim([0,0.09])
plt.tick_params(axis='both', which='major', labelsize=11)
fig.supylabel('Probability', fontsize=14)
fig.supxlabel('Hours', fontsize=14)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.text(0.4, 0.92, 'Leadtime -- 35 knots/hours threshold', fontsize=16,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
plt.text(0.84, 0.5, 'SSP2-4.5 mean = 46 h', fontsize=12,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
plt.text(0.822, 0.58, 'NCEP mean =  53 h', fontsize=12,horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
plt.legend(handles=[line2, line1],fontsize=13)

##### Intensity plots
import scipy.io
mat = scipy.io.loadmat('.mat')
f = mat['FULL_ALL']
fncep = mat['FULL_NCEP']
f2 = mat['FULL_XX']

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math

fig, axs = plt.subplots(1,2, figsize=(15, 5))
fig.subplots_adjust(left=.06, bottom=.13, right=None, top=None, wspace=.09, hspace=.01)

# fig.suptitle('Vertically stacked subplots')
a = (np.nanmean(f,0))

axs[0].set_ylim([-7, 0])
a = (np.nanmean(fncep,0))
axs[0].plot(np.transpose(f2),a,'b',linewidth=4.0)
a = (f)
axs[0].set_xlim([-180, 180])
aa = np.arange(-201,200,10)
for v in np.arange(6):
    axs[0].plot(np.transpose(f2),f[v,:],'r',linewidth=.5, label='SSP2-4.5', alpha=0.3)

a = (np.nanmean(f,0))
axs[0].plot(np.transpose(f2),a,'r',linewidth=4.0)
axs[0].text(0.25, 0.9, '(a) Full North Atlantic', fontsize=16, horizontalalignment='center', verticalalignment='center', transform=axs[0].transAxes)
axs[0].set_ylim([-5, 0])


f = mat['COASTL_FULL']
fncep = mat['COASTL_NCEP']
f2 = mat['FULL_XX']

a = (np.nanmean(f,0))
line1, = axs[1].plot(np.transpose(f2),a,'r',linewidth=3.0, label='SSP2-4.5')
axs[1].set_ylim([-7, 0])
a = fncep
line2, =  axs[1].plot(np.transpose(f2),a,'b',linewidth=3.0, label='NCEP')
a = (f)
axs[1].set_xlim([-180, 180])
aa = np.arange(-201,200,10)

for v in np.arange(6):
    axs[1].plot(np.transpose(f2),f[v,:],'r',linewidth=.5, label='SSP2-4.5', alpha=0.1)


fig.supylabel('Log$_{10}$ (P)', fontsize=14)
fig.supxlabel('Intensity change (knots/24 hours)', fontsize=14)
axs[1].tick_params(axis='both', which='major', labelsize=11)
axs[1].text(0.33, 0.9, '(b) Prior to landfall (48 hours)', fontsize=16,horizontalalignment='center', verticalalignment='center', transform=axs[1].transAxes)
axs[1].set_ylim([-5, 0])
axs[1].legend(handles=[line2, line1],fontsize=13)
plt.savefig('Intenslogplot_SSP2452.eps', format='eps')




