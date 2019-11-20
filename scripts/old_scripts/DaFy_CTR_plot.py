import numpy as np
import pandas as pd 
import glob
import matplotlib.pyplot as plt
import os

DaFy_path = 'C:\\apps\\DaFy_P23'
data_folders = ['C:\\apps\\DaFy_P23\\data\\00L','C:\\apps\\DaFy_P23\\data\\20L','C:\\apps\\DaFy_P23\\data\\11L','C:\\apps\\DaFy_P23\\data\\31L']
potentials_all = [[-0.6, -0.7,-1.7],[-0.6, -0.7,-1.7],[-0.6, -0.7,-1.7],[-0.6, -0.7,-1.7]]
hks_all = [[[0,0],[2,0],[1,1],[1,3]],[[0,0],[2,0],[1,1],[1,3]],[[0,0],[2,0],[1,1],[1,3]],[[0,0],[2,0],[1,1],[1,3]]]
colors = ['r','b','g']

fig = plt.figure(figsize=(5,10))

ax_00 = fig.add_subplot(411)
ax_20 = fig.add_subplot(412)
ax_11 = fig.add_subplot(413)
ax_13 = fig.add_subplot(414)
axs = [ax_00,ax_20,ax_11,ax_13]
#axs = [ax_20]
[ax.set_yscale('log',nonposy='clip') for ax in axs]
[ax.set_xlabel('L(r.s.u)') for ax in axs]
[ax.set_ylabel('Intensity') for ax in axs]
[ax.set_title('{}{}L'.format(*hks_all[0][axs.index(ax)])) for ax in axs]

for data_folder in data_folders:
    data_files = glob.glob(os.path.join(data_folder,'*.xlsx'))
    data = pd.read_excel(data_files[0])
    for file in data_files[1:]:
        data = pd.concat([data,pd.read_excel(file)])
    ii = data_folders.index(data_folder)
    hks = hks_all[ii]
    potentials = potentials_all[ii]
    for potential in potentials:
        #axs[ii].scatter(data[data['potential'].round(1) == potential]['L'],data[data['potential'].round(1) == potential]['peak_intensity'],marker ='.', ls='-',label = 'E = {}V'.format(potential))
        data_temp = data[data['potential'].round(1) == potential]
        #data_temp = data_temp.sort_values('L')
        def _cal_lf(h, k, l, a=3.6148, wl = 0.551):
            d = a/(h**2+k**2+l**2)**0.5
            return 1/(2*(wl/(2*d)*(1-wl**2/4/(d**2))**0.5))
        lf = _cal_lf(data_temp['H'],data_temp['K'],data_temp['L'])
        print(lf)
        axs[ii].plot(data_temp['L'],data_temp['peak_intensity'],'-o',markersize = 1, lw=1, color = colors[potentials.index(potential)],label = 'E = {}V'.format(potential))
        #data_temp.to_excel('E={}V_HK={}.xlsx'.format(potential,hks[data_folders.index(data_folder)]))
[ax.legend() for ax in axs]
plt.tight_layout()
plt.show()
fig.savefig(os.path.join(DaFy_path,'temp','temp_ctr.png'),dpi = 300)
    

