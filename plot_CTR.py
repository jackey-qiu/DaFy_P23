import numpy as np
import sys, os, locate_path
DaFy_path = locate_path.module_path_locator()
sys.path.append(os.path.join(DaFy_path,'util', 'XRD_tools'))
import reciprocal_space_v3 as rsp
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import sys,os
import configparser
from scipy.ndimage import gaussian_filter
from pylab import MultipleLocator, LogLocator, FormatStrFormatter
from util.PlotSetup import *

# which_scans_to_plot = [102,104,118,126,133]
# which_scans_to_plot = [109,111,120,129,135]
# which_scans_to_plot = [114,123,131,137]
which_scans_to_plot = [102,104,144]
# which_scans_to_plot = [109,111,147]
# which_scans_to_plot = [144,155,162]
# which_scans_to_plot = [147,156,164]
# which_scans_to_plot = [158,165]
# which_scans_to_plot = [126,162]
# which_scans_to_plot = [129,164]

# which_scans_to_plot = [97]
config_file_name = 'CV_XRD_plot_i20180678_Jul23_2019.ini'
config_file = os.path.join(DaFy_path, 'config', config_file_name)

#do you want to bin your datapoints
#debug is required to set bin_level>1
bin_level = 1

#specify current density limit, other limits are set automatically
ylim_current_density = [-15, 15]
#crystal reciprocal lattice instance

#################you seldom need to touch the following code lines###############
#extract info from config file
config = configparser.ConfigParser()
config.read(config_file)
global_vals = ['scan_ids', 'scan_number', 'rod_scan', 'data_files', 'hkl', 'colors']
plot_item_index = []
scans_to_plot_original = which_scans_to_plot
for each in global_vals:
    globals()[each] = []

for section in config.sections()[::-1]:
    if section == 'beamtime':
        beamtime =  eval(config.get(section,'beamtime'))
    else:
        scan_number_temp = eval(config.get(section,'scan_number'))
        if scan_number_temp[0] in which_scans_to_plot:
            plot_item_index.append(scans_to_plot_original.index(scan_number_temp[0]))
            del which_scans_to_plot[which_scans_to_plot.index(scan_number_temp[0])]
            for each_global_val in global_vals:
                if each_global_val == 'data_files':
                    globals()['data_files'].append(os.path.join(eval(config.get(section,'data_file_header')),\
                                                                eval(config.get(section,'data_files'))[0]))
                elif each_global_val == 'rod_scan':
                    globals()['rod_scan'].append(bool(eval(config.get(section,'rod_scan'))))
                elif each_global_val == 'hkl':
                    print(config.get(section,'hkl'))
                    globals()['hkl'].append(config.get(section,'hkl'))
                else:
                    globals()[each_global_val].append(eval(config.get(section,each_global_val))[0])
for val in global_vals:
    globals()[val]=[globals()[val][i] for i in plot_item_index]

labels = ["-0.6 V","-0.6 V","-0.8 V","-1 V","-1.7 V","-1.6V"]
colors = ["r","b","g","blue","m",'k']

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
ax.set_yscale("log")
# ax_2=ax.twinx()
ax.set_ylabel('I(au)')
# ax_2.set_ylabel('potential(V)')
start = 0
for i in range(len(data_files)):
    data = np.load(data_files[i])
    H = data["H"][start:]
    K = data["K"][start:]
    L = data["L"][start:]
    potential = data['potential'][start:]
    I = data["peak_intensity"][start:]
    Ierr = data["peak_intensity_error"][start:]
    if rod_scan[i]:
        ax.plot(L,I,".-",color = colors[i], label = scan_ids[i]+'{:3.1f} V'.format(data['potential'][0]))
        plt.legend()
    else:
        ax.plot(potential,I,".-",color = colors[i], label = scan_ids[i]+'{:3.1f} V'.format(data['potential'][0]))
        # ax_2.plot(potential,color ='blue')
hkl_title = '{:3.1f}_{:3.1f}_L'.format(H[0],K[0])
# ax.set_ylim(I.min(),10)
# ax.set_ylim(0.0001,100)
if rod_scan[0]:
    ax.set_xlabel('L')
else:
    ax.set_xlabel('Potential (V/AgCl')
if not rod_scan[0]:
    hkl_title = 'XRV@hkl=({:3.1f},{:3.1f},{:3.1f})'.format(H[0],K[0],L[0])
    ax.set_ylim(I.min(),I.max())
plt.title(hkl_title)
# if rod_scan[0]:
    # plt.legend()
# fig.savefig('./plots/ctr_i20180678/20L.png',dpi=300, bbox_inches='tight')
fig.savefig('20L.png',dpi=300, bbox_inches='tight')
plt.show()
