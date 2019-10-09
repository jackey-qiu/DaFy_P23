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

which_scans_to_plot = [270,280,285,296]
title = '40L'
config_file_name = 'CH5314_prep2_plot.ini'
config_file = os.path.join(DaFy_path, 'config', config_file_name)

#specify this for pot step scan
#do you want to bin your datapoints
#debug is required to set bin_level>1
bin_level = 1

#specify current density limit, other limits are set automatically

#################you seldom need to touch the following code lines###############
#extract info from config file
config = configparser.ConfigParser()
config.read(config_file)
global_vals = ['scan_ids', 'data_files','colors']
for each in global_vals:
    globals()[each] = []

for section in config.sections():
    if section == 'beamtime':
        beamtime =  eval(config.get(section,'beamtime'))
    else:
        scan_number_temp = eval(config.get(section,'scan_number'))
        for each_scan in which_scans_to_plot:
            if each_scan in scan_number_temp:
                which_one = scan_number_temp.index(each_scan)
                for each_global_val in global_vals:
                    if each_global_val == 'data_files':
                        globals()['data_files'].append(os.path.join(eval(config.get(section,'data_file_header')),\
                                                                    eval(config.get(section,'data_files'))[which_one]))
                    else:
                        globals()[each_global_val].append(eval(config.get(section,each_global_val))[which_one])



fig = plt.figure(2,figsize=(8,5))
ax = fig.add_subplot(111)
ax.set_yscale('log', nonposy = 'clip')
ax.set_title(title)

for each in data_files:
    which_one = data_files.index(each)
    data = np.load(data_files[which_one])
    L = data['L']
    I = data['peak_intensity']
    Ierr = data['peak_intensity_error']
    ax.errorbar(L, I, yerr = Ierr, xerr = None, fmt = '.-', label = scan_ids[which_one])
ax.set_xlabel('L (r.l.u)')
ax.set_ylabel('Intensity (a.u)')
plt.legend()
plt.show()

#save figures
path = 'plots/ctr/'
if not os.path.exists(path):
    os.makedirs(path)

fig.savefig(path+'{}.png'.format(title), dpi=300, bbox_inches='tight')

