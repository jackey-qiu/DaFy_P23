import numpy as np
import sys, os, locate_path
DaFy_path = locate_path.module_path_locator()
sys.path.append(os.path.join(DaFy_path,'util', 'XRD_tools'))
import reciprocal_space_v3 as rsp
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import sys,os
import configparser
from scipy.ndimage import gaussian_filter
from pylab import MultipleLocator, LogLocator, FormatStrFormatter
from util.PlotSetup import *

#good data scan numbers [216,224,231,243]
which_scans_to_plot = [229,221,231,243]
# which_scans_to_plot = [229,236,244]
config_file_name = 'CV_XRD_plot_i20180835_Jul18_2019.ini'
config_file_name = 'CV_XRD_plot_i20180835_Jul26_testMPI_2019.ini'
config_file = os.path.join(DaFy_path, 'config', config_file_name)

#do you want to set the max to 0
ref_max_eq_0 = {'strain':0,'size':0,'intensity':0}

#specify this for pot step scan
scan_time = 100 #in seconds

#do you want to bin your datapoints
#debug is required to set bin_level>1
bin_level = 1

#specify current density limit, other limits are set automatically
ylim_current_density = [-15, 15]
#crystal reciprocal lattice instance

crystals = ['Co3O4','CoHO2_9009884']
HKL_normals =[[1,1,1],[0,0,1]]
HKL_para_xs= [[1,1,-2],[1,0,0]]
offset_angles = [0,0]
for each in crystals:
    i = crystals.index(each)
    globals()[each] = rsp.lattice.from_cif(os.path.join(DaFy_path, 'util','cif',"{}.cif".format(each)), HKL_normal=HKL_normals[i],HKL_para_x=HKL_para_xs[i], offset_angle =offset_angles[i])


#################you seldom need to touch the following code lines###############
#extract info from config file
config = configparser.ConfigParser()
config.read(config_file)
global_vals = ['phs', 'scan_ids', 'scan_labels', 'ids_files', 'data_files', 'hkls', 'scan_direction_ranges','colors', 'xtal_lattices', 'plot_pot_steps']
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
                    if each_global_val == 'ids_files':
                        globals()['ids_files'].append(os.path.join(eval(config.get(section,'ids_file_header')),\
                                                                   eval(config.get(section,'ids_files'))[which_one]))
                    elif each_global_val == 'data_files':
                        globals()['data_files'].append(os.path.join(eval(config.get(section,'data_file_header')),\
                                                                    eval(config.get(section,'data_files'))[which_one]))
                    elif each_global_val == 'xtal_lattices':
                        globals()['xtal_lattices'].append(globals()[eval(config.get(section,'xtal_lattices'))[which_one]])
                    else:
                        globals()[each_global_val].append(eval(config.get(section,each_global_val))[which_one])
scan_ids_reordered = ["DaFy_{}".format(each) for each in which_scans_to_plot]
index_reordered =[scan_ids.index(each) for each in scan_ids_reordered]
for each_global in global_vals:
    globals()[each_global]=[globals()[each_global][i] for i in index_reordered]

#matplotlib global setting
matplotlib.rc('xtick', labelsize=10)
matplotlib.rc('ytick', labelsize=10)
plt.rcParams.update({'axes.labelsize': 10})
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.minor.width'] = 1
plt.rcParams['mathtext.default']='regular'
plt.style.use('seaborn-darkgrid')
print(os.path.dirname(os.path.abspath(__file__)))
overplot = 1
plot_vs_RHE = 1
create_ASCII = 0
plot_pot_step = plot_pot_steps[0]
if plot_pot_step:
    xlim = [0, scan_time+10]

scan_info = ScanInfoContainer()
for i in range(len(phs)):
    scan_info.add(scan_ids[i], np.load(data_files[i]), xtal_lattices[i], hkls[i], scan_direction_ranges[i], colors[i], \
                  scan_label=scan_labels[i],pH = phs[i],ids_filename =ids_files[i])
#####################################################################################################################################

x_label = ['E / V$_{RHE}$','Time(s)'][bool(plot_pot_step)]
y_lable_key = ['ip_strain','oop_strain','ip_sigma','oop_sigma','intensity','ip_size','oop_size','opt_ref','current_den']
ylabels = [r'$\Delta\varepsilon_\parallel$  (%)', r'$\Delta\varepsilon_\perp$  (%)',r'FWHM$_\parallel$ / $\AA^{-1}$',\
           r'FWHM$_\perp$ / $\AA^{-1}$',r'Intensity / a.u.',r'$\Delta d_\parallel$ / nm',r'$\Delta d_\perp$ / nm',r'Optical Reflectivity / %',r'j / mAcm$^{-2}$']
y_labels_lib = dict(zip(y_lable_key,ylabels))

fig_ax_container = {'ip_strain':['fig1','ax1a'],\
                    'oop_strain':['fig2','ax2a'],\
                    'ip_sigma': ['fig3','ax3a'],\
                    'oop_sigma':['fig4','ax4a'],\
                    'intensity':['fig5','ax5a'],\
                    'ip_size':['fig6','ax6a'],\
                    'oop_size':['fig7','ax7a'],\
                    'opt_ref':['fig8','ax8a'],\
                    'current_den':['fig9','ax9a'],\
                    'strain_all':['fig_ALL_strain','ax_fig_ALL_strain_ip','ax_fig_ALL_strain_oop'],\
                    'size_all':['fig_ALL_size','ax_fig_ALL_size_ip','ax_fig_ALL_size_oop'],\
                    'all_in_one':['fig_all','ax1_all','ax2_all','ax3_all','ax4_all','ax5_all','ax6_all']}

#use this the setup globally the limits
ax_can_ip_strain, ip_strain_min, ip_strain_max = [],1e10, -1e10
ax_can_oop_strain, oop_strain_min, oop_strain_max = [],1e10, -1e10
ax_can_ip_size, ip_size_min, ip_size_max = [],1e10, -1e10
ax_can_oop_size, oop_size_min, oop_size_max = [],1e10, -1e10
ax_can_current, current_min, current_max = [],1e10, -1e10
ax_can_intensity, intensity_min, intensity_max = [],1e10, -1e10
ax_can_set_yticks_ip_strain,ax_can_set_yticks_current, ax_can_set_yticks_oop_strain,  ax_can_set_yticks_ip_size,  ax_can_set_yticks_oop_size, ax_can_set_yticks_intensity  = [], [], [], [], [], []

#build figures
for each in fig_ax_container:
    if each in ['strain_all','size_all']:
        fig_ax_container[each][0] = plt.figure(figsize=(3*len(scan_ids),5))
        fig_ax_container[each][0].subplots_adjust(wspace=0.02,hspace=0.02)
    elif each == 'all_in_one':
        if len(which_scans_to_plot)>1:
            fig_ax_container[each][0] = plt.figure(figsize=(3*len(scan_ids),12))
            fig_ax_container[each][0].subplots_adjust(wspace=0.02,hspace=0.02)
        else:
            fig_ax_container[each][0] = plt.figure(figsize=(9,5))
    else:
        fig_ax_container[each][0] = plt.figure(figsize=(6,6))
    #fig_ax_container[each][0].tight_layout()

#loop through the datasets
num_datasets = len(scan_ids)

for scan_id in scan_ids:
    scan_no = int(scan_id.split('_')[-1])
    #get high resolution cv data from ids file
    #cv_data = extract_ids_file(scan_info[scan_id].ids_filename)
    #data file saved from DaFy program
    data = scan_info[scan_id].data.f
    pot, current_density = data.potential, data.current_density
    pot_cal = data.potential_cal
    frame_number = data.frame_number
    cv_data = pot, current_density
    Time = data.Time
    if len(Time)==0:
        Time = np.array(range(len(pot)))*float(scan_time/len(pot))
    # print(Time)
    cen_ip, cen_oop = data.cen_ip, data.cen_oop
    strain_ip, sigma_strain_ip = strain_ip_with_uncertainty(cen_ip, \
                                                            scan_info[scan_id].HKL_position,\
                                                            scan_info[scan_id].structure_lattice, 0)

    strain_oop, sigma_strain_oop = strain_oop_with_uncertainty(cen_oop, \
                                                            scan_info[scan_id].HKL_position,\
                                                            scan_info[scan_id].structure_lattice, 0)
    #if scan_id == "DaFy_216":
    #    strain_oop=list(np.array(strain_oop)+0.91)
    FWHM_ip, FWHM_oop = np.array(data.FWHM_ip), np.array(data.FWHM_oop)
    amp_oop, amp_ip = data.amp_oop, data.amp_ip
    pcov_ip, pcov_oop = data.pcov_ip, data.pcov_oop
    try:
        intensity = data.peak_intensity
        if sum(intensity)==0:
            intensity = amp_ip* FWHM_oop
    except:
        intensity = amp_ip*amp_oop
    # intensity = data.peak_intensity
    I = amp_ip * FWHM_oop * FWHM_ip
    # if plot_pot_step:
        # scan_direction_ranges = []
    return_cycle = 0
    if scan_id =='DaFy_231':
        return_cycle =1
    scan_direction_ranges, _, _ = select_cycle_new((pot_cal,pot_cal,frame_number),return_cycle=return_cycle)
    # print(scan_direction_ranges)
    #first point couple from the positive sweep to make the plot a complete circle
    current_axs = []
    point_gap_ip_strain = []
    point_gap_oop_strain = []
    point_gap_ip_size = []
    point_gap_oop_size = []

    # print(scan_direction_ranges)

    for pos in range(len(scan_direction_ranges)):
        if pos == len(scan_direction_ranges) -1:
            break
        elif pos == 1 and plot_pot_step:
            break
        else:
            indx1, indx2 = scan_direction_ranges[pos:pos+2]
            # indx2 += 1
        # print(indx1,indx2)
        fillcolor = 'w' if pos%2 ==1 else scan_info[scan_id].color
        label = scan_info[scan_id].scan_label+[' negative scan',' positive_scan'][int(pos==0)]
        marker = 'o'
        if plot_pot_step:
            marker = ''
        #strain data
        _, pot_temp, Time_temp =  select_cycle_new((pot_cal,Time,frame_number),return_cycle=return_cycle)
        _, pot_temp, strain_ip_temp =  select_cycle_new((pot_cal,strain_ip,frame_number),return_cycle=return_cycle)
        _, pot_temp, sigma_strain_ip_temp =  select_cycle_new((pot_cal,sigma_strain_ip,frame_number),return_cycle=return_cycle)
        _, pot_temp, strain_oop_temp =  select_cycle_new((pot_cal,strain_oop,frame_number),return_cycle=return_cycle)
        _, pot_temp, sigma_strain_oop_temp =  select_cycle_new((pot_cal,sigma_strain_oop,frame_number),return_cycle=return_cycle)

        #grain size
        _, pot_temp, FWHM_ip_temp =  select_cycle_new((pot_cal,FWHM_ip,frame_number),return_cycle=return_cycle)
        _, pot_temp, FWHM_oop_temp =  select_cycle_new((pot_cal,FWHM_oop,frame_number),return_cycle=return_cycle)
        size_ip_temp, size_oop_temp = 0.2*np.pi/np.array(FWHM_ip_temp), 0.2*np.pi/np.array(FWHM_oop_temp)
        #intensity
        _, pot_temp, intensity_temp =  select_cycle_new((pot_cal,intensity,frame_number),return_cycle=return_cycle)
        #x: either potential or time
        x = POT(pot_temp,plot_vs_RHE, scan_info[scan_id].pH) if not plot_pot_step else range(len(Time_temp))

        #init the axis handle
        if type(fig_ax_container['ip_strain'][1])==str:
            fig_ax_container['ip_strain'][1] = fig_ax_container['ip_strain'][0].add_subplot(111)
        if type(fig_ax_container['oop_strain'][1])==str:
            fig_ax_container['oop_strain'][1] = fig_ax_container['oop_strain'][0].add_subplot(111)
        if type(fig_ax_container['ip_size'][1])==str:
            fig_ax_container['ip_size'][1] = fig_ax_container['ip_size'][0].add_subplot(111)
        if type(fig_ax_container['oop_size'][1])==str:
            fig_ax_container['oop_size'][1] = fig_ax_container['oop_size'][0].add_subplot(111)
        if type(fig_ax_container['intensity'][1])==str:
            fig_ax_container['intensity'][1] = fig_ax_container['intensity'][0].add_subplot(111)
        if type(fig_ax_container['current_den'][1])==str:
            fig_ax_container['current_den'][1] = fig_ax_container['current_den'][0].add_subplot(111)
        if num_datasets == 1:# the only dataset occupy the whole figure
            for i in range(1,7):
                fig_ax_container['all_in_one'][i] = fig_ax_container['all_in_one'][0].add_subplot(2,3,i)
        else:#each dataset occupy one column
            for i in range(6):
                fig_ax_container['all_in_one'][i+1] = fig_ax_container['all_in_one'][0].add_subplot(6,num_datasets, scan_ids.index(scan_id)+1+i*num_datasets)
            for i in range(2):
                fig_ax_container['strain_all'][i+1] = fig_ax_container['strain_all'][0].add_subplot(2,num_datasets, scan_ids.index(scan_id)+1+i*num_datasets)
            for i in range(2):
                fig_ax_container['size_all'][i+1] = fig_ax_container['size_all'][0].add_subplot(2,num_datasets, scan_ids.index(scan_id)+1+i*num_datasets)

        ip_strain_axs = [fig_ax_container['ip_strain'][1]]
        oop_strain_axs = [fig_ax_container['oop_strain'][1]]
        ip_size_axs = [fig_ax_container['ip_size'][1]]
        oop_size_axs = [fig_ax_container['oop_size'][1]]
        intensity_axs = [fig_ax_container['intensity'][1]]
        current_axs.append(fig_ax_container['current_den'][1])

        if num_datasets ==1:
            ip_strain_axs.append(fig_ax_container['all_in_one'][1])
            oop_strain_axs.append(fig_ax_container['all_in_one'][2])
            ip_size_axs.append(fig_ax_container['all_in_one'][4])
            oop_size_axs.append(fig_ax_container['all_in_one'][5])
            intensity_axs.append(fig_ax_container['all_in_one'][6])
            current_axs.append(fig_ax_container['all_in_one'][3])
        else:
            ip_strain_axs.append(fig_ax_container['all_in_one'][2])
            oop_strain_axs.append(fig_ax_container['all_in_one'][3])
            ip_strain_axs.append(fig_ax_container['strain_all'][1])
            oop_strain_axs.append(fig_ax_container['strain_all'][2])
            ip_size_axs.append(fig_ax_container['all_in_one'][4])
            oop_size_axs.append(fig_ax_container['all_in_one'][5])
            ip_size_axs.append(fig_ax_container['size_all'][1])
            oop_size_axs.append(fig_ax_container['size_all'][2])
            current_axs.append(fig_ax_container['all_in_one'][1])
            intensity_axs.append(fig_ax_container['all_in_one'][6])
        #get the first point from first pos to feed in second pos
        if pos==0:
            point_gap_ip_strain = [x[0],strain_ip_temp[0]]
            point_gap_oop_strain = [x[0],strain_oop_temp[0]]
            point_gap_ip_size = [x[0],size_ip_temp[0]]
            point_gap_oop_size = [x[0],size_oop_temp[0]]
        else:
        #update the data in second pos
            strain_ip_temp = np.append(strain_ip_temp,point_gap_ip_strain[1])
            strain_oop_temp = np.append(strain_oop_temp,point_gap_oop_strain[1])

            size_ip_temp = np.append(size_ip_temp,point_gap_ip_size[1])
            size_oop_temp = np.append(size_oop_temp,point_gap_oop_size[1])
            x = np.append(x,point_gap_ip_strain[0])

        ip_strain_data_all = [strain_ip_temp]*len(ip_strain_axs)
        oop_strain_data_all = [strain_oop_temp]*len(oop_strain_axs)
        ip_size_data_all = [size_ip_temp]*len(ip_size_axs)
        oop_size_data_all = [size_oop_temp]*len(oop_size_axs)
        intensity_data_all = [intensity_temp]*len(intensity_axs)

        #check the limits for current dataset and update it if necessary
        #here reference point is 0(max value)
        def _update_min_or_not(current_min, current_max,data):
            new_min, new_max = 0, 0
            data = np.array(data)
            if current_min > (data.min() - data.max()):
                new_min = data.min() - data.max()
            else:
                new_min = current_min
            return new_min, new_max
        #here we don't set the reference point (max and min is the original values)
        def _update_min_max_or_not(current_min, current_max, data):
            new_min, new_max = 0, 0
            data = np.array(data)
            if current_min > data.min():
                new_min = data.min()
            else:
                new_min = current_min
            if current_max < data.max():
                new_max = data.max()
            else:
                new_max = current_max
            return new_min, new_max
        update_max_min = [_update_min_max_or_not,_update_min_or_not]
        ip_strain_min, ip_strain_max = update_max_min[ref_max_eq_0['strain']](ip_strain_min,ip_strain_max, strain_ip_temp)
        oop_strain_min, oop_strain_max = update_max_min[ref_max_eq_0['strain']](oop_strain_min,oop_strain_max, strain_oop_temp)
        ip_size_min, ip_size_max = update_max_min[ref_max_eq_0['size']](ip_size_min,ip_size_max, size_ip_temp)
        oop_size_min, oop_size_max = update_max_min[ref_max_eq_0['size']](oop_size_min,oop_size_max, size_oop_temp)
        intensity_min, intensity_max = update_max_min[ref_max_eq_0['intensity']](intensity_min,intensity_max, intensity_temp)
        #now collect all axs for different data (exclude the first ax which is the single plot)
        ax_can_ip_strain = ax_can_ip_strain + ip_strain_axs[1:]
        ax_can_oop_strain = ax_can_oop_strain + oop_strain_axs[1:]
        ax_can_ip_size = ax_can_ip_size + ip_size_axs[1:]
        ax_can_oop_size = ax_can_oop_size + oop_size_axs[1:]
        ax_can_intensity = ax_can_intensity + intensity_axs[1:]

        y_labels_all =[y_labels_lib['ip_strain']]* len(ip_strain_axs)+\
                      [y_labels_lib['oop_strain']]* len(oop_strain_axs)+\
                      [y_labels_lib['ip_size']]* len(ip_size_axs)+\
                      [y_labels_lib['oop_size']]* len(oop_size_axs)+\
                      [y_labels_lib['intensity']]* len(intensity_axs)

        def plot_on_ax(ax,x,data,y_label,ref_max):
            # print(ax)
            #you need this to match the size of data other than strain and size
            if len(x)!=len(data):
                x = x[0:-1]
            ax.plot(x[indx1:indx2],set_max_to_0(data,[indx1,indx2],ref_max),linestyle = 'none', linewidth =1,\
                    color = scan_info[scan_id].color, markerfacecolor = fillcolor,\
                    markeredgecolor = scan_info[scan_id].color,marker = marker, markersize=4,label = label)
            ax.set_ylabel(y_label)
            return None
        for ax, data, y_label,ref in zip(ip_strain_axs+oop_strain_axs+ip_size_axs+oop_size_axs+intensity_axs,\
                            ip_strain_data_all+oop_strain_data_all+ip_size_data_all+oop_size_data_all+intensity_data_all,\
                            y_labels_all,[ref_max_eq_0['strain']]*(len(ip_strain_axs)*2)+\
                            [ref_max_eq_0['size']]*(len(ip_strain_axs)*2)+[ref_max_eq_0['intensity']]*(len(ip_strain_axs)*1)):
            plot_on_ax(ax,x,data,y_label,ref)
            ax.set_xlabel(x_label)
        #now plot current density
        if pos==0:
            for ax in current_axs:
                colors_lib = {0:'sienna',1:'red',2:'green',3:'blue',4:'m',5:'black'}
                #ax.plot(POT(cv_data[0], plot_vs_RHE, scan_info[scan_id].pH), cv_data[1]*1000*8*50, '-',\
                #        color=scan_info[scan_id].color,linewidth=2,label=scan_info[scan_id].scan_label)
                ax.plot(POT(cv_data[0], plot_vs_RHE, scan_info[scan_id].pH), cv_data[1]*(-8)*50, 'x',\
                        color=scan_info[scan_id].color,linewidth=2,label=scan_info[scan_id].scan_label)
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_labels_lib['current_den'])
                #ax.set_ylim(ylim_current_density)
        #now set some x-y lables to '' for commen data range, also set titles
        x_ticks, x_tick_labels = find_tick_ticklables(x, num_ticks =4, endpoint = True, dec_place =1)
        if num_datasets == 1:
            [fig_ax_container['all_in_one'][i].set_xlabel('') for i in [1,2,3]]
            [fig_ax_container['all_in_one'][i].set_xticks(x_ticks) for i in [4,5,6]]
            [fig_ax_container['all_in_one'][i].set_xticklabels(x_tick_labels) for i in [4,5,6]]
            [fig_ax_container['all_in_one'][i].set_xticklabels([]) for i in [1,2,3]]
            ax_can_set_yticks_current = [fig_ax_container['all_in_one'][3]]

        else:
            fig_ax_container['all_in_one'][1].set_title(scan_info[scan_id].scan_label,fontsize = 10)
            fig_ax_container['size_all'][1].set_title(scan_info[scan_id].scan_label,fontsize = 10)
            fig_ax_container['strain_all'][1].set_title(scan_info[scan_id].scan_label,fontsize=10)
            [fig_ax_container['all_in_one'][i].set_xlabel('') for i in [1,2,3,4,5]]
            [fig_ax_container['size_all'][i].set_xlabel('') for i in [1]]
            [fig_ax_container['strain_all'][i].set_xlabel('') for i in [1]]
            #no tick labels
            [fig_ax_container['all_in_one'][i].set_xticklabels([]) for i in [1,2,3,4,5]]
            [fig_ax_container['size_all'][i].set_xticklabels([]) for i in [1]]
            [fig_ax_container['strain_all'][i].set_xticklabels([]) for i in [1]]
            #set the same tick labels
            fig_ax_container['all_in_one'][6].set_xticks(x_ticks)
            fig_ax_container['all_in_one'][6].set_xticklabels(x_tick_labels)
            fig_ax_container['size_all'][2].set_xticks(x_ticks)
            fig_ax_container['size_all'][2].set_xticklabels(x_tick_labels)
            fig_ax_container['strain_all'][2].set_xticks(x_ticks)
            fig_ax_container['strain_all'][2].set_xticklabels(x_tick_labels)
            if scan_id != scan_ids[0]:
                [fig_ax_container['all_in_one'][i].set_ylabel('') for i in [1,2,3,4,5,6]]
                [fig_ax_container['size_all'][i].set_ylabel('') for i in [1,2]]
                [fig_ax_container['strain_all'][i].set_ylabel('') for i in [1,2]]
                [fig_ax_container['all_in_one'][i].set_yticklabels([]) for i in [1,2,3,4,5,6]]
                [fig_ax_container['size_all'][i].set_yticklabels([]) for i in [1,2]]
                [fig_ax_container['strain_all'][i].set_yticklabels([]) for i in [1,2]]
            else:
                ax_can_set_yticks_ip_strain = [fig_ax_container['all_in_one'][2], fig_ax_container['strain_all'][1]]
                ax_can_set_yticks_oop_strain = [fig_ax_container['all_in_one'][3], fig_ax_container['strain_all'][2]]
                ax_can_set_yticks_ip_size = [fig_ax_container['all_in_one'][4], fig_ax_container['size_all'][1]]
                ax_can_set_yticks_oop_size = [fig_ax_container['all_in_one'][5], fig_ax_container['size_all'][2]]
                ax_can_set_yticks_current = [fig_ax_container['all_in_one'][1]]
                ax_can_set_yticks_intensity = [fig_ax_container['all_in_one'][6]]

    #save ascii files
    header = '%s #%d CV with XRD\r\nPotential / V, Current Density / mA/cm^2, Strain ip / %%, Strain oop / %%, d ip / nm, d oop / nm, Intensity (area ip)'%(beamtime, scan_no)
    X = np.array([pot, current_density, strain_ip, strain_oop, (0.2*np.pi/FWHM_ip), (0.2*np.pi/FWHM_oop), amp_ip]).T
    filename = 'data/ascii/%s_%d_CV_XRD.dat'%(beamtime, scan_no)
    np.savetxt(filename, X, newline='\r\n', header=header)

#now let us set the limits
offset_scale = 0.1
[each.set_ylim((ip_strain_min-(ip_strain_max-ip_strain_min)*offset_scale,ip_strain_max+(ip_strain_max-ip_strain_min)*offset_scale)) for each in ax_can_ip_strain]
[each.set_ylim((oop_strain_min-(oop_strain_max-oop_strain_min)*offset_scale,oop_strain_max+(oop_strain_max-oop_strain_min)*offset_scale)) for each in ax_can_oop_strain]
[each.set_ylim((ip_size_min-(ip_size_max-ip_size_min)*offset_scale,ip_size_max+(ip_size_max-ip_size_min)*offset_scale)) for each in ax_can_ip_size]
[each.set_ylim((oop_size_min-(oop_size_max-oop_size_min)*offset_scale,oop_size_max+(oop_size_max-oop_size_min)*offset_scale)) for each in ax_can_oop_size]
[each.set_ylim((intensity_min-(intensity_max-intensity_min)*offset_scale,intensity_max+(intensity_max-intensity_min)*offset_scale)) for each in ax_can_intensity]
# [each.set_ylim((ip_strain_min-(ip_strain_max-ip_strain_min)*offset_scale,ip_strain_max+(ip_strain_max-ip_strain_min)*offset_scale)) for each in ax_can_ip_strain]
# [each.set_ylim((oop_strain_min+oop_strain_min*offset_scale,-oop_strain_min*offset_scale)) for each in ax_can_oop_strain]
# [each.set_ylim((ip_size_min+ip_size_min*offset_scale,-ip_size_min*offset_scale)) for each in ax_can_ip_size]
# [each.set_ylim((oop_size_min+oop_size_min*offset_scale,-oop_size_min*offset_scale)) for each in ax_can_oop_size]
# [each.set_ylim((intensity_min+intensity_min*offset_scale,-intensity_min*offset_scale)) for each in ax_can_intensity]
#for each in ax_can_set_yticks_current:
#    y_ticks, y_tick_labels = find_tick_ticklables(ylim_current_density, num_ticks =6, endpoint = False, dec_place =0)
#    each.set_yticks(y_ticks)
#    each.set_yticklabels(['']+y_tick_labels[1:])

for each in ax_can_set_yticks_ip_strain:
    # y_ticks, y_tick_labels = find_tick_ticklables([ip_strain_min+ip_strain_min*offset_scale,-ip_strain_min*offset_scale], num_ticks =5, endpoint = False, dec_place =3)
    y_ticks, y_tick_labels = find_tick_ticklables([ip_strain_min-(ip_strain_max-ip_strain_min)*offset_scale,ip_strain_max+(ip_strain_max-ip_strain_min)*offset_scale], num_ticks =5, endpoint = False, dec_place =3)
    each.set_yticks(y_ticks)
    each.set_yticklabels(y_tick_labels)
    # print(ip_strain_min,y_ticks,y_tick_labels)

for each in ax_can_set_yticks_oop_strain:
    y_ticks, y_tick_labels = find_tick_ticklables([oop_strain_min-(oop_strain_max-oop_strain_min)*offset_scale,oop_strain_max+(oop_strain_max-oop_strain_min)*offset_scale], num_ticks =5, endpoint = False, dec_place =3)
    each.set_yticks(y_ticks)
    each.set_yticklabels(y_tick_labels)

for each in ax_can_set_yticks_ip_size:
    y_ticks, y_tick_labels = find_tick_ticklables([ip_size_min-(ip_size_max-ip_size_min)*offset_scale,ip_size_max+(ip_size_max-ip_size_min)*offset_scale], num_ticks =5, endpoint = False, dec_place =2)
    each.set_yticks(y_ticks)
    each.set_yticklabels(y_tick_labels)

for each in ax_can_set_yticks_oop_size:
    y_ticks, y_tick_labels = find_tick_ticklables([oop_size_min-(oop_size_max-oop_size_min)*offset_scale,oop_size_max+(oop_size_max-oop_size_min)*offset_scale], num_ticks =5, endpoint = False, dec_place =2)
    each.set_yticks(y_ticks)
    each.set_yticklabels(y_tick_labels)

for each in ax_can_set_yticks_intensity:
    y_ticks, y_tick_labels = find_tick_ticklables([intensity_min-(intensity_max-intensity_min)*offset_scale,intensity_max+(intensity_max-intensity_min)*offset_scale], num_ticks =5, endpoint = False, dec_place =2)
    each.set_yticks(y_ticks)
    each.set_yticklabels(y_tick_labels)
#tight layout for figures
if len(scan_ids)==1:
    fig_ax_container['all_in_one'][0].tight_layout()


#save figures
if len(scan_ids) > 1:
    path = 'plots/many/'
else:
    path =  'plots/' + scan_id + '/'
if not os.path.exists(path):
    os.makedirs(path)
fig_ax_container['ip_strain'][0].savefig(path+'strain_ip.png', dpi=300, bbox_inches='tight')
fig_ax_container['oop_strain'][0].savefig(path+'strain_oop.png', dpi=300, bbox_inches='tight')
fig_ax_container['ip_sigma'][0].savefig(path+'FWHM_ip.png', dpi=300, bbox_inches='tight')
fig_ax_container['oop_sigma'][0].savefig(path+'FWHM_oop.png', dpi=300, bbox_inches='tight')
fig_ax_container['intensity'][0].savefig(path+'Peak_Intensity.png', dpi=300, bbox_inches='tight')
fig_ax_container['ip_size'][0].savefig(path+'Grainsize_ip.png', dpi=300, bbox_inches='tight')
fig_ax_container['oop_size'][0].savefig(path+'Grainsize_oop.png', dpi=300, bbox_inches='tight')
fig_ax_container['current_den'][0].savefig(path+'Pot_CurrentDensity.png', dpi=300, bbox_inches='tight')
fig_ax_container['all_in_one'][0].savefig(path+'Pot_all_together.png', dpi=300, bbox_inches='tight')
fig_ax_container['strain_all'][0].savefig(path+'Pot_all_strain_together.png', dpi=300, bbox_inches='tight')
fig_ax_container['size_all'][0].savefig(path+'Pot_all_size_together.png', dpi=300, bbox_inches='tight')
plt.legend()
plt.show()

