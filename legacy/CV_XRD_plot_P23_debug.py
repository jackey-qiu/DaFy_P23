import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import sys
from scipy.ndimage import gaussian_filter
from pylab import MultipleLocator, LogLocator, FormatStrFormatter
import matplotlib
import os
from P23config import *

matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)
plt.rcParams.update({'axes.labelsize': 15})
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

plot_vs_RHE = 1
create_ASCII = 1
plot_pot_step = 0
scan_time = 100 #in seconds
overplot = 0

class ScanInfo():
    def __init__(self, scan_id, data, structure_lattice, HKL_position, scan_direction_ranges, color, analog_data_filename=None, analog_data_plot_range=[0,-1], ids_filename=None, analog_data_interval=None, scan_label=None, pH=13):
        self.scan_id = scan_id
        self.data = data
        self.structure_lattice = structure_lattice
        self.HKL_position = HKL_position
        self.scan_direction_ranges = scan_direction_ranges
        self.color = color
        self.analog_data_filename = analog_data_filename
        self.analog_data_plot_range = analog_data_plot_range
        self.ids_filename = ids_filename
        self.analog_data_interval = analog_data_interval
        self.scan_label = scan_label
        self.pH = pH

class ScanInfoContainer(dict):
    def __init__(self):
        return
    def add(self, scan_id, *args, **kwargs):
        self[scan_id] = ScanInfo(scan_id, *args, **kwargs)

scan_info = ScanInfoContainer()
##############################fill out this block####################
#datasets based on DaFy
Co3O4_lat = P23config.Co3O4_lat
scan_info.add('DaFy_P23config_724_10', np.load('data/DataBank_724_20190214-112008.npz'), Co3O4_lat, [4, 0, 4], [0, 71, -2], 'red',scan_label='pH 13 0.1 M NaOH',pH = 13,ids_filename ='/home/qiu/data/beamtime/P23_11_18_I20180114/ids/I20180014/010.ids')
scan_info.add('DaFy_P23config_727_10', np.load('data/DataBank_727_20190214-122116.npz'), Co3O4_lat, [4, 0, 4], [0, 71, -2], 'red',scan_label='pH 13 Phosphate',pH = 13,ids_filename ='/home/qiu/data/beamtime/P23_11_18_I20180114/ids/I20180014/014.ids')
scan_info.add('DaFy_P23config_732_10', np.load('data/DataBank_732_20190214-122612.npz'), Co3O4_lat, [4, 0, 4], [0, 101, -2], 'red',scan_label='pH 8 Phosphate',pH = 8,ids_filename ='/home/qiu/data/beamtime/P23_11_18_I20180114/ids/I20180014/017.ids')
scan_ids = ['DaFy_P23config_724_10','DaFy_P23config_727_10','DaFy_P23config_732_10']
#####################################################################
def RHE(E_AgAgCl, pH=13):
    # electrode is 3.4 M KCl
    return 0.205 + E_AgAgCl + 0.059*pH

def POT(E_AgAgCl, plot_vs_RHE=False, pH=13):
    if(plot_vs_RHE):
        return RHE(E_AgAgCl, pH)
    else:
        return E_AgAgCl

def select_cycle(pot_result_couple, bin_level = 1, bin_mode = 'average', return_cycle = 1, pot_step = 0.005, plot_mode='CV'):
    #bin_mode = 'average' or 'select'
    #   'average':average every bin_level data points
    #   'select': select every bin_level data point
    #plot_mode = 'CV' or 'pot_step'
    #   'CV': plot CV data
    #   'pot_step': plot potential step data
    pot, result = pot_result_couple
    pot = np.around(pot,decimals=6)
    index_of_valley_pot = []
    for i in range(1,len(pot)-1):
        if (pot[i]<pot[i-1]) and (pot[i]<pot[i+1]):
            #filter out the noisy flutuation at those same preset potentials
            if abs(pot[i]*2-(pot[i-1]+pot[i+1]))>pot_step:
                index_of_valley_pot.append(i)
        else:
            pass
    if plot_mode == 'pot_step':
        index_of_valley_pot = [0, len(pot)]
        bin_level = 1# you dont want to bin for pot_step dataset
    print(index_of_valley_pot)
    if return_cycle > len(index_of_valley_pot)-2:
        return_cycle = len(index_of_valley_pot)-2
    pot_partial, result_partial = pot[index_of_valley_pot[return_cycle]:index_of_valley_pot[return_cycle+1]],\
                                  result[index_of_valley_pot[return_cycle]:index_of_valley_pot[return_cycle+1]]

    def _bin_array(data, bin_level, bin_mode):
        if bin_mode =='average':
            data_bin =np.array(data[0::bin_level])
            for i in range(1,bin_level):
                temp = np.array(data[i::bin_level])
                if len(temp)<len(data_bin):
                    temp = np.append(temp, np.zeros(len(data_bin) - len(temp))+temp[-1])
                    data_bin = data_bin + temp
            return data_bin/bin_level
        elif bin_mode =='select':
            return np.array(data[0::bin_level])

    if bin_level>1:
        pot_bin = _bin_array(pot_partial, bin_level,bin_mode)
        result_bin = _bin_array(result_partial, bin_level,bin_mode)
        return pot_bin, result_bin
    else:
        return pot_partial, result_partial



def set_max_to_0(data_list,slice_index):
    return np.array(data_list)[slice_index[0]:slice_index[1]]-np.array(data_list).max()

def extract_ids_file(file_path,which_cycle=3):
    data = []
    data_lines =[]
    current_cycle = 0
    with open(file_path) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('primary_data'):
                print(current_cycle)
                current_cycle=current_cycle+1
                if current_cycle == which_cycle:
                    for j in range(i+3,i+3+int(lines[i+2].rstrip())):
                        data.append([float(each) for each in lines[j].rstrip().rsplit()])
                    break
                else:
                    pass
            else:
                pass
    return np.array(data)[:,0], np.array(data)[:,1]

def strain_ip(q_ip, HKL, lattice):
    q_bulk = lattice.q(HKL)
    q_ip_bulk = np.sqrt(q_bulk[0]**2 + q_bulk[1]**2)
    return (q_ip_bulk/q_ip - 1.0)*100.

def strain_ip_with_uncertainty(q_ip, HKL, lattice, uncertainty_q_ip):
    q_bulk = lattice.q(HKL)
    q_ip_bulk = np.sqrt(q_bulk[0]**2 + q_bulk[1]**2)
    _strain_ip = (q_ip_bulk/q_ip - 1.0)*100.
    uncedrtainty_strain_ip = np.abs(q_ip_bulk/q_ip**2*100*uncertainty_q_ip)
    return (_strain_ip, uncedrtainty_strain_ip)

def strain_oop(q_oop, HKL, lattice):
    return (lattice.q(HKL)[2]/q_oop - 1.0)*100.

def strain_oop_with_uncertainty(q_oop, HKL, lattice, uncertainty_q_oop):
    q_oop_bulk = lattice.q(HKL)[2]
    _strain_oop = (q_oop_bulk/q_oop - 1.0)*100.
    uncertainty_strain_oop = np.abs(q_oop_bulk/q_oop**2*100*uncertainty_q_oop)
    return  (_strain_oop, uncertainty_strain_oop)
#####################################################################################################################################

if(1):
    xlim = [0.45, 2.15]
    if len(scan_ids)==1:
        overplot = 1
    else:
        overplot = 0

    if len(scan_ids)>1:
        fig_ALL_strain=plt.figure(figsize=(15,8))
        fig_ALL_size=plt.figure(figsize=(15,8))
        fig_ALL_strain.tight_layout()
        fig_ALL_size.tight_layout()

    set_ylims = 1
    plot_phase_labels = 1
    #fig1 for ip strain
    #fig2 for oop strain
    #fig3 for ip FWHM
    #fig4 for oop FWHM
    #fig5 for intensity
    #fig6 for inplane crystallite size
    #fig7 for out of plane crystallite size
    #fig8 for optical reflectivity
    #fig9 for current density
    figs_name = ['fig1','fig2','fig3','fig4','fig5','fig6','fig7','fig8','fig9']
    axs_name = ['ax1a', 'ax2a','ax3a','ax4a','ax5a','ax6a','ax7a','ax8a','ax9a']
    ylabels = [r'$\Delta\varepsilon_\parallel$  (%)', r'$\Delta\varepsilon_\perp$  (%)',r'FWHM$_\parallel$ / $\AA^{-1}$',\
               r'FWHM$_\perp$ / $\AA^{-1}$',r'Intensity / a.u.',r'$\Delta d_\parallel$ / nm',r'$\Delta d_\perp$ / nm',r'Optical Reflectivity / %']
    for i in range(len(figs_name)):
        globals()[figs_name[i]] = plt.figure()
        globals()[axs_name[i]] = globals()[figs_name[i]].add_subplot(111)
        if i<8:
            globals()[axs_name[i]].set_ylabel(ylabels[i], fontsize=20)
        if plot_pot_step:
            globals()[axs_name[i]].set_xlabel('Time(s)', fontsize=20)
        else:
            if plot_vs_RHE:
                globals()[axs_name[i]].set_xlabel('E / V$_{RHE}$', fontsize=20)
            else:
                globals()[axs_name[i]].set_xlabel('E vs (Ag/AgCl) / V', fontsize=20)

    #figure for current density
    if plot_pot_step:
        if plot_vs_RHE:
            ax9a.set_ylabel('E / V$_{RHE}$', fontsize=20)
        else:
            ax9a.set_ylabel('E vs (Ag/AgCl) / V', fontsize=20)
    else:
        ax9a.set_ylabel(r'Current Density / mAcm$^{-2}$', fontsize=20)

    #figure for all stuff
    if overplot:
        fig_all=plt.figure(figsize=(15,10))
    else:
        fig_all=plt.figure(figsize=(10,15))

    if(set_ylims and not plot_pot_step):
        ax9a.set_ylim([-12, 30])

    x_lines_pos = []#[-0.2, 0.48]
    x_lines_neg = []#[-0.56, 0.08]
    oxidation_potentials = np.array([-0.192, 0.222, 0.562])-0.097
    oxidation_potentials_labels = np.array([r'Co(OH)$_2$/'+'\n'+r'Co$_3$O$_4$', r'Co$_3$O$_4$/'+'\n'+r'CoOOH', r'CoOOH/'+'\n'+r'CoO$_2$'])

    if(xlim[0] > -0.3):
        oxidation_potentials = oxidation_potentials[1:]
        oxidation_potentials_labels = oxidation_potentials_labels[1:]

    axes = [ax1a, ax2a, ax3a, ax4a, ax5a, ax6a, ax7a, ax8a, ax9a]#, ax10a]
    label_pos = [0.95, 0.85, 0.85, 0.85, 0.05, 0.05, 0.55, 0.9, 0.9, 0.05]
    label_pos = [0.05, 0.85, 0.85, 0.85, 0.85, 0.05, 0.55, 0.9, 0.9, 0.05]
    for ax_no, ax in enumerate(axes):
        trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)
        if not plot_pot_step:
            ax.set_xlim(xlim)
        if(0):
            for x in x_lines_pos:
                ax.axvline(x, linestyle='--', color='k')
            for x in x_lines_neg:
                ax.axvline(x, linestyle='-', color='k')
            for i, x in enumerate(oxidation_potentials):
                ax.axvline(x, linestyle='-', color='k', linewidth=2)
                if(plot_phase_labels):
                    ax.text(x, label_pos[ax_no], oxidation_potentials_labels[i], transform=trans, ha='center', backgroundcolor='w')
        if not plot_pot_step:
            majorLocator_x1 = MultipleLocator(0.2)
            minorLocator_x1 = MultipleLocator(0.05)
            ax.xaxis.set_major_locator(majorLocator_x1)
            ax.xaxis.set_minor_locator(minorLocator_x1)

    for fig in [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9]:
        fig.tight_layout()

    for scan_id in scan_ids:
        if len(scan_ids)>1:
            label = scan_info[scan_id].scan_label
            ax_fig_ALL_strain_ip = fig_ALL_strain.add_subplot(2,len(scan_ids),scan_ids.index(scan_id)+1)
            ax_fig_ALL_strain_oop = fig_ALL_strain.add_subplot(2,len(scan_ids),scan_ids.index(scan_id)+1+len(scan_ids))
            ax_fig_ALL_size_ip = fig_ALL_size.add_subplot(2,len(scan_ids),scan_ids.index(scan_id)+1)
            ax_fig_ALL_size_oop = fig_ALL_size.add_subplot(2,len(scan_ids),scan_ids.index(scan_id)+1+len(scan_ids))
            axs = [ax_fig_ALL_strain_ip,ax_fig_ALL_strain_oop,ax_fig_ALL_size_ip,ax_fig_ALL_size_oop]
            y_labels=[r'$\varepsilon_\parallel$  (%)',r'$\varepsilon_\perp$  (%)',r'$d_\parallel$ / nm',r'$d_\perp$ / nm']
            for ax in axs:
                if ax in axs[0::2]:
                    ax.set_title(label)
                if plot_pot_step:
                    if ax in axs[1::2]:
                        ax.set_xlabel('Time(s)', fontsize=15)
                else:
                    if ax in axs[1::2]:
                        if plot_vs_RHE:
                            ax.set_xlabel('E / V$_{RHE}$', fontsize=15)
                        else:
                            ax.set_xlabel('E vs (Ag/AgCl) / V', fontsize=15)
                    if scan_id==scan_ids[0]:
                        ax.set_ylabel(y_labels[axs.index(ax)],fontsize = 15)
        if overplot:
            ax1_all=fig_all.add_subplot(231)
            ax2_all=fig_all.add_subplot(232)
            ax3_all=fig_all.add_subplot(233)
            ax4_all=fig_all.add_subplot(234)
            ax5_all=fig_all.add_subplot(235)
            ax6_all=fig_all.add_subplot(236)
            axs_all=[ax1_all,ax2_all,ax3_all,ax4_all,ax5_all,ax6_all]
            y_labels_CV=[r'$\varepsilon_\parallel$  (%)',r'$\varepsilon_\perp$  (%)',r'Current Density / mAcm$^{-2}$',r'$d_\parallel$ / nm',r'$d_\perp$ / nm',r'Intensity / a.u.']
            for ax in axs_all:
                if plot_pot_step:
                    ax.set_xlabel('Time(s)', fontsize=15)
                else:
                    if plot_vs_RHE: 
                        # if ax==ax6_all:
                        ax.set_xlabel('E / V$_{RHE}$', fontsize=15)
                    else:
                        ax.set_xlabel('E vs (Ag/AgCl) / V', fontsize=15)
                    ax.set_ylabel(y_labels_CV[axs_all.index(ax)],fontsize = 15)
        else:
            axs_all = [fig_all.add_subplot(5,len(scan_ids),i) for i in range(scan_ids.index(scan_id)+1,scan_ids.index(scan_id)+1+3*5,3)]
            y_labels_CV=[r'j/ mAcm$^{-2}$',r'$\Delta\varepsilon_\parallel$  (%)',r'$\Delta\varepsilon_\perp$  (%)',r'$\Delta d_\parallel$ / nm',r'$\Delta d_\perp$ / nm']
            y_ticks = [[-10,-5,0,5,10,15],[-0.20,-0.15,-0.10,-0.05,0.00],[-0.4,-0.3,-0.2,-0.1,0.0],[-2.0,-1.5,-1.0,-0.5,0.0],[-2.0,-1.5,-1.0,-0.5,0.0]]
            for ax in axs_all:
                if scan_id == scan_ids[-1]:
                    ax.set_xlim([0.50,2.0])
                    ax.set_xticks([0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0])
                    ax.set_xticklabels([0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0])
                else:
                    ax.set_xlim([0.70,1.9])
                    ax.set_xticks([0.8,1.0,1.2,1.4,1.6,1.8])
                ax.set_yticks(y_ticks[axs_all.index(ax)])
                ax.set_yticklabels(y_ticks[axs_all.index(ax)])
                if plot_pot_step:
                    ax.set_xlabel('Time(s)', fontsize=15)
                else:
                    if plot_vs_RHE: 
                        if ax==axs_all[-1]:
                            ax.set_xlabel('E / V$_{RHE}$', fontsize=15)
                        else:
                            ax.set_xticklabels([])
                    else:
                        if ax==axs_all[-1]:
                            ax.set_xlabel('E vs (Ag/AgCl) / V', fontsize=15)
                        else:
                            ax.set_xticklabels([])
                    if scan_id == scan_ids[0]:
                        ax.set_ylabel(y_labels_CV[axs_all.index(ax)],fontsize = 16)
                    else:
                        ax.set_yticklabels([])
            ax3_all,ax1_all,ax2_all,ax4_all,ax5_all = axs_all

        scan_no = int(scan_id.split('_')[-1])
        beamtime = scan_id.split('_')[0]
        data = scan_info[scan_id].data
        data_dic = {}
        for key in data.keys():
            data_dic[key] = getattr(data.f, key)
        potential_original = data_dic['potential']
        for key in data_dic.keys():
            if not plot_pot_step and data_dic[key].shape == np.array(potential_original).shape:
                temp,data_dic[key] = select_cycle((potential_original, data_dic[key]),plot_mode='CV')
            elif plot_pot_step and data_dic[key].shape == np.array(potential_original).shape:
                temp,data_dic[key] = select_cycle((potential_original, data_dic[key]),plot_mode='pot_step')
        data = data_dic
        strain_ip, uncertainty_strain_ip = strain_ip_with_uncertainty(data['cen_ip'], scan_info[scan_id].HKL_position, scan_info[scan_id].structure_lattice, 0) #np.sqrt(scans[i]['pcov'][:][1, 1]
        strain_oop, uncertainty_strain_oop = strain_oop_with_uncertainty(data['cen_oop'], scan_info[scan_id].HKL_position, scan_info[scan_id].structure_lattice, 0) #, np.sqrt(scans[i]['pcov'][:][2, 2]))
        intensity = data['amp_ip']*data['amp_oop'] #gaussian_filter(scans[i]['intensity'], sigma=0)*1000 #sigma=20
        I = data['amp_ip']*data['FWHM_oop']*data['FWHM_ip']

        if(scan_info[scan_id].analog_data_filename):
            ana_data = analog_data(analog_data_path_single+scan_info[scan_id].analog_data_filename)
            ana_reflectivity = ana_data.reflectivity/ana_data.laser_intensity
            ana_reflectivity = gaussian_filter(ana_reflectivity, sigma=100)
            ana_reflectivity /= np.max(ana_reflectivity)
            ana_E = ana_data.potential
            ana_j = ana_data.current / (np.pi * 0.2**2)
            #if(scan_id == 'CH5314_056'):
            #    ana_j /= 10
            ana_j = gaussian_filter(ana_j, sigma=100)

        elif(scan_info[scan_id].analog_data_interval):
            ana_data = analog_data.in_interval(analog_data_path_cont, start_date=scan_info[scan_id].analog_data_interval[0], end_date=scan_info[scan_id].analog_data_interval[1])
            ana_reflectivity = ana_data.reflectivity/ana_data.laser_intensity
            ana_reflectivity = gaussian_filter(ana_reflectivity, sigma=100)
            ana_reflectivity /= np.max(ana_reflectivity)
            ana_E = ana_data.potential
            ana_j = ana_data.current / (np.pi * 0.2**2)
            ana_j = gaussian_filter(ana_j, sigma=100)
            ana_epoch = ana_data.epoch

        else:
            ana_reflectivity = []
            ana_E = []
            ana_j = []

        if(create_ASCII):
            np.savetxt('data/%s.dat'%(scan_id), np.array([strain_ip, strain_oop, data['FWHM_ip'], data['FWHM_oop'], data['amp_ip']]).T, '%.5e', ', ', '\r\n', header='Strain ip / %, Strain oop / %, FWHM ip / A^-1, FWHM oop / A^-1, Peak Intensity')

        scan_direction_ranges = scan_info[scan_id].scan_direction_ranges
        for pos in range(len(scan_direction_ranges)):
            if pos == len(scan_direction_ranges)-1:
                # draw red dot at the first position, if not desired put break
                break
                indx1 = 0
                indx2 = 1
                marker = 'o'
                fillcolor = 'r'

            else:
                indx1 = scan_direction_ranges[pos]
                indx2 = scan_direction_ranges[pos+1]+1
                if(1):
                    fillcolor = 'w' if pos%2 == 1 else scan_info[scan_id].color
                    marker = 'o'
                else:
                    fillcolor = scan_info[scan_id].color if pos%2 == 1 else 'r'
                    marker = 'o'

            potential = temp[indx1:indx2]
            if plot_pot_step:
                label=scan_info[scan_id].scan_label
                me_color ='None'
            else:
                if(pos == 0):
                    label=scan_info[scan_id].scan_label+' positive scan'
                    me_color='k'
                else:
                    label = scan_info[scan_id].scan_label+' negative scan'
                    me_color='r'
            if plot_pot_step:
                marker = None
                scale_factor =float(scan_time)/len(data['Time'])
                offset = 0
                if scan_id in ['P23config_zap_742_20','P23config_zap_742_100']:
                    offset = 1
                if scan_id in ['P23config_zap_746_20','P23config_zap_746_100']:
                    offset =2
                ax1a.plot(data['Time'][indx1:indx2]*scale_factor+offset, set_max_to_0(strain_ip,[indx1,indx2]), '-', linewidth=1,  markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label='')
                ax2a.plot(data['Time'][indx1:indx2]*scale_factor+offset, set_max_to_0(strain_oop,[indx1,indx2]), '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label='')
                ax3a.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(data['FWHM_ip'],[indx1,indx2]), '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax4a.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(data['FWHM_oop'],[indx1,indx2]), '-', linewidth=1,  markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax5a.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(data['amp_ip']/data['amp_ip'][0],[indx1,indx2]), '-', linewidth=1,  markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax6a.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(0.2*np.pi/data['FWHM_ip'],[indx1,indx2]), 'o-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax7a.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(0.2*np.pi/data['FWHM_oop'],[indx1,indx2])+0.25, '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label='')
                ax9a.plot(data['Time'][indx1:indx2]*scale_factor+offset, POT(data['potential'][indx1:indx2], plot_vs_RHE, scan_info[scan_id].pH), '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)

                ax1_all.plot(data['Time'][indx1:indx2]*scale_factor+offset, set_max_to_0(strain_ip,[indx1,indx2]), '-', linewidth=1,  markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax2_all.plot(data['Time'][indx1:indx2]*scale_factor+offset, set_max_to_0(strain_oop,[indx1,indx2]), '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax6_all.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(data['amp_ip']/data['amp_ip'][0],[indx1,indx2]), '-', linewidth=1,  markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax4_all.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(0.2*np.pi/data['FWHM_ip'],[indx1,indx2]), 'o-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax5_all.plot(data['Time'][indx1:indx2]*scale_factor+offset,set_max_to_0(0.2*np.pi/data['FWHM_oop'],[indx1,indx2]), '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax3_all.plot(data['Time'][indx1:indx2]*scale_factor+offset, POT(data['potential'][indx1:indx2], plot_vs_RHE, scan_info[scan_id].pH), '-', linewidth=1, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)

            else:
                ax1a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(strain_ip,[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax1_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(strain_ip,[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax1_all.set_ylim([-0.31,0.02])
                try:
                    ax_fig_ALL_strain_ip.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(strain_ip,[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                    ax_fig_ALL_strain_ip.set_ylim([-0.31,0.02])
                except:
                    pass
                ax2a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(strain_oop,[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax2_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(strain_oop,[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax2_all.set_ylim([-0.55,0.05])
                try:
                    ax_fig_ALL_strain_oop.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(strain_oop,[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                    ax_fig_ALL_strain_oop.set_ylim([-0.55,0.05])
                except:
                    pass
                ax3a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(data['FWHM_ip'],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax4a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(data['FWHM_oop'],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax5a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(data['amp_ip']/data['amp_ip'][0],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                try:
                    ax6_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(data['amp_ip']/data['amp_ip'][0],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                except:
                    pass
                ax6a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(0.2*np.pi/data['FWHM_ip'],[indx1,indx2]), 'o-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax4_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(0.2*np.pi/data['FWHM_ip'],[indx1,indx2]), 'o-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax4_all.set_ylim([-2.4,0.08])
                try:
                    ax_fig_ALL_size_ip.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(0.2*np.pi/data['FWHM_ip'],[indx1,indx2]), 'o-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                    ax_fig_ALL_size_ip.set_ylim([-2.4,0.08])
                except:
                    pass
                ax7a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(0.2*np.pi/data['FWHM_oop'],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax5_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(0.2*np.pi/data['FWHM_oop'],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                ax5_all.set_ylim([-2.,0.05])
                try:
                    ax_fig_ALL_size_oop.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH),set_max_to_0(0.2*np.pi/data['FWHM_oop'],[indx1,indx2]), '-', linewidth=1, color=scan_info[scan_id].color, markerfacecolor=fillcolor, markeredgecolor=me_color, marker=marker, label=label)
                    ax_fig_ALL_size_oop.set_ylim([-2.,0.05])
                except:
                    pass
                if scan_info[scan_id].ids_filename==None:
                    ax9a.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH), data['current_density'][indx1:indx2]*8.*50, '-', linewidth=2, color='r', markerfacecolor=fillcolor, markeredgecolor=me_color, marker=None)
                    ax3_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH), data['current_density'][indx1:indx2]*8.*50, '-', linewidth=2, color='r', markerfacecolor=fillcolor, markeredgecolor=me_color, marker=None)
                    ax3_all.plot(POT(potential, plot_vs_RHE, scan_info[scan_id].pH), data['current_density'][indx1:indx2]*8., '--', linewidth=2, color='r', markerfacecolor=fillcolor, markeredgecolor=me_color, marker=None)
                    ax3_all.set_ylim((-5,5))
                else:
                    cv_data = extract_ids_file(scan_info[scan_id].ids_filename)
                    if pos==0:
                        colors_lib = {8:'sienna',10:'red',13:'green'}
                        ax9a.plot(POT(cv_data[0], plot_vs_RHE, scan_info[scan_id].pH), cv_data[1]*1000*8*50, '-',color = colors_lib[scan_info[scan_id].pH],markerfacecolor=fillcolor, linewidth=2,label=scan_info[scan_id].scan_label, markeredgecolor=me_color, marker=None)
                        ax9a.set_ylim((-15,15))
                        ax9a.set_xlim((0.5,1.8))
                        ax9a.annotate(r'$\times$50',xy=(1.2,9.6),xytext=(1.1,7.68),color='r')
                    ax3_all.plot(POT(cv_data[0], plot_vs_RHE, scan_info[scan_id].pH), cv_data[1]*1000*8.*50, '-', label=scan_info[scan_id].scan_label, linewidth=2, color = 'r', markerfacecolor=fillcolor, markeredgecolor=me_color, marker=None)
                    ax3_all.set_ylim((-14,15))
                    ax3_all.annotate(r'$\times$50',xy=(1.2,9.6),xytext=(1.1,7.68),color='r')
                    ax3_all.set_title(scan_info[scan_id].scan_label,fontsize = 15)

        if not overplot:
            fig_all.subplots_adjust(wspace=0.02,hspace=0.02)
        lim_l = scan_info[scan_id].analog_data_plot_range[0]
        lim_r = scan_info[scan_id].analog_data_plot_range[1]

        for ax in [ax1a, ax2a, ax3a, ax4a, ax5a, ax6a, ax7a, ax8a, ax9a]:
            ax.legend(loc=0, frameon=False)

        header = '%s #%d CV with XRD\r\nPotential / V, Current Density / mA/cm^2, Strain ip / %%, Strain oop / %%, d ip / nm, d oop / nm, Intensity (area ip)'%(beamtime, scan_no)
        X = np.array([data['potential'], data['current_density'], strain_ip, strain_oop, (0.2*np.pi/data['FWHM_ip']), (0.2*np.pi/data['FWHM_oop']), data['amp_ip']]).T
        filename = 'data/ascii/%s_%d_CV_XRD.dat'%(beamtime, scan_no)
        np.savetxt(filename, X, newline='\r\n', header=header)

    if len(scan_ids) > 1:
        path = 'plots/many/'
    else:
        path =  'plots/' + scan_id + '/'
        if not os.path.exists(path):
            os.makedirs(path)

    fig1.savefig(path+'strain_ip.png', dpi=300, bbox_inches='tight')
    fig2.savefig(path+'strain_oop.png', dpi=300, bbox_inches='tight')
    fig3.savefig(path+'FWHM_ip.png', dpi=300, bbox_inches='tight')
    fig4.savefig(path+'FWHM_oop.png', dpi=300, bbox_inches='tight')
    fig5.savefig(path+'Peak_Intensity.png', dpi=300, bbox_inches='tight')
    fig6.savefig(path+'Grainsize_ip.png', dpi=300, bbox_inches='tight')
    fig7.savefig(path+'Grainsize_oop.png', dpi=300, bbox_inches='tight')
    fig9.savefig(path+'Pot_CurrentDensity.png', dpi=300, bbox_inches='tight')
    fig_all.savefig(path+'Pot_all_together.png', dpi=300, bbox_inches='tight')
    try:
        fig_ALL_strain.savefig(path+'Pot_all_strain_together.png', dpi=300, bbox_inches='tight')
        fig_ALL_size.savefig(path+'Pot_all_size_together.png', dpi=300, bbox_inches='tight')
    except:
        pass

plt.show()
