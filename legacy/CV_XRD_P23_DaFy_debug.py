import sys,os
sys.path.append('./EnginePool/')
sys.path.append('./FilterPool/')
sys.path.append('./util/')
import time, configparser, ntpath, subprocess
from util import Fio
from FitEnginePool import Reciprocal_Space_Mapping, XRD_Peak_Fitting
from VisualizationEnginePool import show_3_plots, plot_after_fit
from DataFilterPool import create_mask
import matplotlib
matplotlib.use("tkAgg")
import matplotlib.pyplot as plt
import numpy as np
from util.UtilityFunctions import nexus_image_loader,get_UB,gen_find,edf_image_loader
from scipy.interpolate import griddata
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter

t0 = time.time()
#debug setup
debug = 0
debug_img = 0
# L_scan = False
live_image = 1 
if debug:
    live_image = True

#parse all glocal config items to local variables
#You may need to edit the config file accordingly
conf_file = "/home/qiu/apps/DaFy/config_id03.ini"
config = configparser.ConfigParser()
config.read(conf_file)
for section in config.sections():
    for each in config.items(section):
        globals()[each[0]] = eval(each[1])

#data container
data = {}
for key in data_keys:
    data[key]=[]
#id for data saving
scan_id =  'DataBank_{}_{}'.format(scan_no,time.strftime("%Y%m%d-%H%M%S"))

#spec file reader
if name == 'P23':
    spec = Fio.FioFile(fio_path = spec_file)
elif name == 'ID03':
    spec = Fio.SpecFile(spec_file)

#calculate UB matrix
if name == 'P23':
    UB = get_UB(name)(lattice_constants,energy,or0_angles,or1_angles,or0_hkl,or1_hkl)
elif name == 'ID03':
    UB = get_UB(name)(spec_file, scan_no)

#image loader 
if name == 'P23':
    img_loader = nexus_image_loader(fio_path = spec_file, nexus_path = img_path, frame_prefix = frame_prefix)
elif name == 'ID03':
    img_loader = edf_image_loader(spec_file, img_path, is_zap_scan)

#Recip mapping instance
r_map = Reciprocal_Space_Mapping(img = None, E_keV = energy, cen = cen, pixelsize = pix_size, sdd = sdd, UB = UB)

###init a fit###
if name == 'P23':
    img = img_loader.load_frame(scan_no, frame_number = 0)
    size = img.shape
elif name == 'ID03':
    img = img_loader.load_frame(scan_no, frame_no = 0).img
    size = img.shape
###recip mapping###
#get motor angles first#
if name == 'P23':
    motor_angles = spec.extract_motor_angle(motor_angles, scan_no, 0, updated_angles =['mu','delta','gamma','omega_t'])
elif name == 'ID03':
    motor_angles = spec.extract_motor_angle(scan_no, 0, is_zap_scan)

#at this momen, UB is calculated at the beginning, in the future it should be able to parsed from SPEC or FIO file.
r_map.update_img(img, UB = None, motor_angles = motor_angles)
q_oop, q_ip, grid_q_oop, grid_q_ip, grid_intensity_init = r_map.q['q_perp'], r_map.q['q_par'],\
                                                     r_map.q['grid_q_perp'], r_map.q['grid_q_par'], r_map.grid_intensity
# plt.figure()
# plt.pcolormesh(grid_q_ip,grid_q_oop,grid_intensity_init)
# plt.show()
# quit()

#get grid_indices#                            
grid_indices = griddata((q_ip.ravel(), q_oop.ravel()), np.arange(size[0]*size[1]).ravel(), (grid_q_ip, grid_q_oop), method='nearest')
print(guess_ip,guess_oop)
#mask and init peak fitting engine#
mask = create_mask(grid_intensity_init, threshold = 10000, compare_method = 'larger', remove_columns = remove_columns, remove_rows = remove_rows, remove_pix = remove_pix)
peak_fit_engine = XRD_Peak_Fitting(img = grid_intensity_init, mask = mask, q_ip = grid_q_ip, q_oop=grid_q_oop,\
                                   cut_offset=cut_offset,\
                                   data_range_offset = data_range_offset,\
                                   peak_center = [cen[1],cen[0]],\
                                   fit_bounds = {'hor': bounds_ip, 'ver': bounds_oop},\
                                   fit_p0 = {'hor': guess_ip, 'ver': guess_oop})
#Do you want to show live images
if live_image:
    plt.ion()
    fig_3_plot = plt.figure(figsize=(15,5))
 
frame_number_ranges = range(img_loader.total_frame_number)
if debug:
    frame_number_ranges = [debug_img]

for frame_number in frame_number_ranges:
    print("Fitting frame_{} now...".format(frame_number))
    #load image
    img = img_loader.load_frame(scan_no, frame_number)
    if name == 'ID03':
        img = img.img
    #recip mapping
    if l_scan:#this step is time consuming for recal q, but only for L_scan
        motor_angles = spec.extract_motor_angle(motor_angles, scan_no, frame_number, updated_angles =['mu','delta','gamma','omega_t'])
        r_map.update_img(img, UB = None, motor_angles = motor_angles,update_q = True)
        q_oop, q_ip, grid_q_oop, grid_q_ip, img = r_map.q['q_perp'], r_map.q['q_par'],\
                                                  r_map.q['grid_q_perp'], r_map.q['grid_q_par'], r_map.grid_intensity
        grid_indices = griddata((q_ip.ravel(), q_oop.ravel()), np.arange(size[0]*size[1]).ravel(), (grid_q_ip, grid_q_oop), method='nearest')
    else:#if not L_scan no need to recal q, then set update_q to False
        r_map.update_img(img, UB = None, motor_angles = None, update_q = False)
        img = r_map.intensity
        flat_int = img.ravel()
        img = flat_int[grid_indices]
        img = img.reshape(size)
    #reset mask
    mask = create_mask(img, threshold = 10000, compare_method = 'larger', remove_columns = remove_columns, remove_rows = remove_rows, remove_pix = remove_pix)
    peak_fit_engine.mask = mask
    #launch peak fitting
    peak_fit_engine.reset_fit(img)
    #save results
    data = peak_fit_engine.save_data(data)
    #extract potential and current density
    if name == 'P23':
        potential, current = spec.extract_pot_current(scan_no, frame_number)
    elif name == 'ID03':
        potential, current = spec.extract_pot_current(scan_no, frame_number,is_zap_scan)

    data['potential'].append(potential)
    data['current_density'].append(current)
    print(peak_fit_engine.fit_results['hor'][0][0],peak_fit_engine.fit_results['ver'][0][0])
    if live_image:
        fig_3_plot = show_3_plots(fig_3_plot,img*mask,grid_q_ip,grid_q_oop,vmin = vmin_2d_img, vmax =vmax_2d_img, cmap = 'jet', \
                                  fit_data = peak_fit_engine.fit_data, model=peak_fit_engine.model, fit_results=peak_fit_engine.fit_results,\
                                  cut_offset = peak_fit_engine.cut_offset, peak_center = peak_fit_engine.peak_center, title = 'Frame_{}, E ={:04.2f}V'.format(frame_number,potential))
        fig_3_plot.canvas.draw()
        fig_3_plot.tight_layout()
        plt.pause(0.05)
        plt.show()
plt.ioff()
plt.show()
if not debug:
    print('Fit is completed with time elapse of {} seconds!\nCheers:)'.format(time.time()-t0))
    plot_after_fit(data,is_zap_scan)
    #save data
    np.savez('data/%s.npz'%(scan_id), potential=data['potential'], current_density=data['current_density'], Time=data['Time'], pcov_ip=data['pcov_ip'], pcov_oop=data['pcov_oop'],
             cen_ip=data['cen_ip'], FWHM_ip=data['FWHM_ip'], amp_ip=data['amp_ip'], lorfact_ip=data['lfrac_ip'], bg_slope_ip=data['bg_slope_ip'], bg_offset_ip=data['bg_offset_ip'],
             cen_oop=data['cen_oop'], FWHM_oop=data['FWHM_oop'], amp_oop=data['amp_oop'], lorfact_oop=data['lfrac_oop'], bg_slope_oop=data['bg_slope_oop'], bg_offset_oop=data['bg_offset_oop'])
