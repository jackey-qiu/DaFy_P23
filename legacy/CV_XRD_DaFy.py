import sys,os
import locate_path
DaFy_path = locate_path.module_path_locator()
sys.path.append(os.path.join(DaFy_path,'EnginePool'))
sys.path.append(os.path.join(DaFy_path,'FilterPool'))
sys.path.append(os.path.join(DaFy_path,'util'))
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import time,ntpath, subprocess
try:
    import ConfigParser as configparser
except:
    import configparser
from util import Fio
import itertools
from util.XRD_tools import reciprocal_space_v3 as rsp
from FitEnginePool import Reciprocal_Space_Mapping, XRD_Peak_Fitting
from GrainAnalysisEnginePool import cal_strain_and_grain
from VisualizationEnginePool import show_3_plots, plot_after_fit
from DataFilterPool import create_mask
import numpy as np
from util.UtilityFunctions import nexus_image_loader,get_UB,gen_find,edf_image_loader
from scipy.interpolate import griddata
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter
import pandas as pd
#make compatibility of py 2 and py 3#
if (sys.version_info > (3, 0)):
    raw_input = input

#timer#
t0 = time.time()
##########tweak_mode is not working stably now, so keep it off!########
tweak_mode = False

#parse all global config items to local variables
#You may need to edit the config file accordingly(eg, scan_no)
#Choose config_id03.ini for id3 data and config_p23.ini for p23 dataset
conf_file_names = [ "config_p23.ini","config_id03_MA4171.ini","CV_XRD_plot_temp.ini",\
                    "config_id03_ch5314.ini","config_id03_ma3886.ini",\
                    "config_p23_CoOOH.ini","config_p23_CoOOH_Prep3.ini"]
conf_file = os.path.join(DaFy_path, 'config', conf_file_names[-2])
conf_file_plot = os.path.join(DaFy_path, 'config', conf_file_names[2])
config = configparser.ConfigParser()
config.read(conf_file)
for section in config.sections():
    for each in config.items(section):
        globals()[each[0]] = eval(each[1])

for scan_no in scan_nos:
    print('---------------Working on scan_{} dataset now!-------------'.format(scan_no))
    if name == "P23":
        spec_file = os.path.join(spec_file_head,'{}_{:0>5}.fio'.format(frame_prefix,scan_no))
        img_path = os.path.join(img_path_head,'{}_{:0>5}/lmbd'.format(frame_prefix,scan_no))
    elif name == "ID03":
        spec_file = spec_file
        img_path = img_path
    #reciprocal space instance of skin_layer#
    # lattice_skin = rsp.lattice.from_cif(os.path.join(DaFy_path, 'util','cif',"{}".format(film_material_cif)), HKL_normal=[1,1,1],HKL_para_x=[1,1,-2], offset_angle=0)
    lattice_skin = rsp.lattice.from_cif(os.path.join(DaFy_path, 'util','cif',"{}".format(film_material_cif)), HKL_normal=film_hkl_normal,HKL_para_x=film_hkl_x, offset_angle=0)

    ###########data container########
    data = {}
    for key in data_keys:
        data[key]=[]
    if 'peak_intensity' not in data.keys():
        data['peak_intensity'] = []
    #id for data saving
    scan_id =  'DataBank_{}_{}'.format(scan_no,time.strftime("%Y%m%d-%H%M%S"))

    ###########spec file reader#######
    if name == 'P23':
        spec = Fio.FioFile(fio_path = spec_file)
    elif name == 'ID03':
        spec = Fio.SpecFile(spec_file)

    ###########calculate UB matrix#####
    if name == 'P23':
        UB = get_UB(name)(lattice_constants,energy,or0_angles,or1_angles,or0_hkl,or1_hkl)
    elif name == 'ID03':
        UB = get_UB(name)(spec_file, scan_no)

    ###########image loader############
    if name == 'P23':
        img_loader = nexus_image_loader(fio_path = spec_file, nexus_path = img_path, frame_prefix = frame_prefix)
    elif name == 'ID03':
        img_loader = edf_image_loader(spec_file, img_path, is_zap_scan)

    ############init Recip mapping instance##########
    r_map = Reciprocal_Space_Mapping(img = None, E_keV = energy, cen = cen, pixelsize = pix_size, sdd = sdd, UB = UB, boost_mapping=boost_mapping)

    ###############init a fit###########
    ##step 1: load image###
    if name == 'P23':
        img = img_loader.load_frame(scan_no, frame_number = 0)
        size = img.shape
    elif name == 'ID03':
        img = img_loader.load_frame(scan_no, frame_no = 0).img
        size = img.shape

    ##step 2: get motor positions##
    if name == 'P23':
        motor_angles = spec.extract_motor_angle(motor_angles, scan_no, 0, updated_angles =['mu','delta','gamma','omega_t'])
    elif name == 'ID03':
        motor_angles = spec.extract_motor_angle(scan_no, 0, is_zap_scan)

    ##step 3: do recip space mapping##
    r_map.update_img(img, UB = None, motor_angles = motor_angles)
    q_oop, q_ip, grid_q_oop, grid_q_ip, grid_intensity_init = r_map.q['q_perp'], r_map.q['q_par'],\
                                                         r_map.q['grid_q_perp'], r_map.q['grid_q_par'], r_map.grid_intensity
    #get grid_indices for gridded q#
    grid_indices = griddata((q_ip.ravel(), q_oop.ravel()), np.arange(size[0]*size[1]).ravel(), (grid_q_ip, grid_q_oop), method='nearest')

    ##step 4: set mask and init peak fitting engine##
    mask = create_mask(img = grid_intensity_init,img_q_par = grid_q_oop, img_q_ver = grid_q_ip,\
                       threshold = 1000, compare_method = 'larger', remove_columns = remove_columns, remove_rows = remove_rows, remove_pix = remove_pix,\
                       remove_q_range = {'par':remove_q_par,'ver':remove_q_ver}, \
                       remove_partial_range = {'point_couple':line_strike_segments, 'pixel_width':line_strike_width})
    try:
        peak_fit_engine = XRD_Peak_Fitting(img = grid_intensity_init, mask = mask, q_ip = grid_q_ip, q_oop=grid_q_oop,\
                                           cut_offset=cut_offset,\
                                           data_range_offset = data_range_offset,\
                                           peak_center = [cen[1],cen[0]],\
                                           fit_bounds = {'hor': bounds_ip, 'ver': bounds_oop},\
                                           fit_p0 = {'hor': guess_ip, 'ver': guess_oop})
        # plt.show()
        # quit()
    except:
        plt.pcolormesh(grid_q_ip, grid_q_oop, grid_intensity_init*mask)
        plt.show()
        print('Initialization of Peak fitting failed.\n Please check the config file (peak center and boundary).')
        time.sleep(10)
    #Do you want to show live images
    if live_image:
        plt.ion()
        fig_3_plot = plt.figure(figsize=(12,10))
    #how many frames to be processed?
    frame_number_ranges = range(img_loader.total_frame_number)
    if debug:
        frame_number_ranges = [debug_img]

    ##Step 5: loop through all frames##
    frame_to_be_rechecked = []
    for frame_number in frame_number_ranges:
        print("Fitting frame_{} now...".format(frame_number))
        #load image
        try:
            img = img_loader.load_frame(scan_no, frame_number)
        except:
            break
        if name == 'ID03':
            img = img.img
        #recip mapping
        if l_scan:#this step is time consuming for recal q, but only for L_scan
            if name == 'P23':
                motor_angles = spec.extract_motor_angle(motor_angles, scan_no, frame_number, updated_angles =['mu','delta','gamma','omega_t'])
            elif name == 'ID03':
                motor_angles = spec.extract_motor_angle(scan_no, frame_number, is_zap_scan)
            # motor_angles = spec.extract_motor_angle(motor_angles, scan_no, frame_number, updated_angles =['mu','delta','gamma','omega_t'])
            r_map.update_img(img, UB = None, motor_angles = motor_angles,update_q = True)
            q_oop, q_ip, grid_q_oop, grid_q_ip, img = r_map.q['q_perp'], r_map.q['q_par'],\
                                                      r_map.q['grid_q_perp'], r_map.q['grid_q_par'], r_map.grid_intensity
        else:#if not L_scan no need to recal q, then set update_q to False
            r_map.update_img(img, UB = None, motor_angles = None, update_q = False)
            img = r_map.intensity
            flat_int = img.ravel()
            img = flat_int[grid_indices]
            img = img.reshape(size)
        #launch peak fitting
        for i in range(1000):
            check_result = peak_fit_engine.reset_fit(img, use_first_fit_for_pos, check = False, tweak_mode = bool(i))
            #print out the hor and ver post for quick eye check
            print(peak_fit_engine.fit_results['hor'][0][0],peak_fit_engine.fit_results['ver'][0][0])
            #save results for the first loop but update result in the following loop steps used in tweak mode
            if i==0:
                #save results
                data = peak_fit_engine.save_data(data)
                #extract potential and current density
                if name == 'P23':
                    try:
                        potential, current = spec.extract_pot_current(scan_no, frame_number)
                    except:
                        potential, current = data['potential'][-1]*2-data['potential'][-2], data['current_density'][-1]
                elif name == 'ID03':
                    potential, current = spec.extract_pot_current(scan_no, frame_number,is_zap_scan)
                data['potential'].append(potential)
                data['current_density'].append(current)
            else:
                data = peak_fit_engine.update_data(data)
            #calculate stain and crystallite size
            data = cal_strain_and_grain(data,HKL = hkl, lattice = lattice_skin)
            #plot live image if you want
            if live_image:
                fig_3_plot = show_3_plots(fig_3_plot,img,grid_q_ip,grid_q_oop,vmin = vmin_2d_img, vmax =vmax_2d_img, cmap = 'jet', is_zap_scan = is_zap_scan,\
                                        fit_data = peak_fit_engine.fit_data, model=peak_fit_engine.model, fit_results=peak_fit_engine.fit_results_plot,processed_data_container=data,\
                                        cut_offset = peak_fit_engine.cut_offset, peak_center = peak_fit_engine.peak_center, title = 'Frame_{}, E ={:04.2f}V'.format(frame_number,potential))
                fig_3_plot.canvas.draw()
                fig_3_plot.tight_layout()
                plt.pause(0.05)
                plt.show()

            #launch tweak mode(you have 1000 times to tweak for each data point)
            if tweak_mode:
                center_offset = raw_input("left(l)/right(r)/up(u)/down(d) + pix_number.\
                                          \ne.g. l5 move center towards left by 5 pixes.\
                                          \n Enter to quit the looping!\
                                          \nYour input is:") or 'q'
                opt = {"l":"-", "r":"+", "u":"-", "d":"+"}
                if center_offset[0] in ["l", "r"]:
                    peak_fit_engine.peak_center[1] = peak_fit_engine.peak_center[1] + \
                                                    eval(center_offset.replace(center_offset[0],opt[center_offset[0]]))
                elif center_offset[0] in ["u", "d"]:
                    peak_fit_engine.peak_center[0] = peak_fit_engine.peak_center[0] + \
                                                    eval(center_offset.replace(center_offset[0],opt[center_offset[0]]))
                else:
                    break
            else:
                break
    plt.ioff()
    plt.show()
    print("Please recheck following frames:{}".format(frame_to_be_rechecked))
    #Now append info to the plotting config file
    with open(conf_file_plot, 'a+') as f:
        f.write('\n[{}]\n'.format(scan_id))
        f.write("scan_number = [{}]\n".format(scan_no))
        f.write("phs = [{}]\n".format(ph))
        f.write("colors = ['r']\n")
        f.write("xtal_lattices = ['CoHO2_9009884']\n")
        f.write("scan_ids = ['DaFy_{}']\n".format(scan_no))
        f.write("scan_labels = ['pH {}']\n".format(ph))
        f.write("ids_file_header ='/home/qiu/data/beamtime/P23_11_18_I20180114/ids/I20180014'\n")
        f.write("ids_files = ['??.ids']\n")
        f.write("data_file_header = '/home/qiu/apps/DaFy/data/save_data'\n")
        f.write("data_files= ['{}.npz']\n".format(scan_id))
        f.write("hkls = [[1,0,7]]\n")
        f.write("scan_direction_ranges =[[0,111,-2]]\n")
        f.write("plot_pot_steps = [{}]\n".format(int(is_zap_scan)))
    if not debug:
        print('Fit is completed with time elapse of {} seconds!\nCheers:)'.format(time.time()-t0))
        #plot results at the end for single dataset only
        if len(scan_nos)==1:
            plot_after_fit(data,is_zap_scan)
        #save data
        np.savez('data/%s.npz'%(scan_id),potential=data['potential'], \
                current_density=data['current_density'], Time=data['Time'], pcov_ip=data['pcov_ip'], pcov_oop=data['pcov_oop'],\
                cen_ip=data['cen_ip'], FWHM_ip=data['FWHM_ip'], amp_ip=data['amp_ip'], lorfact_ip=data['lfrac_ip'],\
                bg_slope_ip=data['bg_slope_ip'], bg_offset_ip=data['bg_offset_ip'], cen_oop=data['cen_oop'], \
                FWHM_oop=data['FWHM_oop'], amp_oop=data['amp_oop'], lorfact_oop=data['lfrac_oop'], \
                bg_slope_oop=data['bg_slope_oop'], bg_offset_oop=data['bg_offset_oop'])
