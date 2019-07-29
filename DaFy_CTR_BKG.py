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
from FitEnginePool import Reciprocal_Space_Mapping, XRD_Peak_Fitting,background_subtraction_single_img
from GrainAnalysisEnginePool import cal_strain_and_grain
from VisualizationEnginePool import show_3_plots, plot_after_fit
from DataFilterPool import create_mask_bkg, extract_subset_of_zap_scan
import numpy as np
from util.UtilityFunctions import nexus_image_loader,get_UB,gen_find,edf_image_loader, update_data, pop_last_item, tiff_image_loader
from scipy.interpolate import griddata
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter
import pandas as pd
try:
    from mpi4py import MPI
except:
    print('mpi4py not installed, use single processor!')

#make compatibility of py 2 and py 3#
if (sys.version_info > (3, 0)):
    raw_input = input

#timer#
t0 = time.time()
mpi_run = False
if mpi_run:
    comm=MPI.COMM_WORLD
    size_cluster=comm.Get_size()
    rank=comm.Get_rank()
    if rank == 0:
        print('Starting MPI run with {} processors in total!'.format(size_cluster))            
else:
    size_cluster = 1
    rank = 0
##########tweak_mode to manually change the peak center!############
tweak_mode = False

#parse all global config items to local variables
#You may need to edit the config file accordingly(eg, scan_no)
#Choose config_id03_template.ini for id3 data and config_p23_template.ini for p23 dataset
#'CV_XRD' config is the config for running this program
#'ploter' config is the file which will be written with config info for later plotting step
#'bkg_sub' config is the config file for background substraction

conf_file_names = {'XRD':'config_p23_i20180678.ini',\
                   'ploter':'CV_XRD_plot_i20180678_Jul29_2019.ini'}
conf_file = os.path.join(DaFy_path, 'config', conf_file_names['XRD'])
conf_file_plot = os.path.join(DaFy_path, 'config', conf_file_names['ploter'])
config = configparser.ConfigParser()
config.read(conf_file)
for section in config.sections():
    for each in config.items(section):
        globals()[each[0]] = eval(each[1])

if mpi_run:
    tweak_mode, live_image, debug = False, False, False

#fig handles for live images in tweak mode
plt.ion()
fig_bkg_plot = plt.figure(figsize=(9,6))

for scan_no in scan_nos:
    if rank == 0:
        print('---------------Working on scan_{} dataset now!-------------'.format(scan_no))
    if name == "P23":
        spec_file = os.path.join(spec_file_head,'{}_{:0>5}.fio'.format(frame_prefix,scan_no))
        img_path = img_path_head
    elif name in ['13IDC','ID03']:
        spec_file = spec_file
        img_path = img_path
    #reciprocal space instance of skin_layer#
    # lattice_skin = rsp.lattice.from_cif(os.path.join(DaFy_path, 'util','cif',"{}".format(film_material_cif)), HKL_normal=[1,1,1],HKL_para_x=[1,1,-2], offset_angle=0)

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
        #spec = Fio.FioFile(fio_path = spec_file)
        spec = None
    elif name == 'ID03':
        spec = Fio.SpecFile(spec_file)
    elif name == '13IDC':
        spec = Fio.SpecFile_APS(spec_file)

    ###########image loader############
    if name == 'P23':
        img_loader = nexus_image_loader(fio_path = spec_file, nexus_path = img_path, frame_prefix = frame_prefix, scan_number = scan_no)
    elif name == 'ID03':
        img_loader = edf_image_loader(spec_file, img_path, is_zap_scan)
    elif name == '13IDC':
        img_loader = tiff_image_loader(spec_file,img_path)
    ############init Recip mapping instance##########

    ###############init a fit###########
    ##step 1: load image###
    if name == 'P23':
        img = img_loader.load_frame_new(frame_number = 0)
        size = img.shape
    elif name == 'ID03':
        img = img_loader.load_frame(scan_no, frame_no = 0).img
        size = img.shape
    elif name == '13IDC':
        img = img_loader.load_frame(scan_no, frame_no = 0)
        size = img.shape

    ##step 2: get motor positions##
    if name == 'P23':
        motor_angles = img_loader.extract_motor_angles(0, constant_motors =constant_motors)
    elif name == 'ID03':
        motor_angles = spec.extract_motor_angle(scan_no, 0, is_zap_scan)
    elif name == '13IDC':
        motor_angles = spec.extract_motor_angle(scan_no, 0)

    ############background subtraction#####
    bkg_sub = background_subtraction_single_img(conf_file, sections = ['Integration_setup','Correction_pars','Spec_info'])
    # geo_info = spec.extract_header_info(scan_no)
    # bkg_sub.update_geom(geo_info)
    bkg_sub.update_motor_angles(motor_angles)

    ##step 3: do recip space mapping##
    #get grid_indices for gridded q#

    ##step 4: set mask and init peak fitting engine##
    mask = create_mask_bkg(img = img,\
                       threshold = 10000, compare_method = 'larger', remove_columns = remove_columns, remove_rows = remove_rows, remove_pix = remove_pix,\
                       remove_xy_range = {'par':remove_xy_par,'ver':remove_xy_ver}, \
                       remove_partial_range = {'point_couple':line_strike_segments, 'pixel_width':line_strike_width})

    #how many frames to be processed?
    def find_boundary(n_process,n_jobs,rank):
        step_len=int(n_jobs/n_process)
        remainder=int(n_jobs%n_process)
        left,right=0,0
        if rank<=remainder-1:
            left=rank*(step_len+1)
            right=(rank+1)*(step_len+1)-1
        elif rank>remainder-1:
            left=remainder*(step_len+1)+(rank-remainder)*step_len
            right=remainder*(step_len+1)+(rank-remainder+1)*step_len-1
        return left,right

    if size_cluster == 1:
        frame_number_ranges = range(img_loader.total_frame_number)
    else:
        lf, rg = find_boundary(size_cluster, img_loader.total_frame_number,rank)
        frame_number_ranges = range(lf, rg)


    ##Step 5: loop through all frames##
    frame_to_be_rechecked = []
    if size_cluster > 1:
        comm.Barrier()
    t1 = time.time()
    if debug:
        frame_number_ranges = [debug_img]
    process_through = False

    for frame_number in frame_number_ranges:
        if rank == 0:
            print("Work on frame_{} now...".format(frame_number))
        #load image
        try:
            img = img_loader.load_frame_new(frame_number)
        except:
            break
        
        t2_a = time.time()
        if name == 'ID03':
            img = img.img
        #recip mapping
        #launch peak fitting
        if name == 'P23':
            #motor_angles = spec.extract_motor_angle(motor_angles, scan_no, frame_number, updated_angles =['mu','delta','gamma','omega_t'])
            motor_angles = img_loader.extract_motor_angles(frame_number, constant_motors =constant_motors)
        elif name == 'ID03':
            motor_angles = spec.extract_motor_angle(scan_no, frame_number, is_zap_scan)
            
        # img = img/motor_angles['mon']/motor_angles['transm']
        
        bkg_sub.update_motor_angles(motor_angles)
        #Now do backgrond subtraction for the peak
        if name in ['13IDC']:
            L = spec.extract_L(scan_no, frame_number)
        elif name == 'ID03':
            L = spec.extract_L(scan_no, frame_number,is_zap_scan)
        elif name == 'P23':
            H,K,L = img_loader.extract_HKL(frame_number)
            potential, current = img_loader.extract_pot_current(frame_number)
        data = update_data(data,keys = ['potential', 'current_density','H','K','L', 'frame_number', 'phi', 'chi','mu','delta','gamma','omega_t','mon','transm'], \
                           values =[potential, current, H, K, L, frame_number,motor_angles['phi'],motor_angles['chi'],motor_angles['mu'],\
                           motor_angles['delta'],motor_angles['gamma'],motor_angles['omega_t'],motor_angles['mon'], motor_angles['transm']])

        tweak = True
        pre_tweak_motion = ''
        repeat_last = False
        while tweak:
            print(pre_tweak_motion)
            if bkg_sub.int_direct =='y':
                check_result = bkg_sub.fit_background(fig_bkg_plot,img*mask, data, plot_live = plot_live, check=True,check_level = check_level)
            elif bkg_sub.int_direct =='x':
                check_result = bkg_sub.fit_background(fig_bkg_plot,img*mask, data, plot_live = plot_live, check=True, check_level = check_level)
            if process_through:
                tweak = False
                data['peak_intensity'].append(bkg_sub.fit_results['I'])
                data['peak_intensity_error'].append(bkg_sub.fit_results['Ierr'])
            else:
                if not repeat_last:
                    tweak_motion_str = raw_input(", splited stream of string\n\
                                                  ud:up or down\n\
                                                  lr:left or right\n\
                                                  cw:column width\n\
                                                  rw:row width\n\
                                                  pw:peak width\n\
                                                  ps:peak shift\n\
                                                  od:polynomial order\n\
                                                  sf:ss_factor, smaller lower bkg\n\
                                                  r:repeat last motion\n\
                                                  #r:repeat motion for rest points\n\
                                                  ft:fit function(ah, sh, stq, atq)\n\
                                                  qw:quit and write date\n\
                                                  rm:remove current date and quit\n\
                                                  Your input is:") or 'qw'
                    if tweak_motion_str =='#r':
                        tweak_return = 'process_through'
                    else:
                        tweak_return = bkg_sub.update_var_for_tweak(tweak_motion_str)
                else:
                    tweak_return = bkg_sub.update_var_for_tweak(pre_tweak_motion)

                if tweak_return == 'qw':
                    data['peak_intensity'].append(bkg_sub.fit_results['I'])
                    data['peak_intensity_error'].append(bkg_sub.fit_results['Ierr'])
                    tweak = False
                    # repeat_last = False
                elif tweak_return == 'rm':
                    data = pop_last_item(data,keys = ['potential', 'current_density','H','K','L', 'frame_number', 'phi', 'chi','mu','delta','gamma','omega_t','mon','transm'])
                    tweak = False
                    # repeat_last = False
                elif tweak_return == 'tweak':
                    if tweak_motion_str!='r':
                        pre_tweak_motion = tweak_motion_str
                    tweak = True
                    repeat_last = False
                elif tweak_return == 'repeat_last':
                    repeat_last = True
                    teak = True
                elif tweak_return == 'process_through':
                    data['peak_intensity'].append(bkg_sub.fit_results['I'])
                    data['peak_intensity_error'].append(bkg_sub.fit_results['Ierr'])
                    tweak = False
                    process_through = True


                #save results for the first loop but update result in the following loop steps used in tweak mode
                #extract potential and current density
            fig_bkg_plot.canvas.draw()
            fig_bkg_plot.tight_layout()
            plt.pause(0.05)
            plt.show()

    if debug or len(scan_nos)==1:
        plt.ioff()
        plt.show()
    
    print(rank, 'Fit is completed with time elapse of {} seconds!\nAnd it takes {} seconds to set up. Cheers:)'.format(time.time()-t0,t1-t0))
    if size_cluster>1:
        comm.Barrier()
        data = comm.gather(data,root = 0)
        if rank == 0:
            data_list = list(data)
            data_complete = {}
            for each_key in data[0]:
                data_complete[each_key]=[]
            for each_data in data_list:
                for each_key in each_data:
                    data_complete[each_key] = list(data_complete[each_key])+list(each_data[each_key])
            data = data_complete

    if rank == 0:
        #Now append info to the plotting config file
        with open(conf_file_plot, 'a+') as f:
            f.write('\n[{}]\n'.format(scan_id))
            f.write("scan_number = [{}]\n".format(scan_no))
            f.write("HKL = ({:3.1f} {:3.1f} {:3.1f})\n".format(data['H'][0],data['K'][0],data['L'][0]))
            f.write("rod_scan ={}\n".format(rod_scan))
            f.write("colors = ['r']\n")
            f.write("scan_ids = ['DaFy_{}']\n".format(scan_no))
            f.write("data_file_header = '{}/data'\n".format(DaFy_path))
            f.write("data_files= ['{}.npz']\n".format(scan_id))
        if not debug:
            print('Fit is completed with time elapse of {} seconds!\nCheers:)'.format(time.time()-t0))
        #plot results at the end for single dataset only
        #save data
        np.savez('data/%s.npz'%(scan_id),frame_number=data['frame_number'], potential = data['potential'], current = data['current_density'], H=data['H'], K=data['K'], L=data['L'], \
                 peak_intensity = data['peak_intensity'], peak_intensity_error = data['peak_intensity_error'],\
                 phi = data['phi'], chi = data['chi'], mu=data['mu'], delta = data['delta'],\
                 gamma=data['gamma'],omega_t =data['omega_t'],mon=data['mon'],transm=data['transm'])
    else:
        pass
