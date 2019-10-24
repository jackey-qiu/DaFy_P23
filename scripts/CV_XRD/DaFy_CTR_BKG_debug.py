import sys,os,itertools
import numpy as np
import pandas as pd
try:
    from . import locate_path
except:
    import locate_path
try:
    from mpi4py import MPI
except:
    print('mpi4py not installed, use single processor!')
try:
    import ConfigParser as configparser
except:
    import configparser
script_path = locate_path.module_path_locator()
DaFy_path = os.path.dirname(os.path.dirname(script_path))
sys.path.append(DaFy_path)
sys.path.append(os.path.join(DaFy_path,'EnginePool'))
sys.path.append(os.path.join(DaFy_path,'FilterPool'))
sys.path.append(os.path.join(DaFy_path,'util'))
import matplotlib,time
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from GrainAnalysisEnginePool import cal_strain_and_grain
from VisualizationEnginePool import plot_bkg_fit
from DataFilterPool import create_mask, merge_data_bkg, update_data_bkg, merge_data_image_loader, make_data_config_file
from FitEnginePool import fit_pot_profile
from FitEnginePool import Reciprocal_Space_Mapping
from FitEnginePool import XRD_Peak_Fitting
from FitEnginePool import background_subtraction_single_img
from util.XRD_tools import reciprocal_space_v3 as rsp
from util.UtilityFunctions import pop_last_item
from util.UtilityFunctions import image_generator_bkg
from util.UtilityFunctions import scan_generator
from util.UtilityFunctions import nexus_image_loader
from util.UtilityFunctions import find_boundary
from util.UtilityFunctions import extract_global_vars_from_string
from util.UtilityFunctions import extract_vars_from_config
from util.UtilityFunctions import get_console_size
from util.UtilityFunctions import make_tweak_string
from util.UtilityFunctions import tweak_integration

def run():
    #make compatibility of py 2 and py 3#
    if (sys.version_info > (3, 0)):
        raw_input = input

    #config files
    conf_file_names = {'CTR':'config_p23_ctr_new_I20180678.ini'}
    conf_file = os.path.join(DaFy_path, 'config', conf_file_names['CTR'])

    #extract global vars from config
    kwarg_global = extract_vars_from_config(conf_file, section_var ='Global')
    for each in kwarg_global:
        globals()[each] = kwarg_global[each]

    #pars lib for everything else
    kwarg_visulization = extract_vars_from_config(conf_file, section_var ='Visulization')
    kwarg_data = extract_vars_from_config(conf_file, section_var ='Data_Storage')
    kwarg_image = extract_vars_from_config(conf_file, section_var = 'Image_Loader')
    kwarg_mask = extract_vars_from_config(conf_file,section_var = 'Mask')

    #recal clip_boundary and cen
    clip_boundary = {"ver":[cen[0]-clip_width['ver'],cen[0]+clip_width['ver']+1],
                     "hor":[cen[1]-clip_width['hor'],cen[1]+clip_width['hor']+1]}     
    cen_clip = [clip_width['ver'],clip_width['hor']]
    
    if live_image:
        plt.ion()
        fig = plt.figure(figsize=(9,6))

    #data file
    data = {}
    for key in data_keys:
        data[key]=[]
    # print(data)
    #init peak fit, bkg subtraction and reciprocal space and image loader instance
    bkg_sub = background_subtraction_single_img(cen_clip, conf_file, sections = ['Background_Subtraction'])
    img_loader = nexus_image_loader(clip_boundary = clip_boundary, kwarg = kwarg_image)
    create_mask_new = create_mask(kwarg = kwarg_mask)

    #build generator funcs
    _scans = scan_generator(scans = scan_nos)
    _images = image_generator_bkg(_scans,img_loader,create_mask_new)

    i = 0
    scan_number = img_loader.scan_number
    scan_numbers_all = []
    process_through = False

    for img in _images:
        if img_loader.scan_number!=scan_number:
            i = 0
            scan_number = img_loader.scan_number
            scan_numbers_all.append(scan_number)
        data = merge_data_image_loader(data, img_loader)
        if update_width:
            _ = bkg_sub.find_peak_width(img, img_no = img_loader.frame_number, initial_c_width=400, initial_r_width = 400)
        else:
            if i == 0:
                _ = bkg_sub.find_peak_width(img, img_no = img_loader.frame_number, initial_c_width=400, initial_r_width = 400)
        bkg_sub.fit_background(None, img, data, plot_live = True, freeze_sf = True)
        data = merge_data_bkg(data, bkg_sub)
        i = i+1

        #make nice looking status bar
        finish_percent = (i+1)/float(img_loader.total_frame_number)
        column_size = int(get_console_size()[0])-22
        output_text =list('{}{}{}{}{}'.format('BEGIN(0)','='*int(finish_percent*column_size),'==>',' '*int((1-finish_percent)*column_size),'>|END('+str(img_loader.total_frame_number)+')'))
        print(''.join(output_text),end="\r")

        pre_tweak_motion = 'qw'
        repeat_last = False
        if live_image:
            tweak = True
        else:
            tweak = False

        while tweak:
            if not process_through:
                print('pre_tweak_motion:',pre_tweak_motion)
                fig = plot_bkg_fit(fig, data, bkg_sub)
                plt.pause(.05)
                tweak_motion_str = make_tweak_string()
                all_return = tweak_integration(bkg_sub,tweak_motion_str,pre_tweak_motion)
                bkg_sub, tweak, tweak_return, repeat_last, pre_tweak_motion, process_through = all_return
                check_result = bkg_sub.fit_background(None,img, data, freeze_sf = True)
                data = update_data_bkg(data, bkg_sub)
                fig = plot_bkg_fit(fig, data, bkg_sub)
                plt.pause(.05)
                if tweak_return == 'rm':
                    data = pop_last_item(data)
            else:
                tweak = False

    df = pd.DataFrame(data)
    scan_labels='CTR_data_'
    for scan in scan_numbers_all:
        if scan == scan_numbers_all[-1]:
            scan_labels = scan_labels + 'scan{}'.format(scan)
        else:
            scan_labels = scan_labels + 'scan{}_'.format(scan)
    df.to_excel(os.path.join(DaFy_path,'data','{}.xlsx'.format(scan_labels)))

    fig = plt.figure(figsize=(9,6))
    fig = plot_bkg_fit(fig, df, bkg_sub, True)
    fig.savefig(os.path.join(DaFy_path,'temp','temp_ctr{}.png'.format(scan_labels)),dpi = 300)
    plt.pause(9.05)
    return df

if __name__ == "__main__":
    run()
