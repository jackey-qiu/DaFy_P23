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

def run():
    #make compatibility of py 2 and py 3#
    if (sys.version_info > (3, 0)):
        raw_input = input

    #config files
    conf_file_names = {'CTR':'config_p23_ctr_new.ini'}
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
    process_through = False
    for img in _images:
        if img_loader.scan_number!=scan_number:
            i = 0
            scan_number = img_loader.scan_number
        data = merge_data_image_loader(data, img_loader)
        bkg_sub.fit_background(None, img, data, plot_live = True, freeze_sf = True)
        data = merge_data_bkg(data, bkg_sub)
        # fig = plot_bkg_fit(fig, data, bkg_sub)
        # fig.canvas.draw()
        # fig.tight_layout()
        # plt.pause(0.05)
        # plt.show()
        i = i+1
        finish_percent = (i+1)/float(img_loader.total_frame_number)
        column_size = get_console_size()[0]-22
        print('{}{}{}{}{}'.format('BEGIN(0)','='*int(finish_percent*column_size),'=->',' '*int((1-finish_percent)*column_size),'>|END('+str(img_loader.total_frame_number)+')'),end="\r")
        time.sleep(0.003)

        tweak = True
        if not live_image:
            tweak = False
        pre_tweak_motion = ''
        repeat_last = False
        ii = 0
        #if not rod_scan:
        #    ii = 1
        while tweak:
            print('pre_tweak_motion:',pre_tweak_motion)
            # print('ss factor before:',bkg_sub.ss_factor)
            if not process_through:
                fig = plot_bkg_fit(fig, data, bkg_sub)
                fig.canvas.draw()
                fig.tight_layout()
                plt.pause(0.05)
                plt.show()
            if update_width:
                check_result = bkg_sub.fit_background(fig, img, data, plot_live = True, freeze_sf = True)
                if ii ==0:
                    temp_c_width= bkg_sub.find_peak_width(img, initial_c_width=400, initial_r_width = None, direction = 'vertical')
                    if temp_c_width!=None:
                        if temp_c_width-bkg_sub.col_width>50 and frame_number_ranges.index(frame_number)>10:#bragg peak area, not very robust
                            bkg_sub.col_width = bkg_sub.col_width
                        else:
                            bkg_sub.col_width = temp_c_width

                    check_result = bkg_sub.fit_background(fig, img, data, plot_live = True, freeze_sf = True)
                ii  = ii+1
            # else:
                # check_result = bkg_sub.fit_background(fig,img, data, plot_live =True, freeze_sf = True)
            if process_through:
                tweak = False
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
                    tweak = False
                elif tweak_return == 'rm':
                    data = pop_last_item(data,keys = ['potential', 'current','H','K','L', 'image_no', 'phi', 'chi','mu','delta','gamma','omega_t','mon','transm'])
                    tweak = False
                elif tweak_return == 'tweak':
                    if tweak_motion_str!='r':
                        pre_tweak_motion = tweak_motion_str
                    tweak = True
                    repeat_last = False
                elif tweak_return == 'repeat_last':
                    repeat_last = True
                    teak = True
                elif tweak_return == 'process_through':
                    tweak = False
                    process_through = True
            # print('ss factor after:',bkg_sub.ss_factor)
            check_result = bkg_sub.fit_background(None,img, data, plot_live =True, freeze_sf = True)
            data = update_data_bkg(data, bkg_sub)
            if not process_through:
                fig = plot_bkg_fit(fig, data, bkg_sub)
                fig.canvas.draw()
                fig.tight_layout()
                plt.pause(.05)
                plt.show()
    fig = plt.figure(figsize=(9,6))
    fig = plot_bkg_fit(fig, data, bkg_sub)
    fig.canvas.draw()
    fig.tight_layout()
    plt.pause(9.05)
    plt.show()
    df = pd.DataFrame(data)
    df.to_excel('test_ctr.xlsx')

if __name__ == "__main__":
    run()
