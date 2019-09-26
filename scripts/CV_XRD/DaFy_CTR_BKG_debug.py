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
        #make nice looking status bar
        finish_percent = (i+1)/float(img_loader.total_frame_number)
        column_size = get_console_size()[0]-22
        left_abnormal = int((img_loader.abnormal_range[0]+1)/float(img_loader.total_frame_number)*column_size+8)
        right_abnormal = int((img_loader.abnormal_range[1]+1)/float(img_loader.total_frame_number)*column_size+8)
        output_text =list('{}{}{}{}{}'.format('BEGIN(0)','='*int(finish_percent*column_size),'==>',' '*int((1-finish_percent)*column_size),'>|END('+str(img_loader.total_frame_number)+')'))
        for index_text in range(len(output_text)):
            if output_text[index_text]!=' ' and index_text>left_abnormal and index_text<right_abnormal:
                output_text[index_text] = 'x'
            else:
                pass
        print(''.join(output_text),end="\r")


        ii = 0
        #if not rod_scan:
        #    ii = 1

        def make_tweak_string():
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
            return tweak_motion_str

        def tweak_integration(integration_object, tweak_motion_str, pre_tweak_motion):
            repeat_last = ''
            tweak = ''
            process_through = False
            tweak_return = integration_object.update_var_for_tweak(tweak_motion_str)
            if tweak_return in ['process_through','qw','rm']:
                tweak = False
            else:
                tweak = True
            if tweak_return == 'tweak':
                repeat_last = False
            elif tweak_return == 'repeat_last':
                repeat_last = True
            if tweak_return != 'repeat_last':
                #all parameters are updated in previous step, so you just 'qw' to repeat. 
                pre_tweak_motion = tweak_motion_str
            else:
                pre_tweak_motion = pre_tweak_motion
            if tweak_return == 'process_through':
                process_through = True
            return integration_object, tweak, tweak_return, repeat_last, pre_tweak_motion, process_through

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
                if repeat_last:
                    tweak_return = bkg_sub.update_var_for_tweak(pre_tweak_motion)
                else:
                    tweak_motion_str = make_tweak_string()
                    print('tweak_motion_str:',tweak_motion_str)
                    all_return = tweak_integration(bkg_sub,tweak_motion_str,pre_tweak_motion)
                    bkg_sub, tweak, tweak_return, repeat_last, pre_tweak_motion, process_through = all_return 
                    print('pre_tweak_motion:',pre_tweak_motion)
                check_result = bkg_sub.fit_background(None,img, data, freeze_sf = True)
                data = update_data_bkg(data, bkg_sub) 
                fig = plot_bkg_fit(fig, data, bkg_sub)
                plt.pause(.05)
            else:
                tweak = False
            if tweak_return == 'rm':
                data = pop_last_item(data)

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
    
    fig = plt.figure(figsize=(9,6))
    fig = plot_bkg_fit(fig, data, bkg_sub)
    plt.pause(9.05)

    df = pd.DataFrame(data)
    df.to_excel('test_ctr.xlsx')

if __name__ == "__main__":
    run()
