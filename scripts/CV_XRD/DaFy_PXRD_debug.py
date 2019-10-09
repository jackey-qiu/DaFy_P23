import sys,os,itertools
import numpy as np
import pandas as pd
import copy
import scipy
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
from VisualizationEnginePool import plot_pxrd_profile,plot_pxrd_profile_time_scan
from DataFilterPool import create_mask, merge_data_bkg, update_data_bkg, merge_data_image_loader, make_data_config_file
from FitEnginePool import fit_pot_profile,backcor
from util.XRD_tools import reciprocal_space_v3 as rsp
from util.UtilityFunctions import image_generator_bkg
from util.UtilityFunctions import scan_generator
from util.UtilityFunctions import nexus_image_loader
from util.UtilityFunctions import extract_vars_from_config
from util.UtilityFunctions import show_status_bar_2

def run():
    #make compatibility of py 2 and py 3#
    if (sys.version_info > (3, 0)):
        raw_input = input

    #config files
    conf_file_names = {'CTR':'config_p23_pxrd_new.ini'}
    conf_file = os.path.join(DaFy_path, 'config', conf_file_names['CTR'])

    #extract global vars from config
    kwarg_global = extract_vars_from_config(conf_file, section_var ='Global')
    for each in kwarg_global:
        globals()[each] = kwarg_global[each]
    if time_scan:
        plot_pxrd = plot_pxrd_profile_time_scan
    else:
        plot_pxrd = plot_pxrd_profile

    #pars lib for everything else
    kwarg_image = extract_vars_from_config(conf_file, section_var = 'Image_Loader')
    kwarg_mask = extract_vars_from_config(conf_file,section_var = 'Mask')
    kwarg_bkg = extract_vars_from_config(conf_file,section_var = 'Background_Subtraction')

    #recal clip_boundary and cen(you need to remove the edges)
    ver_offset = clip_width['ver']
    hor_offset = clip_width['hor']
    clip_boundary = {"ver":[ver_offset,dim_detector[0]-ver_offset],"hor":[hor_offset,dim_detector[1]-hor_offset]}
    cen_clip = [cen[0]-ver_offset,cen[1]-hor_offset]

    #init peak fit, bkg subtraction and reciprocal space and image loader instance
    img_loader = nexus_image_loader(clip_boundary = clip_boundary, kwarg = kwarg_image)
    create_mask_new = create_mask(kwarg = kwarg_mask)

    #live image or not
    plt.ion()
    fig = plt.figure(figsize=(15,8))

    #build generator funcs
    _scans = scan_generator(scans = scan_nos)
    _images = image_generator_bkg(_scans,img_loader,create_mask_new)

    scan_number = img_loader.scan_number#scan_number = None at the very beginning
    int_intensity = {}
    for img in _images:
        if img_loader.scan_number not in int_intensity:
            int_intensity[img_loader.scan_number] = {'2theta':[],'intensity':[]}
            if time_scan:
                int_intensity[img_loader.scan_number]['frame_number']=[]
                int_intensity[img_loader.scan_number]['potential']=[]
                for i in range(len(delta_segment_time_scan)):
                    int_intensity[img_loader.scan_number]['intensity_peak{}'.format(i+1)]=[]

        if scan_number == None:
            scan_number = img_loader.scan_number
        else:
            if img_loader.scan_number!=scan_number:
                #save data
                int_intensity_pd = pd.DataFrame(int_intensity[scan_number]).sort_values(by = '2theta')
                #smooth the intensity
                inten_smoothed = scipy.signal.savgol_filter(int_intensity_pd['intensity'],101,3)
                int_intensity_pd['intensity_smoothed']=inten_smoothed
                int_intensity_pd.to_excel(os.path.join(DaFy_path,'data','pxrd_scan{}.xlsx'.format(scan_number)))
                fig = plot_pxrd(fig,int_intensity[scan_number],None,plot_final = True)
                scan_number = img_loader.scan_number

        delta = img_loader.motor_angles['delta']
        delta_range = list(np.round(delta + np.arctan((cen[0] - np.array(range(dim_detector[0]-ver_offset*2)))*ps/sd)/np.pi*180, 3))
        if not(min(delta_range)>delta_range_plot[1] or max(delta_range)<delta_range_plot[0]):
            int_range = np.sum(img[:,:],axis = 1)

            #bkg subtraction
            n=np.array(range(len(int_range)))
            #by default the contamination rate is 25%
            #the algorithm may fail if the peak cover >40% of the cut profile
            bkg_n = int(len(n)/4)
            y_sorted = list(copy.deepcopy(int_range))
            y_sorted.sort()
            std_bkg =np.array(y_sorted[0:bkg_n*3]).std()/(max(y_sorted)-min(y_sorted))
            int_range[np.argmin(int_range)] =  int_range[np.argmin(int_range)+1]
            int_range_bkg, *discard = backcor(range(len(int_range)),int_range,\
                                              ord_cus=kwarg_bkg['ord_cus_s'],\
                                              s=std_bkg*kwarg_bkg['ss_factor'],fct=kwarg_bkg['fct'])
            int_range = list(int_range-int_range_bkg)
            int_range_bkg = list(int_range_bkg)

            #append results
            for j in delta_range:
                if j in int_intensity[img_loader.scan_number]['2theta']:
                    jj = int_intensity[img_loader.scan_number]['2theta'].index(j)
                    int_intensity[img_loader.scan_number]['intensity'][jj] = 0.5*int_intensity[img_loader.scan_number]['intensity'][jj] + 0.5*int_range[delta_range.index(j)]
                else:
                    int_intensity[img_loader.scan_number]['2theta'].append(j)
                    int_intensity[img_loader.scan_number]['intensity'].append(int_range[delta_range.index(j)])

            #append results for timescan
            if time_scan:
                k=0
                for each_segment in delta_segment_time_scan:
                    k = k+1
                    index_left = np.argmin(np.abs(np.array(delta_range)-each_segment[0]))
                    index_right = np.argmin(np.abs(np.array(delta_range)-each_segment[1]))
                    int_intensity[img_loader.scan_number]['intensity_peak{}'.format(k)].append(np.array(int_range)[min([index_left,index_right]):max([index_left,index_right])].sum())
                int_intensity[img_loader.scan_number]['potential'].append(img_loader.current)
                int_intensity[img_loader.scan_number]['frame_number'].append(img_loader.frame_number)
        #make nice looking status bar
        show_status_bar_2(img_loader, column_size_offset = 22)
        #plot after each frame
        if live_image:
            fig = plot_pxrd(fig,int_intensity[img_loader.scan_number],img,delta_range,int_range, int_range_bkg)
    #save data
    int_intensity_pd = pd.DataFrame(int_intensity[scan_number]).sort_values(by = '2theta')
    #smooth the intensity
    inten_smoothed = scipy.signal.savgol_filter(int_intensity_pd['intensity'],101,3)
    int_intensity_pd['intensity_smoothed']=inten_smoothed
    int_intensity_pd.to_excel(os.path.join(DaFy_path,'data','pxrd_scan{}.xlsx'.format(scan_number)))
    #plot in the end
    fig = plot_pxrd(fig,int_intensity[scan_number],None,plot_final = True)

if __name__ == "__main__":
    run()
