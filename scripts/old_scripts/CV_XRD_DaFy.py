import sys,os,itertools
import numpy as np
import pandas as pd
try:
    from . import locate_path
except:
    import locate_path
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
from VisualizationEnginePool import show_all_plots_new
from DataFilterPool import create_mask, merge_data, make_data_config_file
from FitEnginePool import fit_pot_profile
from FitEnginePool import Reciprocal_Space_Mapping
from FitEnginePool import XRD_Peak_Fitting
from FitEnginePool import background_subtraction_single_img
from util.XRD_tools import reciprocal_space_v3 as rsp
from util.UtilityFunctions import image_generator
from util.UtilityFunctions import scan_generator
from util.UtilityFunctions import nexus_image_loader
from util.UtilityFunctions import find_boundary
from util.UtilityFunctions import extract_global_vars_from_string
from util.UtilityFunctions import extract_vars_from_config
from util.UtilityFunctions import show_status_bar

def run():
    #make compatibility of py 2 and py 3#
    if (sys.version_info > (3, 0)):
        raw_input = input

    print("CV_XRD is running now ...")
    #config files
    conf_file_names = {'CV_XRD':'config_p23_template_new.ini'}

    conf_file = os.path.join(DaFy_path, 'config', conf_file_names['CV_XRD'])
    #conf_file_plot = os.path.join(DaFy_path, 'config', conf_file_names['ploter'])
    #conf_file_mpi = os.path.join(DaFy_path, 'config', conf_file_names['mpi'])

    #extract global vars from config
    kwarg_global = extract_vars_from_config(conf_file, section_var ='Global')
    for each in kwarg_global:
        globals()[each] = kwarg_global[each]

    #pars lib for everything else
    kwarg_visulization = extract_vars_from_config(conf_file, section_var ='Visulization')
    kwarg_film = extract_vars_from_config(conf_file, section_var ='Film_Lattice')
    kwarg_data = extract_vars_from_config(conf_file, section_var ='Data_Storage')
    kwarg_peak_fit = extract_vars_from_config(conf_file, section_var = 'Peak_Fit')
    kwarg_rsp = extract_vars_from_config(conf_file, section_var = 'Reciprocal_Mapping')
    kwarg_image = extract_vars_from_config(conf_file, section_var = 'Image_Loader')
    kwarg_mask = extract_vars_from_config(conf_file,section_var = 'Mask')

    #recal clip_boundary and cen
    clip_boundary = {"ver":[cen[0]-clip_width['ver'],cen[0]+clip_width['ver']+1],
                     "hor":[cen[1]-clip_width['hor'],cen[1]+clip_width['hor']+1]}
    cen_clip = [clip_width['ver'],clip_width['hor']]
    plt.ion()
    fig = plt.figure(figsize=(14,10))

    #data file
    data = {}
    for key in data_keys:
        data[key]=[]

    #init peak fit, bkg subtraction and reciprocal space and image loader instance
    peak_fitting_instance = XRD_Peak_Fitting(img = None, cen=cen_clip, kwarg = kwarg_peak_fit)
    bkg_sub_instance = background_subtraction_single_img(cen_clip, conf_file, sections = ['Background_Subtraction'])
    rsp_instance = Reciprocal_Space_Mapping(img =None, cen=cen_clip, kwarg = kwarg_rsp)
    img_loader = nexus_image_loader(clip_boundary = clip_boundary, kwarg = kwarg_image)
    create_mask_new = create_mask(kwarg = kwarg_mask)
    lattice_skin = rsp.lattice.from_cif(os.path.join(DaFy_path, 'util','cif',"{}".format(kwarg_film['film_material_cif'])),
                                        HKL_normal=kwarg_film['film_hkl_normal'],\
                                        HKL_para_x=kwarg_film['film_hkl_x'],\
                                        E_keV=rsp_instance.e_kev, offset_angle=0)

    #build generator funcs
    _scans = scan_generator(scans = scan_nos)
    _images = image_generator(_scans,img_loader,rsp_instance,peak_fitting_instance,create_mask_new)
    i = 0
    scan_number = img_loader.scan_number
    for img in _images:
        if img_loader.scan_number!=scan_number:
            i = 0
            scan_number = img_loader.scan_number
        peak_fitting_instance.reset_fit(img, check = False)
        bkg_sub_instance.fit_background(None, img, plot_live = False, freeze_sf = True)
        data = merge_data(data, img_loader, peak_fitting_instance, bkg_sub_instance, kwarg_global)
        data = cal_strain_and_grain(data,HKL = kwarg_film['film_hkl_bragg_peak'][0], lattice = lattice_skin)

        #make nice looking status bar
        show_status_bar(img_loader,column_size_offset = 22)

        i=i+1
        if live_image:
            fig = show_all_plots_new(fig=fig,fit_engine_instance=peak_fitting_instance,\
                                    bkg_int=bkg_sub_instance,\
                                    processed_data_container=data,\
                                    title='Frame_{}, E ={:04.2f}V'.format(i,data['potential'][-1]),\
                                    kwarg =kwarg_visulization)
            fig.canvas.draw()
            fig.tight_layout()
            plt.pause(0.05)
            plt.show()
    df = pd.DataFrame(data)
    df.to_csv('temp_current_result.csv')
    df.to_excel('temp_current_result.xlsx',columns =kwarg_global['data_keys'])
    make_data_config_file(
                          DaFy_path = DaFy_path,
                          data_folder = os.path.join(DaFy_path,'data'),
                          data = df,
                          film_material_cif=kwarg_film['film_material_cif'],
                          hkls=kwarg_film['film_hkl_bragg_peak'],
                          pot_step = kwarg_peak_fit['pot_step_scan'],
                          beamline=name,
                          beamtime=beamtime,
                          kwarg = kwarg_data
                          )
    if not live_image:
        fig = show_all_plots_new(fig=fig,fit_engine_instance=peak_fitting_instance,\
                                bkg_int=bkg_sub_instance,\
                                processed_data_container=data,\
                                title='Frame_{}, E ={:04.2f}V'.format(i,data['potential'][-1]),\
                                kwarg =kwarg_visulization)
        fig.canvas.draw()
        fig.tight_layout()
        plt.pause(500.05)
        plt.show()
        print('\nRun is finished!')

if __name__ == "__main__":
    run()

