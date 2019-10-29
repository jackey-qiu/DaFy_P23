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
#make compatibility of py 2 and py 3#
if (sys.version_info > (3, 0)):
    raw_input = input

class run_app(object):
    def __init__(self):
        self.stop = True
        #self.conf_file_names = {'CTR':'config_p23_ctr_new.ini'}
        #config files
        #self.conf_file = os.path.join(DaFy_path, 'config', self.conf_file_names['CTR'])
        self.conf_file = None

    def run(self, config = None):
        #extract global vars from config
        if config == None:
            #self.kwarg_global = extract_vars_from_config(self.conf_file, section_var ='Global')
            pass
        else:
            self.conf_file = config
            self.kwarg_global = extract_vars_from_config(self.conf_file, section_var ='Global')
        for each in self.kwarg_global:
            setattr(self,each,self.kwarg_global[each])

        #pars lib for everything else
        self.kwarg_visulization = extract_vars_from_config(self.conf_file, section_var ='Visulization')
        self.kwarg_data = extract_vars_from_config(self.conf_file, section_var ='Data_Storage')
        self.kwarg_image = extract_vars_from_config(self.conf_file, section_var = 'Image_Loader')
        self.kwarg_mask = extract_vars_from_config(self.conf_file,section_var = 'Mask')

        #recal clip_boundary and cen
        self.clip_boundary = {"ver":[self.cen[0]-self.clip_width['ver'],self.cen[0]+self.clip_width['ver']+1],
                        "hor":[self.cen[1]-self.clip_width['hor'],self.cen[1]+self.clip_width['hor']+1]}     
        self.cen_clip = [self.clip_width['ver'],self.clip_width['hor']]

        self.img = None
        #data file
        self.data = {}
        for key in self.data_keys:
            self.data[key]=[]
        # print(data)
        #init peak fit, bkg subtraction and reciprocal space and image loader instance
        self.bkg_sub = background_subtraction_single_img(self.cen_clip, self.conf_file, sections = ['Background_Subtraction'])
        self.img_loader = nexus_image_loader(clip_boundary = self.clip_boundary, kwarg = self.kwarg_image)
        self.create_mask_new = create_mask(kwarg = self.kwarg_mask)
        self.setup_frames()

    def setup_frames(self):
        #build generator funcs
        self._scans = scan_generator(scans = self.scan_nos)
        self._images = image_generator_bkg(self._scans,self.img_loader,self.create_mask_new)

    def run_script(self):
        try:
            img = next(self._images)
            self.img = img
            self.data = merge_data_image_loader(self.data, self.img_loader)
            self.bkg_sub.fit_background(None, img, self.data, plot_live = True, freeze_sf = True)
            self.data = merge_data_bkg(self.data, self.bkg_sub)
            return True
        except:
            pass
            #return False

    def run_update(self):
        self.bkg_sub.fit_background(None, self.img, self.data, plot_live = True, freeze_sf = True)
        self.data = update_data_bkg(self.data, self.bkg_sub)

if __name__ == "__main__":
    run_app()
