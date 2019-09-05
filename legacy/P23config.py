# Settings and definitions for MA3886/CH5314
import sys
sys.path.append('/home/qiu/apps/SXRD/XRD_tools/')
import reciprocal_space_v3 as rsp
import numpy as np
# from pyspec import spec
from Fio import Fiofile
import os
import p23_tools_debug as p23
# import id03_tools as id03

class P23config:
    sdd = 700 # Sample-detector distance / mm
    #note cenx is index of horizontal direction towards right
    #ceny is index of vertical direction towards bottom
    #The corresponding num array index is actually [ceny,cenx]
    cenx = 280#247
    ceny = 624#745
    # cen = (cenx, ceny)
    E_keV = 18.739
    pix_size = (0.055, 0.055)

    data_path = '/home/qiu/data/beamtime/P23_11_18_I20180114/raw'
    # spec_filename = data_path + 'startup/FirstTest_00713.fio'
    # img_path = data_path + 'FirstTest_00713/lmbd'
    # spec_file = Fiofile(spec_filename)
    # DI = p23.DetectorImg(E_keV, cen, (0.055, 0.055), sdd, spec_filename=spec_filename,edf_path = img_path/raw)
    cont_labview_path = '/home/reikowski/data/2018_05_MA3886/Continuous Acquisition_ma3886/'
    single_labview_path = '/home/reikowski/data/2018_05_MA3886/CV_ma3886/'

    # data_path2 = '/home/qiu/data/CH5314/'
    # spec_filename2 = data_path2 + 'ch5314_sixcvertical.spec'
    # img_path2 = data_path2 + 'ch5314_img/'
    # spec_file2 = spec.SpecDataFile(spec_filename2)
    # DI2 = id03.DetectorImg(E_keV, cen, (0.055, 0.055), sdd, spec_filename2, img_path2)
    # cont_labview_path2 = '/home/qiu/data/CH5314/Continuous Acquisition_ch5314/'
    # single_labview_path2 = '/home/qiu/data/CH5314/CV_ch5314/'

    ivium_path = data_path + 'ec/'
    Au111_lat = rsp.lattice(a=2.8837, b=2.8837, c=7.0636, alpha=90, beta=90, gamma=120, basis=[['Au', 0,0,0],['Au', 2./3.,1./3.,1./3.],['Au', 1./3.,2./3.,2./3.]], HKL_normal = [0,0,1], HKL_para_x = [1,0,0], offset_angle=0)
    Co3O4_lat = rsp.lattice.from_cif('/home/qiu/apps/SXRD/XRD_simulator/cif/Co/Co3O4.cif', HKL_normal = [1,1,1], HKL_para_x=[1,1,-2], offset_angle=0+15)
    CoOOH_lat = rsp.lattice.from_cif('/home/qiu/apps/SXRD/XRD_simulator/cif/Co/CoHO2_9009884.cif', HKL_normal = [0,0,1], HKL_para_x=[1,0,0], offset_angle=15)
