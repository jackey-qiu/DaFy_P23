import numpy as np
from Reciprocal_space_tools.HKLVlieg import Crystal, printPos, UBCalculator, VliegAngles, printPos_prim, vliegDiffracAngles
from nexusformat.nexus import *
import h5py
import fnmatch, os
import re,os
from scipy import misc
from PyMca5.PyMcaIO import specfilewrapper
# from pyspec import spec
from PyMca5.PyMcaIO import EdfFile
from math import sqrt
from skimage.feature import blob_log
import matplotlib.pyplot as plt
from numpy import linspace, loadtxt, ones, convolve
import pandas as pd
import collections
from random import randint
from itertools import count
izip = zip

def remove_abnormality(mon, left_offset,right_offset):
    mon = np.array(mon)
    max_index = np.argmax(mon)
    print('kick of data points from {} to {}'.format(max_index-left_offset, max_index+right_offset))
    return [i for i in range(len(mon)) if (i<(max_index - left_offset) or i>(max_index+right_offset))]



# 3. Lets define some use-case specific UDF(User Defined Functions)

def moving_average(data, window_size):
    """ Computes moving average using discrete linear convolution of two one dimensional sequences.
    Args:
    -----
            data (pandas.Series): independent variable
            window_size (int): rolling window size

    Returns:
    --------
            ndarray of linear convolution

    References:
    ------------
    [1] Wikipedia, "Convolution", http://en.wikipedia.org/wiki/Convolution.
    [2] API Reference: https://docs.scipy.org/doc/numpy/reference/generated/numpy.convolve.html

    """
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(data, window, 'same')


def explain_anomalies(y, window_size, sigma=1.0):
    """ Helps in exploring the anamolies using stationary standard deviation
    Args:
    -----
        y (pandas.Series): independent variable
        window_size (int): rolling window size
        sigma (int): value for standard deviation
    Returns:
    --------
        a dict (dict of 'standard_deviation': int, 'anomalies_dict': (index: value))
        containing information about the points indentified as anomalies
    """
    avg = moving_average(y, window_size).tolist()
    residual = y - avg
    # Calculate the variation in the distribution of the residual
    std = np.std(residual)
    return {'standard_deviation': round(std, 3),
            'anomalies_dict': collections.OrderedDict([(index, y_i) for
                                                       index, y_i, avg_i in izip(count(), y, avg)
              if (y_i > avg_i + (sigma*std)) | (y_i < avg_i - (sigma*std))])}


def explain_anomalies_rolling_std(y, window_size, sigma=1.0):
    """ Helps in exploring the anamolies using rolling standard deviation
    Args:
    -----
        y (pandas.Series): independent variable
        window_size (int): rolling window size
        sigma (int): value for standard deviation

    Returns:
    --------
        a dict (dict of 'standard_deviation': int, 'anomalies_dict': (index: value))
        containing information about the points indentified as anomalies
    """
    avg = moving_average(y, window_size)
    avg_list = avg.tolist()
    residual = y - avg
    # Calculate the variation in the distribution of the residual
    # testing_std = pd.rolling_std(residual, window_size)
    testing_std = pd.rolling(residual, window_size).std()
    testing_std_as_df = pd.DataFrame(testing_std)
    rolling_std = testing_std_as_df.replace(np.nan,
                                  testing_std_as_df.ix[window_size - 1]).round(3).iloc[:,0].tolist()
    std = np.std(residual)
    return {'stationary standard_deviation': round(std, 3),
            'anomalies_dict': collections.OrderedDict([(index, y_i)
                                                       for index, y_i, avg_i, rs_i in izip(count(),
                                                                                           y, avg_list, rolling_std)
              if (y_i > avg_i + (sigma * rs_i)) | (y_i < avg_i - (sigma * rs_i))])}


# This function is repsonsible for displaying how the function performs on the given dataset.
def plot_results(x, y, window_size, sigma_value=1,
                 text_xlabel="X Axis", text_ylabel="Y Axis", applying_rolling_std=False):
    """ Helps in generating the plot and flagging the anamolies.
        Supports both moving and stationary standard deviation. Use the 'applying_rolling_std' to switch
        between the two.
    Args:
    -----
        x (pandas.Series): dependent variable
        y (pandas.Series): independent variable
        window_size (int): rolling window size
        sigma_value (int): value for standard deviation
        text_xlabel (str): label for annotating the X Axis
        text_ylabel (str): label for annotatin the Y Axis
        applying_rolling_std (boolean): True/False for using rolling vs stationary standard deviation
    """
    plt.figure(figsize=(15, 8))
    plt.plot(x, y, "k.")
    y_av = moving_average(y, window_size)
    plt.plot(x, y_av, color='green')
    plt.xlim(0, 1000)
    plt.xlabel(text_xlabel)
    plt.ylabel(text_ylabel)

    # Query for the anomalies and plot the same
    events = {}
    if applying_rolling_std:
        events = explain_anomalies_rolling_std(y, window_size=window_size, sigma=sigma_value)
    else:
        events = explain_anomalies(y, window_size=window_size, sigma=sigma_value)

    x_anomaly = np.fromiter(events['anomalies_dict'].keys(), dtype=int, count=len(events['anomalies_dict']))
    y_anomaly = np.fromiter(events['anomalies_dict'].values(), dtype=float,
                                            count=len(events['anomalies_dict']))
    plt.plot(x_anomaly, y_anomaly, "r*", markersize=12)

    # add grid and lines and enable the plot
    plt.grid(True)
    plt.show()

def check_abnormality(data, mask, tolerance_factor = 3):
    data, mask = np.array(data), np.array(mask)
    last_accpt_pt = np.where(mask==1)[0][-1]
    #intrinsic_offset = np.array(len(data)- last_accpt_pt)
    data_filtered = data[np.where(mask==1)]
    if np.abs(data[-1]-data_filtered[-1])>tolerance_factor:
        return False
    else:
        return True

def check_abnormality_old(data, mask,pot, tolerance_factor = 3):
    data, mask, pot = np.array(data), np.array(mask), np.array(pot)
    last_accpt_pt = np.where(mask==1)[0][-1]
    pot_diff = 1*(np.diff(pot)<0)
    intrinsic_offset = sum(pot_diff[last_accpt_pt-1:])
    #intrinsic_offset = np.array(len(data)- last_accpt_pt)
    data_filtered = data[np.where(mask==1)]
    diff_mean = np.abs(np.diff(data_filtered)).mean()
    if abs(np.abs(data[-1]-data_filtered[-1])+diff_mean*intrinsic_offset)>diff_mean*tolerance_factor:
        return False
    else:
        return True


def peak_checker(img_path='/media/qiu/JACKEY/1704_APS_13IDC/images/sb1_32mM_CaCl2_Zr_1/S013/sb1_32mM_CaCl2_Zr_1_S013_00207.tif',max_sigma= 30, num_sigma = 10, threshold = 0.1):
    img = misc.imread(img_path)
    img[np.isnan(img)] = 0
    img=img/img.max()
    blobs_log = blob_log(img, max_sigma=max_sigma, num_sigma=num_sigma, threshold=threshold)
    # Compute radii in the 3rd column.
    blobs_log[:, 2] = blobs_log[:, 2] * sqrt(2)
    fig, ax = plt.subplots(1, 1, figsize=(9, 3), sharex=True, sharey=True)
    # ax = axes.ravel()
    ax.imshow(img, interpolation='nearest')
    for blob in blobs_log:
        y, x, r = blob
        print(x,y,r)
        c = plt.Circle((x, y), r, color='g', linewidth=2, fill=False)
        ax.add_patch(c)
    ax.set_axis_off()
    plt.tight_layout()
    plt.show()

def gen_find_old(filepat, top):
    #find all filenames in a dic tree that match a shell wildcard pattern
    files = []
    for path, dirlist, filelist  in os.walk(top):
        for name in fnmatch.filter(filelist, filepat):
            files.append(os.path.join(path,name))
    return files

def gen_find(filepat, top):
    #find all filenames in a dic tree that match a shell wildcard pattern
    files = [f for f in os.listdir(top) if re.match(filepat, f)]
    return files

def extract_arg(config, section, local_lib):
    arg_list=[]
    temp = dict(config.items(section))
    keys=temp.keys()
    keys.sort()
    for key in keys:
        try:
            arg_list.append(eval(temp[key],None,local_lib))
        except:
            arg_list.append(None)
    return arg_list

def collect_args(local_lib,tag):
    args_list = [None]*20
    for key in local_lib.keys():
        items_key = key.rsplit("_")
        if tag == items_key[-1]:
            args_list[int(items_key[-2][3:])-1]=local_lib[key]
    return [value for value in args_list if value!= None]

def rename(newname):
    def decorator(f):
        f.__name__ = newname
        return f
    return decorator

def update_data(data,keys,values):
    for key, value in zip(keys,values):
        data[key].append(value)
    return data

def pop_last_item(data,keys):
    for key in keys:
        data[key].pop()
    return data

def extract_potential(pt_no = 2000, time_step = [10, 50, 10], pot_step = [0.2, 0.5, 0.8]):
    potential_container= []
    frames_per_time = float(pt_no)/sum(time_step)
    for i in range(len(pot_step)):
        potential_container = potential_container + [pot_step[i]]*int(time_step[i]*frames_per_time)
    for i in range(abs(len(potential_container)-pt_no)):
        potential_container.pop()
    return potential_container

def cal_UB_p23(lattice_constants=[2.8837,2.8837,7.0636,90,90,120],energy=18.739,or0_angles=[0.4,15.4,22.43,-30.9,0.,0.],or1_angles=[0.4,7.61,13.63,-38.,0.,0.],or0_hkl=[1.0009,1.-1.0009,4.0359],or1_hkl=[0.0,-0.5045,2.5225]):
    substrate=Crystal(lattice_constants[0:3],lattice_constants[3:])
    # or0_angles=[0.4,15.4,22.43,-30.9,0.,0.]
    # or1_angles=[0.4,7.61,13.63,-38.,0.,0.]
    # or0_hkl=[1.0009,1.-1.0009,4.0359]
    # or1_hkl=[0.0,-0.5045,2.5225]
    # energy=18.739
    ub_substrate=UBCalculator(substrate,energy)
    or0_angles=np.deg2rad(or0_angles)
    or1_angles=np.deg2rad(or1_angles)
    or0_angles = vliegDiffracAngles(or0_angles)
    or1_angles = vliegDiffracAngles(or1_angles)
    ub_substrate.setPrimaryReflection(or0_angles,or0_hkl)
    ub_substrate.setSecondayReflection(or1_angles,or1_hkl)
    ub_substrate.calculateU()
    return ub_substrate.getUB()

def cal_UB_id03(spec_file_name, scan_no):
    scan = specfilewrapper.Specfile(spec_file_name).select('{}.1'.format(scan_no))
    return np.array(scan.header('G')[2].split(' ')[-9:],dtype = np.float)

def get_UB(name = 'P23'):
    if name in ['P23','DIFFABS']:
        return cal_UB_p23
    elif name =='ID03':
        return cal_UB_id03

def normalize_img_intensity(img, q_grid,mask_img, mask_profile, cen, offset, direction = 'horizontal'):
    if direction == 'horizontal':
        cut_mask = np.sum(mask_img[cen-offset:cen+offset+1,:], axis=0)
        cut_img = np.sum(img[cen-offset:cen+offset+1,:],axis=0)
        cut_img = cut_img/cut_mask
        cut_img = cut_img[mask_profile]
        cut_q = q_grid[0,mask_profile]
        #now remove nan values
        cut_img_nan = cut_img == np.nan
        return cut_img[~cut_img_nan],cut_q[~cut_img_nan]
    elif direction == 'vertical':
        cut_mask = np.sum(mask_img[:,cen-offset:cen+offset+1], axis=1)
        cut_img = np.sum(img[:,cen-offset:cen+offset+1],axis=1)
        cut_img = cut_img/cut_mask
        cut_img = cut_img[mask_profile]
        cut_q = q_grid[mask_profile,0]
        #now remove nan values
        cut_img_nan = cut_img == np.nan
        return cut_img[~cut_img_nan],cut_q[~cut_img_nan]

def update_bounds(original_bounds, guess_values,offset=0.1, which=0):
    original_bounds_left, original_bounds_right = original_bounds
    guess_value = guess_values[which]
    if guess_value > original_bounds_right[which]:
        original_bounds_right[which] = guess_value + offset
    elif guess_value < original_bounds_left[which]:
        original_bounds_left[which] = guess_value - offset
    return [original_bounds_left,original-bounds_right]

def debug_output(local_lib,var_list,debug):
    if debug:
        for var in var_list:
            print(var,local_lib[var])
    return None

class nexus_image_loader_diffabs(object):
    def __init__(self,nexus_path='/home/qiu/data/beamtime/P23_11_18_I20180114/raw/FirstTest_00666/lmbd',frame_prefix='FirstTest'):
        self.nexus_path=nexus_path
        self.frame_prefix=frame_prefix
        # self.get_frame_number()
        self.current_scan_no=None
        self.current_scan_img_no=None
        self.current_energy = None

    def get_frame_number(self, frame_prefix,scan_no):
        img_name='{}_{}.nxs'.format(frame_prefix,current_scan_no)
        img_path=os.path.join(self.nexus_path,img_name)
        # data=nxload(img_path)
        data=h5py.File(img_path,'r')
        access_path="scan_{}/scan_data/data_57".format(scan_number)
        total_img_number = data[access_path].shape[0]
        self.total_frame_number = total_img_number
        return self.total_frame_number

    def load_frame(self,scan_number,frame_number,offset={},flip=True):
        img_name='{}_{}.nxs'.format(self.frame_prefix,scan_number)
        self.current_scan_no, self.current_scan_img_no = scan_number, frame_number
        img_path=os.path.join(self.nexus_path,img_name)
        # data=nxload(img_path)
        data=h5py.File(img_path,'r')
        access_path="scan_{}/scan_data/data_57".format(scan_number)
        img=np.array(data[access_path][frame_number])
        self.total_frame_number = data[access_path].shape[0]
        self.extract_motor_angles(data, scan_number, frame_number,offset)
        self.update_energy(data, scan_number, frame_number)
        if flip:
            return np.flip(img.T,1)
        else:
            return img

    def extract_motor_angles(self, scan, scan_number,frame_number, offset ={'delta':0}):
        motors={}
        common_head_header = "scan_{}/DIFFABS/D13-1-CX1__EX__DIF.1-{}__#1/raw_value"
        common_head_counter ="scan_{}/scan_data/{}"
        motor_names = ['phi', 'chi', 'delta', 'gamma', 'mu', 'omega_t']
        motor_map_lib = {'phi':["header","PHI_E"], "chi":["header","CHI_E"], \
                         'delta': ["counter", "actuator_1_4"], "gamma":["header","GAMMA"],\
                         "mu":['header','MU'],'omega_t':['counter','actuator_1_3']}
        for key in motor_map_lib:
            pot, map_key = motor_map_lib[key]
            _offset = 0
            if key in offset.keys():
                _offset = offset[key]
            if pot == "header":
                # print(key,common_head_header.format(scan_number,map_key))
                motors[key] =np.array(scan[common_head_header.format(scan_number,map_key)])[0]+_offset
            elif pot == "counter":
                # print(key, map_key,common_head_counter.format(scan_number,map_key))
                motors[key] = np.array(scan[common_head_counter.format(scan_number,map_key)])[frame_number]+_offset
        self.motor_angles = motors

    def update_energy(self, scan, scan_number, frame_number):
        e_counter ="scan_{}/scan_data/actuator_1_1".format(scan_number)
        energy = np.array(scan[e_counter])[frame_number]
        self.current_energy = energy
        return energy

    def load_frame_from_path(self,img_path,frame_number = 0,flip=True):
        try:
            #if one frame one nxs file
            data=nxload(img_path)
            img=np.array(data.entry.instrument.detector.data.nxdata[0])
        except:
            #if all frames in one nxs file
            data=nxload(img_path)
            img=np.array(data.entry.instrument.detector.data[frame_number])
        if flip:
            return np.flip(img.T,1)
        else:
            return img

    def show_frame(self,scan_number,frame_number,one_frame_in_one_nxs=True,flip=True):
        img=self.load_frame(scan_number,frame_number,one_frame_in_one_nxs,flip)
        fig,ax=pyplot.subplots()
        pyplot.imshow(img,cmap='jet')
        if flip:
            pyplot.colorbar(extend='both',orientation='vertical')
        else:
            pyplot.colorbar(extend='both',orientation='horizontal')
        pyplot.clim(0,205)
        # pyplot.show()
        return img

    def find_dead_pix(self,scan_number=666,img_end=100):
        dead_pix_container=self.load_frame(scan_number,0)==self.load_frame(scan_number,1)
        dead_pix_container=np.where(dead_pix_container==True)
        dead_pix_container=zip(tuple(dead_pix_container[0]),tuple(dead_pix_container[1]))
        img0= self.load_frame(scan_number,0)
        print(len(dead_pix_container))
        for i in range(2,img_end):
            print('Processing img_',i)
            img = self.load_frame(scan_number,i)
            temp= img != img0
            temp= np.where(temp==True)
            temp= zip(tuple(temp[0]),tuple(temp[1]))
            for each in temp:
                if each in dead_pix_container:
                    dead_pix_container.remove(each)
        return dead_pix_container

class nexus_image_loader(object):
    def __init__(self,fio_path='/home/qiu/data/beamtime/P23_11_18_I20180114/raw/startup/FirstTest_00666.fio',nexus_path='/home/qiu/data/beamtime/P23_11_18_I20180114/raw/FirstTest_00666/lmbd',frame_prefix='FirstTest', scan_number = 126):
        self.fio_path=fio_path
        self.nexus_path=nexus_path
        self.frame_prefix=frame_prefix
        #load nexus data only once here
        img_name='{}_{:0>5}.nxs'.format(frame_prefix,scan_number)
        img_path=os.path.join(self.nexus_path,img_name)
        self.nexus_data = nxload(img_path)
        self.get_frame_number()

    def get_frame_number(self):
        #total_img_number = len(os.listdir(self.nexus_path))
        #if total_img_number ==1:
            #img_name = os.listdir(self.nexus_path)[0]
            #img_path=os.path.join(self.nexus_path,img_name)
            #data=nxload(img_path)
            #try:
            #    total_img_number = len(np.array(self.nexus_data.entry.instrument.detector.data))
            #except:
            #    total_img_number = len(np.array(self.nexus_data.scan.data.atten))
            #total_img_number = len(np.array(self.nexus_data.scan.data.atten))
            #print(total_img_number)
        total_img_number = len(np.array(self.nexus_data.scan.data.eh_c01))
        self.total_frame_number = total_img_number
        return total_img_number

    def load_frame(self,scan_number,frame_number,flip=True):
        try:
            #if one frame one nxs file
            img_name='{}_{:0>5}_{:0>5}.nxs'.format(self.frame_prefix,scan_number,frame_number)
            img_path=os.path.join(self.nexus_path,img_name)
            data=nxload(img_path)
            img=np.array(data.entry.instrument.detector.data.nxdata[0])
        except:
            #if all frames in one nxs file
            img_name='{}_{:0>5}.nxs'.format(self.frame_prefix,scan_number)
            img_path=os.path.join(self.nexus_path,img_name)
            data=nxload(img_path)
            img=np.array(data.entry.instrument.detector.data[frame_number])
        if flip:
            return np.flip(img.T,1)
        else:
            return img

    def load_frame_new(self,frame_number,flip=True, clip_boundary = {'ver':[0,10000],'hor':[0,10000]}):
        #if one frame one nxs file
        #img_name='{}_{:0>5}.nxs'.format(self.frame_prefix,scan_number)
        #img_path=os.path.join(self.nexus_path,img_name)
        #data=nxload(img_path)
        #img=np.array(data.entry.instrument.detector.data.nxdata[0])
        img=np.array(self.nexus_data.scan.data.lmbd)[frame_number]
        if flip:
            img = np.flip(img.T,1)
        img = img[clip_boundary['ver'][0]:clip_boundary['ver'][1],clip_boundary['hor'][0]:clip_boundary['hor'][1]]
        return img
            
    def extract_beam_mon_ct(self,mon_path = 'scan/data/eh_c01'):
        return np.array(self.nexus_data['scan/data/eh_c01'])


    def extract_motor_angles(self, frame_number, constant_motors ={'omega_t':0.5, 'phi':0, 'chi':90}):
        #img_name='{}_{:0>5}.nxs'.format(self.frame_prefix,scan_number)
        #img_path=os.path.join(self.nexus_path,img_name)
        #data=nxload(img_path)
        motors={}
        motor_names = ['phi', 'chi', 'delta', 'gamma', 'mu', 'omega_t']
        for motor in constant_motors:
            motors[motor] = constant_motors[motor]
        for motor in motor_names:
            if motor not in motors.keys():
                fetch_path = 'scan/data/{}'.format(motor)
                motors[motor] = np.array(self.nexus_data[fetch_path])[frame_number] 
        motors['mon'] = np.array(self.nexus_data['scan/data/eh_c01'])[frame_number]
        try:
            motors['transm']=1./np.array(self.nexus_data['scan/data/atten'])[frame_number]
        except:
            motors['transm']=np.array(self.nexus_data['scan/data/lmbd_countsroi1'])[frame_number]/np.array(self.nexus_data['scan/data/lmbd_countsroi1_atten'])[frame_number]

        self.motor_angles = motors
        return motors
        
    def update_motor_angles_in_data(self,data):
        for motor in self.motor_angles:
            data[motor].append(self.motor_angles[motor])
        return data

    def extract_pot_current(self, frame_number):
        pot = np.array(self.nexus_data['scan/data/voltage2'])[frame_number]
        cur = np.array(self.nexus_data['scan/data/voltage1'])[frame_number]
        return pot, cur

    def extract_pot_profile(self):
        return np.array(self.nexus_data['scan/data/voltage2'])
        
    def extract_HKL(self, frame_number):
        H = np.array(self.nexus_data['scan/data/diffractometer_h'])[frame_number]
        K = np.array(self.nexus_data['scan/data/diffractometer_k'])[frame_number]
        L = np.array(self.nexus_data['scan/data/diffractometer_l'])[frame_number]
        #cur = np.array(self.nexus_data['scan/data/voltage1'])[frame_number]
        return H, K, L
        
    def load_frame_from_path(self,img_path,frame_number = 0,flip=True):
        try:
            #if one frame one nxs file
            data=nxload(img_path)
            img=np.array(data.entry.instrument.detector.data.nxdata[0])
        except:
            #if all frames in one nxs file
            data=nxload(img_path)
            img=np.array(data.entry.instrument.detector.data[frame_number])
        if flip:
            return np.flip(img.T,1)
        else:
            return img

    def show_frame(self,scan_number,frame_number,one_frame_in_one_nxs=True,flip=True):
        img=self.load_frame(scan_number,frame_number,one_frame_in_one_nxs,flip)
        fig,ax=pyplot.subplots()
        pyplot.imshow(img,cmap='jet')
        if flip:
            pyplot.colorbar(extend='both',orientation='vertical')
        else:
            pyplot.colorbar(extend='both',orientation='horizontal')
        pyplot.clim(0,205)
        # pyplot.show()
        return img

    def find_dead_pix(self,scan_number=666,img_end=100):
        dead_pix_container=self.load_frame(scan_number,0)==self.load_frame(scan_number,1)
        dead_pix_container=np.where(dead_pix_container==True)
        dead_pix_container=zip(tuple(dead_pix_container[0]),tuple(dead_pix_container[1]))
        img0= self.load_frame(scan_number,0)
        print(len(dead_pix_container))
        for i in range(2,img_end):
            print('Processing img_',i)
            img = self.load_frame(scan_number,i)
            temp= img != img0
            temp= np.where(temp==True)
            temp= zip(tuple(temp[0]),tuple(temp[1]))
            for each in temp:
                if each in dead_pix_container:
                    dead_pix_container.remove(each)
        return dead_pix_container


class DetImage:
    def __init__(self, img, motors, counters, header=None):
        self.img = img
        self.motors = motors
        self.counters = counters
        self.header = header

#extract info from spec file
class pyspec(object):
    def __init__(self, spec_name = "/home/qiu/data/CH5314/ch5314_sixcvertical.spec"):
        self.spec_file = spec_name
        self.data= {}
        self.headers = {}
        self.comments = {}
        self.data_label = {}
        self.extract_info()

    def extract_info(self, ):
        with open(self.spec_file) as f_spec:
            spec_lines=f_spec.readlines()
            current_scan = " "
            for i in range(len(spec_lines)):
                if spec_lines[i].startswith("#S"):
                    current_scan = int(spec_lines[i].rsplit()[1])
                    self.data[current_scan] = []
                    self.headers[current_scan] = [spec_lines[i]]
                    self.comments[current_scan] = []
                elif spec_lines[i].startswith("#C"):
                    if current_scan!=" ":
                        self.comments[current_scan].append(spec_lines[i])
                elif spec_lines[i].startswith("#L"):
                    if current_scan!=" ":
                        self.data_label[current_scan] = spec_lines[i].rstrip().rsplit()[1:]
                elif spec_lines[i].startswith("#"):
                    if current_scan!=" ":
                        self.headers[current_scan].append(spec_lines[i])
                elif spec_lines[i] == "\n":
                    pass
                else:
                    if current_scan!=" ":
                        try:
                            self.data[current_scan].append([float(k) for k in spec_lines[i].rstrip().rsplit()])
                        except:
                            print('Check scan_{} in the spec file'.format(current_scan))
'''
    Loads the images from a spec scan.
'''
class edf_image_loader_old:
    def __init__(self, spec_filename, image_foldername,is_zap_scan):
        self.spec_file = spec.SpecDataFile(spec_filename)
        self.scan_selector = specfilewrapper.Specfile(spec_filename)
        self.spec_filename = spec_filename
        self.image_foldername = image_foldername
        self.total_frame_number = None
        self.is_zap_scan = is_zap_scan

    def get_frame_number(self, scan_no):
        return self.spec_file[scan_no].data.shape[0]

    def get_trans_mon_factor(self,scan,frame_no):
        mon, trans = None, None
        if self.is_zap_scan:
            mon = scan.datacol('zap_mon')[frame_no]
            trans = scan.datacol('zap_transm')[frame_no]
        else:
            counter_prefix = ''
            if 'ccoscan' in scan.header('S')[0]:
                counter_prefix = 'zap_'
            mon = scan.datacol('%smon'%(counter_prefix))[frame_no]
            trans = scan.datacol('%stransm'%(counter_prefix))[frame_no]
        self.transm = trans
        self.mon = mon
        return mon*trans

    def load_frame(self, scan_no, frame_no, gz_compressed=True, normalize=False, monitor_name=None, monitor_names=None, remove_rows=None, remove_cols=None):
        self.total_frame_number = self.get_frame_number(scan_no)
        mon_trans_factor = self.get_trans_mon_factor(self.scan_selector.select('{}.1'.format(scan_no)),frame_no)
        header = dict()
        for line in self.spec_file[scan_no].header.split('\n'):
            header[line.split(' ')[0]] = line[len(line.split(' ')[0]):]
        comments = dict()
        for line in self.spec_file[scan_no].comments.split('\n'):
            comments[line.split(':')[0].strip()] = line[len(line.split(':')[0])+1:].strip()
        frames = self.total_frame_number

        if('ccoscan' in header['#S'] or 'zapline' in header['#S']):
            img_folder = comments['#C DIRECTORY'].split('/')[-1]+'/'
            if(img_folder == '/'):
                img_folder = comments['#C DIRECTORY'].split('/')[-2]+'/'
            zap_scan_no = int(comments['#C ZAP SCAN NUMBER'])
            radix = comments['#C RADIX']
            filename = radix + '_mpx-x4_%s_0000_0000.edf'%(str(zap_scan_no).zfill(4))
            multiframe_edf_frame_no = frame_no
        else:
            first_frame = int(header['#UCCD'].split('#r')[-1].split('.')[0])
            img_folder = header['#UCCD'].replace('//','/').split('/')[-2]+'/' #in the beginning the path of MA3886 had '//' for some reason
            filename_template = header['#UCCD'].split('/')[-1].split('#r')[0] + '#r.' + header['#UCCD'].split('/')[-1].split('.')[-1]
            filename = filename_template.replace('#n', str(scan_no).zfill(3)).replace('#p', str(frame_no).zfill(3)).replace('#r', str(first_frame+frame_no).zfill(3))
            multiframe_edf_frame_no = 0

        # automatically detect if frame is compressed or not
        if(not os.path.exists(os.path.join(self.image_foldername, img_folder, filename))):
            filename = filename.replace('.edf', '.edf.gz')

        #if(gz_compressed):
        #    filename = filename.replace('.edf', '.edf.gz')

        if(frame_no >= frames):
            raise IndexError("Frame number does not exist.")

        img_filename = os.path.join(self.image_foldername, img_folder, filename)
        edf = EdfFile.EdfFile(img_filename, 'r')


        edf_header = edf.GetHeader(multiframe_edf_frame_no)
        motors = dict()
        motor_mne = edf_header['motor_mne'].split()
        motor_pos = edf_header['motor_pos'].split()
        for i in range(len(motor_mne)):
            motors[motor_mne[i]] = float(motor_pos[i])
        counters = dict()
        if(not ('ccoscan' in header['#S'] or 'zapline' in header['#S'])):
            counter_mne = edf_header['counter_mne'].split()
            counter_pos = edf_header['counter_pos'].split()
            for i in range(len(counter_mne)):
                counters[counter_mne[i]] = float(counter_pos[i])

        the_img = np.array(edf.GetData(multiframe_edf_frame_no), dtype='float')
        if(remove_rows != None):
            the_img = np.delete(the_img, remove_rows, axis=0)
        if(remove_cols != None):
            the_img = np.delete(the_img, remove_cols, axis=1)

        if(normalize):
            if(monitor_name):
                mon_count = float(getattr(self.spec_file[scan_no], monitor_name)[frame_no])
                the_img /= mon_count
            if(monitor_names):
                for mon_name in monitor_names:
                    mon_count = float(getattr(self.spec_file[scan_no], mon_name)[frame_no])
                    the_img /= mon_count

        return DetImage(the_img/mon_trans_factor, motors, counters, header)

    def load_all_frames(self, scan_no, gz_compressed=True, normalize=False, monitor_name=None, remove_rows=None, remove_cols=None):
        frame_no = self.get_no_frames(scan_no)
        frames = np.zeros(frame_no)
        for i in range(frame_no):
            frames[i] = self.load_frame(scan_no, i, gz_compressed, normalize, monitor_name, remove_rows, remove_cols)
        return frame

class tiff_image_loader:
    def __init__(self, spec_filename, image_foldername):
        self.scan_selector = specfilewrapper.Specfile(spec_filename)
        self.spec_filename = spec_filename
        self.image_foldername = image_foldername
        self.total_frame_number = None
        # print(spec_filename)

    def get_frame_number(self, scan_no):
        return self.scan_selector.select('{}.1'.format(scan_no)).data().shape[1]

    def get_trans_mon_factor(self,scan,frame_no):
        mon, trans = None, None

        mon = scan.datacol('io')[frame_no]
        trans = scan.datacol('transm')[frame_no]
        # self.transm = trans
        self.transm = 1
        self.mon = mon
        return mon*trans

    def load_frame(self, scan_no, frame_no):
        self.total_frame_number = self.get_frame_number(scan_no)
        mon_trans_factor = self.get_trans_mon_factor(self.scan_selector.select('{}.1'.format(scan_no)),frame_no)
        scan_folder = 'S{:0>3}'.format(scan_no)
        img_name = '{}_{}_{:0>5}.tif'.format(os.path.basename(self.image_foldername),\
                                         scan_folder, frame_no)
        img_path = os.path.join(self.image_foldername, scan_folder, img_name)
        img=misc.imread(img_path)
        return img/mon_trans_factor

class edf_image_loader:
    def __init__(self, spec_filename, image_foldername,is_zap_scan):
        self.spec_file = pyspec(spec_filename)
        self.scan_selector = specfilewrapper.Specfile(spec_filename)
        self.spec_filename = spec_filename
        self.image_foldername = image_foldername
        self.total_frame_number = None
        self.is_zap_scan = is_zap_scan

    def get_frame_number(self, scan_no):
        return np.array(self.spec_file.data[scan_no]).shape[0]

    def get_trans_mon_factor(self,scan,frame_no):
        mon, trans = None, None
        if self.is_zap_scan:
            mon = scan.datacol('zap_mon')[frame_no]
            trans = scan.datacol('zap_transm')[frame_no]
        else:
            counter_prefix = ''
            if 'ccoscan' in scan.header('S')[0]:
                counter_prefix = 'zap_'
            try:
                mon = scan.datacol('%smon'%(counter_prefix))[frame_no]
            except:
                # mon = scan.datacol('%sattn'%(counter_prefix))[frame_no]
                mon = 1
            trans = scan.datacol('%stransm'%(counter_prefix))[frame_no]
        self.transm = trans
        self.mon = mon
        return mon*trans

    def load_frame(self, scan_no, frame_no, gz_compressed=False, normalize=False, monitor_name=None, monitor_names=None, remove_rows=None, remove_cols=None):
        self.total_frame_number = self.get_frame_number(scan_no)
        mon_trans_factor = self.get_trans_mon_factor(self.scan_selector.select('{}.1'.format(scan_no)),frame_no)
        header = dict()
        for line in self.spec_file.headers[scan_no]:
            # print(self.spec_file.headers[scan_no])
            # print(line.split()[0])
            header[line.split()[0]] = line[len(line.split()[0]):]
        comments = dict()
        for line in self.spec_file.comments[scan_no]:
            comments[line.split(':')[0].strip()] = line[len(line.split(':')[0])+1:].strip()
        frames = self.total_frame_number

        if('ccoscan' in header['#S'] or 'zapline' in header['#S']):
            img_folder = comments['#C DIRECTORY'].split('/')[-1]
            if(img_folder in ['','/']):
                img_folder = comments['#C DIRECTORY'].split('/')[-2]
            zap_scan_no = int(comments['#C ZAP SCAN NUMBER'])
            radix = comments['#C RADIX']
            filename = radix + '_mpx-x4_%s_0000_0000.edf'%(str(zap_scan_no).zfill(4))
            multiframe_edf_frame_no = frame_no
        else:
            first_frame = int(header['#UCCD'].split('#r')[-1].split('.')[0])
            img_folder = header['#UCCD'].replace('//','/').split('/')[-2]+'/' #in the beginning the path of MA3886 had '//' for some reason
            filename_template = header['#UCCD'].split('/')[-1].split('#r')[0] + '#r.' + header['#UCCD'].split('/')[-1].split('.')[-1]
            filename = filename_template.replace('#n', str(scan_no).zfill(3)).replace('#p', str(frame_no).zfill(3)).replace('#r', str(first_frame+frame_no).zfill(3))
            multiframe_edf_frame_no = 0

        # automatically detect if frame is compressed or not
        # if(not os.path.exists(os.path.join(self.image_foldername, img_folder, filename))):
            # filename = filename.replace('.edf', '.edf.gz')

        #if(gz_compressed):
        #    filename = filename.replace('.edf', '.edf.gz')

        if(frame_no >= frames):
            raise IndexError("Frame number does not exist.")

        img_filename = os.path.join(self.image_foldername, img_folder, filename.rstrip())
        edf = EdfFile.EdfFile(img_filename, 'r')
        # print(self.image_foldername,img_folder,filename)
        # edf = EdfFile.EdfFile("/Users/cqiu/app/ma4171/ma4171_img/ma4171_mpx02/ma4171_mpx_247_000_1977.edf", 'r')


        edf_header = edf.GetHeader(multiframe_edf_frame_no)
        motors = dict()
        motor_mne = edf_header['motor_mne'].split()
        motor_pos = edf_header['motor_pos'].split()
        for i in range(len(motor_mne)):
            motors[motor_mne[i]] = float(motor_pos[i])
        counters = dict()
        if(not ('ccoscan' in header['#S'] or 'zapline' in header['#S'])):
            counter_mne = edf_header['counter_mne'].split()
            counter_pos = edf_header['counter_pos'].split()
            for i in range(len(counter_mne)):
                counters[counter_mne[i]] = float(counter_pos[i])

        the_img = np.array(edf.GetData(multiframe_edf_frame_no), dtype='float')
        if(remove_rows != None):
            the_img = np.delete(the_img, remove_rows, axis=0)
        if(remove_cols != None):
            the_img = np.delete(the_img, remove_cols, axis=1)

        if(normalize):
            if(monitor_name):
                mon_count = self.spec_file.data[scan_no][frame_no][self.spec_file.data_label.index(monitor_name)]
                the_img /= mon_count
            if(monitor_names):
                for mon_name in monitor_names:
                    mon_count = self.spec_file.data[scan_no][frame_no][self.spec_file.data_label.index(mon_name)]
                    the_img /= mon_count

        return DetImage(the_img/mon_trans_factor, motors, counters, header)

    def load_all_frames(self, scan_no, gz_compressed=True, normalize=False, monitor_name=None, remove_rows=None, remove_cols=None):
        frame_no = self.get_no_frames(scan_no)
        frames = np.zeros(frame_no)
        for i in range(frame_no):
            frames[i] = self.load_frame(scan_no, i, gz_compressed, normalize, monitor_name, remove_rows, remove_cols)
        return frame

