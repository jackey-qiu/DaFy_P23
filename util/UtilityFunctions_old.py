import numpy as np
from Reciprocal_space_tools.HKLVlieg import Crystal, printPos, UBCalculator, VliegAngles, printPos_prim, vliegDiffracAngles
from nexusformat.nexus import *
import fnmatch, os
import re,os
from PyMca5.PyMcaIO import specfilewrapper 
from pyspec import spec
from PyMca5.PyMcaIO import EdfFile

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
    if name == 'P23':
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

class nexus_image_loader(object):
    def __init__(self,fio_path='/home/qiu/data/beamtime/P23_11_18_I20180114/raw/startup/FirstTest_00666.fio',nexus_path='/home/qiu/data/beamtime/P23_11_18_I20180114/raw/FirstTest_00666/lmbd',frame_prefix='FirstTest'):
        self.fio_path=fio_path
        self.nexus_path=nexus_path
        self.frame_prefix=frame_prefix
        self.get_frame_number()

    def get_frame_number(self):
        total_img_number = len(os.listdir(self.nexus_path))
        if total_img_number ==1:
            img_name = os.listdir(self.nexus_path)[0]
            img_path=os.path.join(self.nexus_path,img_name)
            data=nxload(img_path)
            total_img_number = len(np.array(data.entry.instrument.detector.data))
        self.total_frame_number = total_img_number
        return self.total_frame_number

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

'''
    Loads the images from a spec scan.
'''
class edf_image_loader:
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
        for i in xrange(len(motor_mne)):
            motors[motor_mne[i]] = float(motor_pos[i])
        counters = dict()
        if(not ('ccoscan' in header['#S'] or 'zapline' in header['#S'])):
            counter_mne = edf_header['counter_mne'].split()
            counter_pos = edf_header['counter_pos'].split()
            for i in xrange(len(counter_mne)):
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
        for i in xrange(frame_no):
            frames[i] = self.load_frame(scan_no, i, gz_compressed, normalize, monitor_name, remove_rows, remove_cols)
        return frames
