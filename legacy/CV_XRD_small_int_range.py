import sys
from numpy import dtype
sys.path.append('/home/qiu/apps/2017_MA3589_data_analysis/')
sys.path.append('/home/qiu/files/surf_diff/XRD_tools/')
sys.path.append('..')
# from MA3589 import *
from MA3886 import *
#import XRD_tools.id03_tools as id03
import numpy as np
from scipy.interpolate import griddata
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter

import matplotlib.pyplot as plt
from pyspec import spec

class ScanInfo():
    def __init__(self, scan_id, scan_no, beamtime_config, specfile, DI, int_lim_ip, int_lim_oop, fit_lim_ip=None, fit_lim_oop=None, lattice=MA3886.Co3O4_lat, HKL=[4,0,4]):
        self.scan_id = scan_id
        self.beamtime_config = beamtime_config
        self.specfile = specfile
        self.DI = DI
        self.scan_no = scan_no
        self.int_lim_ip = int_lim_ip
        self.int_lim_oop = int_lim_oop
        self.fit_lim_ip = fit_lim_ip
        self.fit_lim_oop = fit_lim_oop
        if self.fit_lim_ip == None:
            self.fit_lim_ip = [0,515]
        if self.fit_lim_oop == None:
            self.fit_lim_oop = [0,515]
        self.lattice=lattice
        self.HKL = HKL

class ScanInfoContainer(dict):
    def __init__(self):
        return
    def add(self, scan_id, *args, **kwargs):
        self[scan_id] = ScanInfo(scan_id, *args, **kwargs)

class fit_data():
    def __init__(self):
        self.potential = []
        self.current_density = []
        self.Time = []
        self.pcov_ip = []
        self.pcov_oop = []
        self.cen_ip = []
        self.FWHM_ip = []
        self.amp_ip = []
        self.lfrac_ip = []
        self.bg_slope_ip = []
        self.bg_offset_ip = []
        self.cen_oop = []
        self.FWHM_oop = []
        self.amp_oop = []
        self.lfrac_oop = []
        self.bg_slope_oop = []
        self.bg_offset_oop = []
    def print_at(self, index):
        for attr, value in self.__dict__.iteritems():
            print attr, '=', value[index]




def gauss(x, x0, sig, amp):
    return amp*np.exp(-(x-x0)**2/2./sig**2)

def lor(x, x0, FWHM, amp):
    return amp*FWHM/((x-x0)**2+FWHM**2/4)

def pvoigt2(x, x0, FWHM, amp, lorfact):
    w = FWHM/2.
    return amp*(lorfact/(1+((x-x0)/w)**2)+(1.-lorfact)*np.exp(-np.log(2)*((x-x0)/w)**2))

def pvoigt(x, x0, FWHM, area, lfrac):
    return area / FWHM / ( lfrac*np.pi/2 + (1-lfrac)*np.sqrt(np.pi/4/np.log(2)) ) * ( lfrac / (1 + 4*((x-x0)/FWHM)**2) + (1-lfrac)*np.exp(-4*np.log(2)*((x-x0)/FWHM)**2) )

def model2(x, x0, FWHM, amp, bg_slope, bg_offset):
    return lor(x, x0, FWHM, amp) + x*bg_slope*0 + bg_offset

def model3(x, x0, FWHM, amp, bg_slope, bg_offset):
    sig = FWHM/2.35482
    return gauss(x, x0, sig, amp) + x*bg_slope*0 + bg_offset

def model(x, x0, FWHM, area, lfrac, bg_slope, bg_offset):
    return pvoigt(x, x0, FWHM, area, lfrac) + x*bg_slope + bg_offset


scan_info = ScanInfoContainer()
################################################################################################################################
# MA3589
################################################################################################################################

#scan_info.add(scan_no=161, int_lim_ip=[50, 160], int_lim_oop=[150, 250], fit_lim_ip=[50,160], fit_lim_oop=[180, 450])
# scan_info.add('MA3589_161', scan_no=161, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 110], int_lim_oop=[190, 220], fit_lim_ip=[50,160], fit_lim_oop=None)


######################################################################################
# Co3O4 sample
######################################################################################

# Co3O4(3-33)
# scan_info.add('MA3589_388', scan_no=388, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_396', scan_no=396, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_399', scan_no=399, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_409', scan_no=409, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_414', scan_no=414, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_417', scan_no=417, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_426', scan_no=426, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])

# Co3O4(044)
# scan_info.add('MA3589_392', scan_no=392, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_397', scan_no=397, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_410', scan_no=410, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_415', scan_no=415, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_418', scan_no=418, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_427', scan_no=427, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])

# exclude Au from fitting range
# scan_info.add('MA3589_392', scan_no=392, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_397', scan_no=397, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_410', scan_no=410, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_415', scan_no=415, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_418', scan_no=418, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_427', scan_no=427, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])


# CoOOH(015) -> nothing visible
#393 atten too high
# scan_info.add('MA3589_394', scan_no=394, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_400', scan_no=400, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])

# CoO2(012) -> nothing visible
# scan_info.add('MA3589_395', scan_no=395, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_401', scan_no=401, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])

# Co(OH)2(-112) -> at sufficiently negative potential a rod appears
# scan_info.add('MA3589_416', scan_no=416, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_419', scan_no=419, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])
# scan_info.add('MA3589_420', scan_no=420, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])

# scan_info.add('MA3589_424', scan_no=424, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[20, 250], int_lim_oop=[50, 350], fit_lim_ip=[[10, 250]], fit_lim_oop=[10, 500])

######################################################################################
# CoOOH sample
######################################################################################

# CoOOH(012)
# scan_info.add('MA3589_162', scan_no=162, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_163', scan_no=163, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_183', scan_no=183, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_185', scan_no=185, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_263', scan_no=263, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_276', scan_no=276, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_277', scan_no=277, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_278', scan_no=278, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_309', scan_no=309, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])
# scan_info.add('MA3589_310', scan_no=310, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[150, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[170, 500])

# CoOOH(017)
# scan_info.add('MA3589_264', scan_no=264, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[100, 150], int_lim_oop=[205, 250], fit_lim_ip=[[10, 500]], fit_lim_oop=[10, 500])

# CoOOH(003)
# scan_info.add('MA3589_258', scan_no=258, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[95, 150], int_lim_oop=[205, 250], fit_lim_ip=[[10, 200]], fit_lim_oop=[200, 500])

# Co3O4(113)
# intensity too low for small integration range
# scan_info.add('MA3589_255', scan_no=255, beamtime_config=MA3589, specfile=MA3589.spec_file, DI=MA3589.DI, int_lim_ip=[70, 150], int_lim_oop=[150, 250], fit_lim_ip=[[50, 170]], fit_lim_oop=[140, 250])

################################################################################################################################
# MA3886 - Co3O4
################################################################################################################################
scan_info.add('MA3886_443', scan_no=443, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3886_448', scan_no=448, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3886_449', scan_no=449, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3886_450', scan_no=450, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])

# scan_info.add('MA3886_451', scan_no=451, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 90],[110,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3886_197', scan_no=197, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 85],[115,250]], fit_lim_oop=[10, 500])
# scan_info.add('MA3886_254', scan_no=254, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[20, 85],[115,250]], fit_lim_oop=[10, 500])

# pH dependence
# scan_info.add('MA3886_273', scan_no=273, beamtime_config=MA3886, specfile=MA3886.spec_file, DI=MA3886.DI, int_lim_ip=[110, 250], int_lim_oop=[50, 350], fit_lim_ip=[[0, 85],[115,500]], fit_lim_oop=[10, 500])




################################################################################################################################
# CH5314 - CoOOH
################################################################################################################################
scan_info.add('CH5314_056', scan_no=56, beamtime_config=MA3886, specfile=MA3886.spec_file2, DI=MA3886.DI2, int_lim_ip=[50, 150], int_lim_oop=[50, 350], fit_lim_ip=[[20,250]], fit_lim_oop=[50, 400], lattice=MA3886.CoOOH_lat, HKL=[0,1,7])
# scan_info.add('CH5314_060', scan_no=60, beamtime_config=MA3886, specfile=MA3886.spec_file2, DI=MA3886.DI2, int_lim_ip=[50, 150], int_lim_oop=[50, 350], fit_lim_ip=[[20,250]], fit_lim_oop=[50, 400], lattice=MA3589.CoOOH_lat, HKL=[0,1,7])



scan_id =  'MA3886_443'#'MA3886_450'#CH5314_060' #'MA3886_443'#
is_zap_scan = 0
save_all_frames = 0

f = scan_info[scan_id].specfile
DI = scan_info[scan_id].DI
scan_no = scan_info[scan_id].scan_no


if('zap' in f[scan_info[scan_id].scan_no].header):
    is_zap_scan = True

data = fit_data()

plt.ion()
fig = plt.figure()
fig_cut_ip = plt.figure()
fig_cut_oop = plt.figure()
#fig_d_oop = plt.figure()



# since the detector does not move during the scan the grid_intensity always puts the same pixels in the same bin
# this can be calculated once which should speed up things alot.
if(is_zap_scan):
    number_of_frames = len(f[scan_no].mottime)
    potential = f[scan_no].zap_Eana
    current = f[scan_no].zap_Iana

else:
    number_of_frames = len(f[scan_no].Time)
    potential = f[scan_no].Eana
    current = f[scan_no].Iana



Intensity, q = DI.get_q(scan_no, 0, gz_compressed=True)
q = np.array(q)
if('MA3589' in scan_id):
    q = MA3589.q_correction(q, 2 if scan_no > 350 else 1)
q_ip = np.sqrt(q[0]**2 + q[1]**2)
q_oop = q[2]


if(scan_no == 258):
    # TODO:
    # this is a bit tricky to implement in q_ip and q_oop
    # if this is to be used a bit more thought has to be put into this
    for row in xrange(q_ip.shape[1]):
        last_qip = 1000
        for col in xrange(q_ip.shape[0]):
            if(q_ip[row,col] < last_qip):
                last_qip = q_ip[row,col]
                q_ip[row,col] *= -1
            else:
                break


# place image on regular grid
size = Intensity.shape
grid_q_oop, grid_q_ip = np.mgrid[np.max(q_oop):np.min(q_oop):(1.j*size[1]), np.min(q_ip):np.max(q_ip):(1.j*size[0])]
grid_indices = griddata((q_ip.ravel(), q_oop.ravel()), np.arange(size[0]*size[1]).ravel(), (grid_q_ip, grid_q_oop), method='nearest')
# print grid_indices
# quit()

for frame_no in [0]:#xrange(number_of_frames):
    print 'Frame number :', frame_no, '/', number_of_frames
#for frame_no in xrange(1,2):#(70,95):

    Intensity, q = DI.get_q(scan_no, frame_no, gz_compressed=False)

    flat_int = Intensity.ravel()
    grid_intensity = flat_int[grid_indices]
    grid_intensity = grid_intensity.reshape(size)

    mask = np.ones(grid_intensity.shape)
    mask_error_value = 0
    if(scan_id in ['MA3886_443', 'MA3886_448', 'MA3886_449', 'MA3886_450','MA3886_451']):
        vmax_2d_img = 0.02
        vmin_2d_img = 0
        # mask dead pixels
        mask[152, 146] = mask_error_value
        mask[176, 131:133] = mask_error_value
        mask[172, 143:145] = mask_error_value
        mask[172, 147:150] = mask_error_value
        mask[173, 145:147] = mask_error_value
        mask[173, 148] = mask_error_value
        mask[191, 145] = mask_error_value
        mask[192, 145:147] = mask_error_value
        mask[193, 144] = mask_error_value
        mask[193, 146:] = mask_error_value
        mask[203, 145:148] = mask_error_value
        mask[204, 146:148] = mask_error_value
        mask[203, 124] = mask_error_value
        mask[202, 126] = mask_error_value
        mask[243, 140] = mask_error_value
        mask[245, 133] = mask_error_value
    elif(scan_id in ['MA3886_197', 'MA3886_254', 'MA3886_273']):
        vmax_2d_img = 0.005
        vmin_2d_img = 0
        # mask dead pixels
        mask[152, 146] = mask_error_value
        mask[176, 131:133] = mask_error_value
        mask[172, 143:145] = mask_error_value
        mask[172, 147:150] = mask_error_value
        mask[173, 145:147] = mask_error_value
        mask[173, 148] = mask_error_value
        mask[191, 145] = mask_error_value
        mask[192, 145:147] = mask_error_value
        mask[193, 144] = mask_error_value
        mask[193, 146:] = mask_error_value
        mask[203, 145:148] = mask_error_value
        mask[204, 146:148] = mask_error_value
        mask[203, 124] = mask_error_value
        mask[202, 126] = mask_error_value
        mask[243, 140] = mask_error_value
        mask[245, 133] = mask_error_value

        mask[189:197, 16:25] = mask_error_value
        mask[205:215, 37:50] = mask_error_value
        mask[213:235, 50:68] = mask_error_value
        mask[226:251, 68:91] = mask_error_value

        mask[400:500, 300:400] = mask_error_value



    elif(scan_id in ['CH5314_056', 'CH5314_060']):
        vmax_2d_img = 0.03
        vmin_2d_img = -0.001
        # mask dead pixels
        mask[254:257, 94] = mask_error_value
        mask[254:257, 106] = mask_error_value
        mask[245, 131] = mask_error_value
        mask[243, 138] = mask_error_value
        mask[261, 139] = mask_error_value
        mask[271, 138] = mask_error_value
        mask[279, 139] = mask_error_value
        mask[278, 140] = mask_error_value
        mask[268, 114] = mask_error_value
        mask[204:206, 143:145] = mask_error_value
        mask[194, 142:] = mask_error_value
        mask[192, 143] = mask_error_value
        mask[193, 143:145] = mask_error_value
        mask[177, 128:130] = mask_error_value
        mask[173, 140:147] = mask_error_value
        mask[174, 142:146] = mask_error_value
        mask[92, 99] = mask_error_value
        mask[79:81, 102] = mask_error_value
        mask[77, 103:105] = mask_error_value
        #mask[255:261, :] = mask_error_value



    elif(scan_no in [392, 397, 410, 415, 418, 427]):
        vmax_2d_img = 0.2
        # mask dead pixels
        mask[176, 131] = mask_error_value
        mask[176, 132] = mask_error_value
        mask[173, 145] = mask_error_value
        mask[173, 148] = mask_error_value
        mask[191, 145] = mask_error_value
        mask[192, 145] = mask_error_value
        mask[193, 144] = mask_error_value
        mask[193, 146:] = mask_error_value
        mask[203, 145:148] = mask_error_value
        mask[204, 146:148] = mask_error_value
        mask[203, 124] = mask_error_value
        mask[202, 126] = mask_error_value
    elif(scan_no in [264]):
        vmax_2d_img = 0.02
        # mask dead pixels
        mask[255:258, 94] = mask_error_value
        mask[255:258, 106] = mask_error_value
        mask[246, 131] = mask_error_value
        mask[244, 138] = mask_error_value
        mask[262, 139] = mask_error_value
        mask[272, 138] = mask_error_value
        mask[280, 139] = mask_error_value
        mask[281, 150] = mask_error_value
        mask[288, 159] = mask_error_value
        mask[204, 142:145] = mask_error_value
        mask[205, 143:145] = mask_error_value
    elif(scan_no in [162, 163, 183, 185, 263, 276, 277, 278, 309, 310]):
        vmax_2d_img = 0.4
        print 'No dead pixel mask set for this scan.'
    elif(scan_no in [258]):
        vmax_2d_img = 0.1
        print 'No dead pixel mask set for this scan.'
    elif(scan_no in [255]):
        vmax_2d_img = 0.002
        print 'No dead pixel mask set for this scan.'
    else:
        vmax_2d_img = 0.02
        print 'No dead pixel mask set for this scan.'

    bool_mask = np.array(mask, dtype=np.bool)
    grid_intensity[~bool_mask] = 0


    if(0):
        plt.ioff()
        plt.figure()
        plt.imshow(grid_intensity-(mask-1), vmin=vmin_2d_img, vmax=vmax_2d_img, interpolation='None')

        #plt.figure()
        #plt.imshow(Intensity, vmin=0, vmax=vmax_2d_img, interpolation='None')

        dd = 5
        plt.figure()
        plt.plot(np.sum(grid_intensity[190-dd:190+dd, :], axis=0))
        plt.plot([0,515], [0.2, 0.2])

        plt.figure()
        plt.plot(np.sum(grid_intensity[:, 138-dd:138+dd], axis=1))
        plt.plot([0,515], [0.19, 0.19])

        plt.show()

    cut_ip = np.sum(grid_intensity[scan_info[scan_id].int_lim_oop[0]:scan_info[scan_id].int_lim_oop[1], :], axis=0)
    cut_oop = np.sum(grid_intensity[:, scan_info[scan_id].int_lim_ip[0]:scan_info[scan_id].int_lim_ip[1]], axis=1)


    #cut_ip = gaussian_filter(cut_ip, sigma=10)
    #cut_oop = gaussian_filter(cut_oop, sigma=10)

    if(is_zap_scan):
        #(x0, FWHM, area, lfrac, bg_slope, bg_offset)
        bounds = ((-1, 0.005, 0, 0, -1e3, 0),
                  (10, 0.2, 1e3, 1, 1e3, 1e3))
        if(frame_no == 0):
            guess_ip = (grid_q_ip[0,scan_info[scan_id].beamtime_config.cenx], 0.01, 50, 0.5, 0, 500)
            guess_oop = (grid_q_oop[scan_info[scan_id].beamtime_config.ceny, 0], 0.01, 50, 0.5, 0, 500)
        else:
            guess_ip = popt_ip
            guess_oop = popt_oop
    else:
        bounds = ((-1, 0.005, 0, 0, -1e3, 0),
                  (10, 0.2, 1e3, 1, 1e3, 1e3))
        #bounds = ((-1, 0.01, 0, 0, -1e-100, -1e-100),
        #          (10, 0.2, 1, 1, 1e-100, 1e-100))

        if(frame_no == 0):
            guess_ip = (grid_q_ip[0,scan_info[scan_id].beamtime_config.cenx], 0.01, 0.1, 0.5, 0, 0)
            guess_oop = (grid_q_oop[scan_info[scan_id].beamtime_config.ceny, 0], 0.01, 0.1, 0.5, 0, 0)
            if(scan_id == 'MA3886_273'):
                guess_ip = (grid_q_ip[0,scan_info[scan_id].beamtime_config.cenx], 0.1, 1e-4, 0, 0, 0)
                guess_oop = (grid_q_oop[scan_info[scan_id].beamtime_config.ceny, 0], 0.1, 1e-4, 0, 0, 0)
        else:
            guess_ip = popt_ip
            guess_oop = popt_oop

    # create mask to exclude Au peak from the fit
    mask_ip = np.zeros(cut_ip.shape, dtype=np.bool)
    for lims in scan_info[scan_id].fit_lim_ip:
        for j in xrange(lims[0], lims[1]):
            mask_ip[j] = 1



    # first fit to determine positions
    if(0):
        try:
            #popt_ip, pcov_ip = opt.curve_fit(model, grid_q_ip[0,scan_info[scan_id].fit_lim_ip[0]:scan_info[scan_id].fit_lim_ip[1]], cut_ip[scan_info[scan_id].fit_lim_ip[0]:scan_info[scan_id].fit_lim_ip[1]], p0=guess_ip, bounds=bounds)
            popt_ip_1, pcov_ip_1 = opt.curve_fit(model, grid_q_ip[0,mask_ip], cut_ip[mask_ip], p0=guess_ip, bounds=bounds)
        except:
            popt_ip_1, pcov_ip_1 = (np.zeros(len(guess_ip)), None)

        try:
            popt_oop_1, pcov_oop_1 = opt.curve_fit(model, grid_q_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1], 0], cut_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1]], p0=guess_oop, bounds=bounds)
        except:
            popt_oop_1, pcov_oop_1 = (np.zeros(len(guess_ip)), None)
    else:
        popt_ip_1, pcov_ip_1 = opt.curve_fit(model, grid_q_ip[0,mask_ip], cut_ip[mask_ip], p0=guess_ip, bounds=bounds)
        popt_oop_1, pcov_oop_1 = opt.curve_fit(model, grid_q_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1], 0], cut_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1]], p0=guess_oop, bounds=bounds, max_nfev=10000)

    print popt_ip_1
    print popt_oop_1



    cen_q_ip = np.argmin(np.abs(grid_q_ip[0, :]-popt_ip_1[0]))
    cen_q_oop = np.argmin(np.abs(grid_q_oop[:, 0]-popt_oop_1[0]))
    d_ip_pix = 5#5#2#1
    d_oop_pix = 5#5#2#1


    ####################################################################################
    # TODO: fix this, mask has to be applied properly
    ###################################################################################

    # normalize cut intensity to the number of good pixels
    cut_ip_contributions = np.sum(mask[cen_q_oop-d_oop_pix:cen_q_oop+d_oop_pix+1, :], axis=0)
    cut_ip = np.sum(grid_intensity[cen_q_oop-d_oop_pix:cen_q_oop+d_oop_pix+1, :], axis=0)
    cut_ip /= cut_ip_contributions
    cut_ip = cut_ip[mask_ip]
    cut_ip_q = grid_q_ip[0,mask_ip]
    cut_ip_mask = cut_ip == np.nan
    cut_ip = cut_ip[~cut_ip_mask]
    cut_ip_q = cut_ip_q[~cut_ip_mask]

    cut_oop_contributions = np.sum(mask[:, cen_q_ip-d_ip_pix:cen_q_ip+d_ip_pix+1], axis=1)
    print cut_oop_contributions
    cut_oop = np.sum(grid_intensity[:, cen_q_ip-d_ip_pix:cen_q_ip+d_ip_pix+1], axis=1)
    print cut_oop
    cut_oop /= cut_oop_contributions
    cut_oop = cut_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1]]
    cut_oop_q = grid_q_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1], 0]
    cut_oop_mask = cut_oop == np.nan
    cut_oop = cut_oop[~cut_oop_mask]
    cut_oop_q = cut_oop_q[~cut_oop_mask]

    popt_ip = popt_ip_1
    pcov_ip = pcov_ip_1
    popt_oop = popt_oop_1
    pcov_oop = pcov_oop_1

    #print cut_oop, cut_oop_q



    # second fit with small integration windows
    fit_2_times = 1
    if(fit_2_times):

        if(0):
            try:
                popt_ip, pcov_ip = opt.curve_fit(model, grid_q_ip[0,mask_ip], cut_ip[mask_ip], p0=guess_ip, bounds=bounds)
            except:
                popt_ip, pcov_ip = (np.zeros(len(guess_ip)), None)

            try:
                popt_oop, pcov_oop = opt.curve_fit(model, grid_q_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1], 0], cut_oop[scan_info[scan_id].fit_lim_oop[0]:scan_info[scan_id].fit_lim_oop[1]], p0=guess_oop, bounds=bounds)
            except:
                popt_oop, pcov_oop = (np.zeros(len(guess_ip)), None)
        else:
            popt_ip, pcov_ip = opt.curve_fit(model, cut_ip_q, cut_ip, p0=guess_ip, bounds=bounds)
            popt_oop, pcov_oop = opt.curve_fit(model, cut_oop_q, cut_oop, p0=guess_oop, bounds=bounds)


    # 2D-frame
    ###################################################################################
    fig.clear()
    ax = fig.add_subplot(111)
    #ax.pcolormesh(q_ip, q_oop, Intensity, vmin=0, vmax=0.2)
    ax.pcolormesh(grid_q_ip, grid_q_oop, grid_intensity, vmin=0, vmax=vmax_2d_img)
    ax.plot([grid_q_ip[0,0], grid_q_ip[0,-1]], [grid_q_oop[scan_info[scan_id].int_lim_oop[0],0], grid_q_oop[scan_info[scan_id].int_lim_oop[0],0]], '-', color='grey')
    ax.plot([grid_q_ip[0,0], grid_q_ip[0,-1]], [grid_q_oop[scan_info[scan_id].int_lim_oop[1],0], grid_q_oop[scan_info[scan_id].int_lim_oop[1],0]], '-', color='grey')
    ax.plot([grid_q_ip[0,scan_info[scan_id].int_lim_ip[0]], grid_q_ip[0,scan_info[scan_id].int_lim_ip[0]]], [grid_q_oop[0,0], grid_q_oop[-1,0]], '-', color='grey')
    ax.plot([grid_q_ip[0,scan_info[scan_id].int_lim_ip[1]], grid_q_ip[0,scan_info[scan_id].int_lim_ip[1]]], [grid_q_oop[0,0], grid_q_oop[-1,0]], '-', color='grey')

    if(fit_2_times):
        ax.plot([grid_q_ip[0,0], grid_q_ip[0,-1]], [grid_q_oop[cen_q_oop-d_oop_pix,0], grid_q_oop[cen_q_oop-d_oop_pix,0]], '-', color='r')
        ax.plot([grid_q_ip[0,0], grid_q_ip[0,-1]], [grid_q_oop[cen_q_oop+d_oop_pix,0], grid_q_oop[cen_q_oop+d_oop_pix,0]], '-', color='r')
        ax.plot([grid_q_ip[0, cen_q_ip-d_ip_pix], grid_q_ip[0, cen_q_ip-d_ip_pix]], [grid_q_oop[0,0], grid_q_oop[-1,0]], '-', color='r')
        ax.plot([grid_q_ip[0, cen_q_ip+d_ip_pix], grid_q_ip[0, cen_q_ip+d_ip_pix]], [grid_q_oop[0,0], grid_q_oop[-1,0]], '-', color='r')

    ax.set_xlabel(r'$q_\parallel$ / $\AA^{-1}$', fontsize=20)
    ax.set_ylabel(r'$q_\perp$ / $\AA^{-1}$', fontsize=20)
    fig.text(0.15, 0.95, 'E = '+str(potential[frame_no])+' V', fontsize=20)
    #ax.plot([2.549], [3.605], 'kx', mew=5, ms=20)

    #ax.set_xlim([2.5, 2.56])
    #ax.set_ylim([3.28, 3.36])


    fig.tight_layout()
    fig.canvas.draw()
    if(save_all_frames):
        fig.savefig('plots/%s/all_imgs/map_2d_%s.png'%(scan_id, str(frame_no).zfill(3)), dpi=300, bbox_inches='tight')

    # Cut ip
    ###################################################################################
    fig_cut_ip.clear()
    ax_cut_ip = fig_cut_ip.add_subplot(111)
    ax_cut_ip.plot(cut_ip_q, cut_ip)
    ax_cut_ip.plot(grid_q_ip[0, :], model(grid_q_ip[0,:],*popt_ip))
    #ax_cut_ip.plot(grid_q_ip[0, :], model(grid_q_ip[0,:],*popt_ip_1)/popt_ip_1[2]*popt_ip[2])
    ax_cut_ip.plot([grid_q_ip[0,scan_info[scan_id].int_lim_ip[0]], grid_q_ip[0,scan_info[scan_id].int_lim_ip[0]]], [np.min(cut_ip), np.max(cut_ip)], '-', color='grey')
    ax_cut_ip.plot([grid_q_ip[0,scan_info[scan_id].int_lim_ip[1]], grid_q_ip[0,scan_info[scan_id].int_lim_ip[1]]], [np.min(cut_ip), np.max(cut_ip)], '-', color='grey')
    for pos in np.array(scan_info[scan_id].fit_lim_ip).flatten():
        ax_cut_ip.plot([grid_q_ip[0,pos], grid_q_ip[0,pos]], [np.min(cut_ip), np.max(cut_ip)], '-', color='r')

    #ax_cut_ip.plot([grid_q_ip[0,scan_info[scan_id].fit_lim_ip[0]], grid_q_ip[0,scan_info[scan_id].fit_lim_ip[0]]], [np.min(cut_ip), np.max(cut_ip)], '-', color='r')
    #ax_cut_ip.plot([grid_q_ip[0,scan_info[scan_id].fit_lim_ip[1]], grid_q_ip[0,scan_info[scan_id].fit_lim_ip[1]]], [np.min(cut_ip), np.max(cut_ip)], '-', color='r')
    ax_cut_ip.set_xlabel(r'$q_\parallel$ / $\AA^{-1}$', fontsize=20)
    ax_cut_ip.set_ylabel(r'Intensity / a.u.', fontsize=20)
    fig_cut_ip.tight_layout()
    fig_cut_ip.canvas.draw()
    if(save_all_frames):
        fig_cut_ip.savefig('plots/%s/all_imgs/cut_ip_%s.png'%(scan_id, str(frame_no).zfill(3)), dpi=300, bbox_inches='tight')

    # Cut oop
    ###################################################################################
    fig_cut_oop.clear()
    ax_cut_oop = fig_cut_oop.add_subplot(111)
    ax_cut_oop.plot(cut_oop_q, cut_oop)
    ax_cut_oop.plot(grid_q_oop[:, 0], model(grid_q_oop[:, 0],*popt_oop))
    #ax_cut_oop.plot(grid_q_oop[:, 0], model(grid_q_oop[:, 0],*popt_oop_1))
    ax_cut_oop.plot([grid_q_oop[scan_info[scan_id].int_lim_oop[0], 0], grid_q_oop[scan_info[scan_id].int_lim_oop[0], 0]], [np.min(cut_oop), np.max(cut_oop)], '-', color='grey')
    ax_cut_oop.plot([grid_q_oop[scan_info[scan_id].int_lim_oop[1], 0], grid_q_oop[scan_info[scan_id].int_lim_oop[1], 0]], [np.min(cut_oop), np.max(cut_oop)], '-', color='grey')
    ax_cut_oop.plot([grid_q_oop[scan_info[scan_id].fit_lim_oop[0], 0], grid_q_oop[scan_info[scan_id].fit_lim_oop[0], 0]], [np.min(cut_oop), np.max(cut_oop)], '-', color='r')
    ax_cut_oop.plot([grid_q_oop[scan_info[scan_id].fit_lim_oop[1], 0], grid_q_oop[scan_info[scan_id].fit_lim_oop[1], 0]], [np.min(cut_oop), np.max(cut_oop)], '-', color='r')
    ax_cut_oop.set_xlabel(r'$q_\perp$ / $\AA^{-1}$', fontsize=20)
    ax_cut_oop.set_ylabel(r'Intensity / a.u.', fontsize=20)
    fig_cut_oop.tight_layout()
    fig_cut_oop.canvas.draw()
    if(save_all_frames):
        fig_cut_oop.savefig('plots/%s/all_imgs/cut_oop_%s.png'%(scan_id, str(frame_no).zfill(3)), dpi=300, bbox_inches='tight')


    data.potential.append(potential[frame_no])
    data.current_density.append(current[frame_no]/(np.pi*(0.2)**2))
    if(is_zap_scan):
        data.Time.append(f[scan_no].mottime)
    else:
        data.Time.append(f[scan_no].Epoch[frame_no]-f[scan_no].Epoch[0])

    data.pcov_ip.append(pcov_ip)
    data.pcov_oop.append(pcov_oop)

    data.cen_ip.append(popt_ip[0])
    data.FWHM_ip.append(popt_ip[1])
    data.amp_ip.append(popt_ip[2])
    data.lfrac_ip.append(popt_ip[3])
    data.bg_slope_ip.append(popt_ip[4])
    data.bg_offset_ip.append(popt_ip[5])

    data.cen_oop.append(popt_oop[0])
    data.FWHM_oop.append(popt_oop[1])
    data.amp_oop.append(popt_oop[2])
    data.lfrac_oop.append(popt_oop[3])
    data.bg_slope_oop.append(popt_oop[4])
    data.bg_offset_oop.append(popt_oop[5])

    plt.show()

    if(0):
        # d oop
        #########################################################################################
        fig_d_oop.clear()
        ax_d_oop = fig_d_oop.add_subplot(111)
        ax_d_oop.set_xlabel(r'E / V$_Ag/AgCl$', fontsize=20)
        ax_d_oop.set_ylabel(r'd$_\perp$ / nm', fontsize=20)
        d = 2*np.pi / np.array(data.FWHM_oop) / 10
        E = np.array(data.potential)
        ax_d_oop.plot(E, d)

        plt.show()


    if(0):
        import time
        time.sleep(100)


np.savez('data/%s.npz'%(scan_id), potential=data.potential, current_density=data.current_density, Time=data.Time, pcov_ip=data.pcov_ip, pcov_oop=data.pcov_oop,
         cen_ip=data.cen_ip, FWHM_ip=data.FWHM_ip, amp_ip=data.amp_ip, lorfact_ip=data.lfrac_ip, bg_slope_ip=data.bg_slope_ip, bg_offset_ip=data.bg_offset_ip,
         cen_oop=data.cen_oop, FWHM_oop=data.FWHM_oop, amp_oop=data.amp_oop, lorfact_oop=data.lfrac_oop, bg_slope_oop=data.bg_slope_oop, bg_offset_oop=data.bg_offset_oop)

plt.ioff()
plt.show()
