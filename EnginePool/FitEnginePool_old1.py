import scipy.optimize as opt
import matplotlib
matplotlib.use("tkAgg")
import matplotlib.pyplot as plt
from util.UtilityFunctions import collect_args
import sys, copy
from numpy import dtype
from scipy.interpolate import griddata
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter
# import p23_tools_debug as p23
import subprocess
from DataFilterPool import *
from util.Reciprocal_space_tools.HKLVlieg import Crystal, printPos, UBCalculator, VliegAngles, printPos_prim, vliegDiffracAngles
from PyMca5.PyMcaPhysics import SixCircle

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

def calculate_UB_matrix_p23(lattice_constants, energy, or0_angles, or1_angles,or0_hkl,or1_hkl):
    return p23.cal_UB(lattice_constants, energy, or0_angles, or1_angles,or0_hkl, or1_hkl)

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


class XRD_Peak_Fitting(object):
    def __init__(self, img, mask, q_ip, q_oop, cut_offset={'hor':[50,10],'ver':[50,10]}, data_range_offset = {'hor':[100,50],'ver':[100,50]},peak_center= None,model=model,fit_p0=[], fit_bounds={'hor':[],'ver':[]}):
        self.img = img
        self.mask = mask
        self.model = model
        self.first_frame = True
        self.q_ip = q_ip
        self.q_oop = q_oop
        self.cut_offset = cut_offset
        self.data_range_offset = data_range_offset
        self.peak_center = peak_center
        self.fit_data = {'hor':{'x':[],'y':[]},'ver':{'x':[],'y':[]}}
        self.fit_results = {'hor':[],'ver':[]}
        self.fit_results_plot = {'hor':[],'ver':[]}
        self.fit_p0 = fit_p0
        self.fit_p0_2 = fit_p0
        self.fit_bounds = fit_bounds
        self.fit()

    def reset_fit(self,img, use_first_fit_for_pos=True, check = False, level = 0.05, tweak_mode = False):
        self.first_frame = False
        self.img = img
        check_result = self.fit(use_first_fit_for_pos, check, level, tweak_mode)
        return check_result

    def update_bounds(self, cen_oop, cen_ip):
        if cen_oop>self.fit_bounds['ver'][1][0]:
            self.fit_bounds['ver'][1][0]=cen_oop+0.1
        elif cen_oop<self.fit_bounds['ver'][0][0]:
            self.fit_bounds['ver'][0][0]=cen_oop-0.1

        if cen_ip>self.fit_bounds['hor'][1][0]:
            self.fit_bounds['hor'][1][0]=cen_ip+0.1
        elif cen_ip<self.fit_bounds['hor'][0][0]:
            self.fit_bounds['hor'][0][0]=cen_ip-0.1

    def cut_profile_from_2D_img_around_center(self, img, cut_offset = {'hor':10, 'ver':20}, data_range_offset = {'hor':50, 'ver':50}, center_index = None, sum_result = True):

        def _cut_profile_from_2D_img(img, cut_range, cut_direction, sum_result=True):
            if cut_direction=='horizontal':
                if sum_result:
                    return np.sum(img[cut_range[0]:cut_range[1],:],axis = 0)
                else:
                    return np.average(img[cut_range[0]:cut_range[1],:],axis = 0)

            elif cut_direction=='vertical':
                if sum_result:
                    return np.sum(img[:,cut_range[0]:cut_range[1]],axis = 1)
                else:
                    return np.average(img[:,cut_range[0]:cut_range[1]],axis = 1)
        size = img.shape
        f = lambda x, y: [max([x-y,0]), x+y]
        if center_index == None:
            center_index = [int(each/2) for each in size]
        data_range = {'hor':f(center_index[1],data_range_offset['hor']),'ver':f(center_index[0], data_range_offset['ver'])}
        cut_range = {'hor': f(center_index[0],cut_offset['hor']),'ver': f(center_index[1], cut_offset['ver'])}
        # print(data_range['hor'])
        cut = {'hor':_cut_profile_from_2D_img(img, cut_range['hor'], cut_direction ='horizontal', sum_result = sum_result)[data_range['hor'][0]:data_range['hor'][1]],
                'ver':_cut_profile_from_2D_img(img, cut_range['ver'], cut_direction ='vertical', sum_result = sum_result)[data_range['ver'][0]:data_range['ver'][1]]}
        return cut

    def fit(self, use_first_fit_for_pos = True, check=False, level = 0.05, tweak_mode = False):
        self.img[self.mask ==0]=0
        self.fit_data = {'hor':{'x':[],'y':[]},'ver':{'x':[],'y':[]}}
        fit_ip_0, fom_ip_0 = None, None
        fit_oop_0, fom_oop_0 = None, None
        peak_locating_step = True
        if tweak_mode:
            i_range = [len(self.cut_offset['hor'])-1]
        else:
            i_range = range(len(self.cut_offset['hor']))
        #first cut with large window
        for i in i_range:
            fit_p0 = {0:self.fit_p0,1:self.fit_p0_2}
            center_index = {0:self.peak_center, 1:self.peak_center, 2:self.peak_center}
            cut_offset_temp =dict([(key,value[i]) for key, value in self.cut_offset.items()])
            data_range_offset_temp =dict([(key,value[i]) for key, value in self.data_range_offset.items()])
            cut = self.cut_profile_from_2D_img_around_center(self.img,cut_offset_temp,data_range_offset_temp, center_index[i], sum_result = True)
            cut_mask = self.cut_profile_from_2D_img_around_center(self.mask,cut_offset_temp,data_range_offset_temp, center_index[i], sum_result = False)
            cut_ip_q = self.cut_profile_from_2D_img_around_center(self.q_ip,cut_offset_temp, data_range_offset_temp, center_index[i], sum_result = False)['hor']
            cut_oop_q = self.cut_profile_from_2D_img_around_center(self.q_oop,cut_offset_temp, data_range_offset_temp, center_index[i], sum_result = False)['ver']
            #normalize the dead pixel contribution
            cut = {'hor':cut['hor']/cut_mask['hor'],'ver':cut['ver']/cut_mask['ver']}

            #check nan and inf values
            index_used = {'hor':~(np.isnan(cut['hor'])|np.isinf(cut['hor'])), 'ver':~(np.isnan(cut['ver'])|np.isinf(cut['ver']))}
            if self.first_frame and i==0:
                self.fit_p0['ver'][0] = self.q_oop[self.peak_center[0],0]
                self.fit_p0['hor'][0] = self.q_ip[0,self.peak_center[1]]
            fit_ip,fom_ip = opt.curve_fit(f=self.model, xdata=cut_ip_q[index_used['hor']], ydata=cut['hor'][index_used['hor']], p0 = fit_p0[i]['hor'], bounds = self.fit_bounds['hor'], max_nfev = 10000)
            fit_oop,fom_oop = opt.curve_fit(f=self.model, xdata=cut_oop_q[index_used['ver']], ydata=cut['ver'][index_used['ver']], p0 = fit_p0[i]['ver'], bounds = self.fit_bounds['ver'], max_nfev = 10000)

            self.fit_data['hor']['x'].append(cut_ip_q[index_used['hor']])
            self.fit_data['hor']['y'].append(cut['hor'][index_used['hor']])
            self.fit_data['ver']['x'].append(cut_oop_q[index_used['ver']])
            self.fit_data['ver']['y'].append(cut['ver'][index_used['ver']])
            if tweak_mode:
                self.fit_data['hor']['x'].append(cut_ip_q[index_used['hor']])
                self.fit_data['hor']['y'].append(cut['hor'][index_used['hor']])
                self.fit_data['ver']['x'].append(cut_oop_q[index_used['ver']])
                self.fit_data['ver']['y'].append(cut['ver'][index_used['ver']])

            #update the peak center, but not change the other par values
            if i==0:
                self.fit_p0['ver'] = fit_oop
                self.fit_p0['hor'] = fit_ip
                fit_ip_0, fom_ip_0 = fit_ip, fom_ip
                fit_oop_0, fom_oop_0 = fit_oop, fom_oop
            else:
                self.fit_p0_2['ver'] = fit_oop
                self.fit_p0_2['hor'] = fit_ip
            self.update_bounds(fit_oop[0],fit_ip[0])
            peak_center_ = [np.argmin(np.abs(self.q_oop[:,0]-fit_oop[0])),np.argmin(np.abs(self.q_ip[0,:]-fit_ip[0]))]
            self.peak_center = peak_center_
            # if check and abs(np.array(self.peak_center)-np.array(peak_center_)).sum()<2000:
                # self.peak_center = peak_center_
            # elif check and abs(np.array(self.peak_center)-np.array(peak_center_)).sum()>=2000:
                # if i == 0:
                    # peak_locating_step = False
            print(self.peak_center,'at',fit_oop[0],fit_ip[0])
            # if check and abs(np.array(self.peak_center)-np.array(peak_center_)).sum()<2000:
                # self.peak_center = peak_center_
            # elif check and abs(np.array(self.peak_center)-np.array(peak_center_)).sum()>=2000:
                # if i == 0:
                    # peak_locating_step = False
            # quit()

        # quit()
        #finish the fit and update the fit par values
        # self.update_bounds(fit_oop[0],fit_ip[0])
        self.fit_results_plot['hor'] = [copy.deepcopy(fit_ip), copy.deepcopy(fom_ip)]
        self.fit_results_plot['ver'] = [copy.deepcopy(fit_oop), copy.deepcopy(fom_oop)]
        # if use_first_fit_for_pos and peak_locating_step:
            # fit_ip[0], fom_ip[0] = fit_ip_0[0], fom_ip_0[0]
            # fit_oop[0], fom_oop[0] = fit_oop_0[0], fom_oop_0[0]
        if use_first_fit_for_pos and peak_locating_step and (not tweak_mode):
            fit_ip[0], fom_ip[0] = fit_ip_0[0], fom_ip_0[0]
            fit_oop[0], fom_oop[0] = fit_oop_0[0], fom_oop_0[0]
        def _check(old, new, level = level):
            check_result = bool((abs((np.array(old)[0:2] - np.array(new)[0:2])/np.array(old)[0:2])>level).sum())
            # print(check_result)
            return check_result
        if check:
            if (_check(self.fit_results['hor'][0],fit_ip) | _check(self.fit_results['ver'][0],fit_oop))==False:
                self.fit_results['hor'] = [fit_ip, fom_ip]
                self.fit_results['ver'] = [fit_oop, fom_oop]
                return True
            else:
                return False
        else:
            self.fit_results['hor'] = [fit_ip, fom_ip]
            self.fit_results['ver'] = [fit_oop, fom_oop]
            return True
        # print fit_ip,fit_oop

    def save_data(self,data):
        # data = container
        pcov_ip = self.fit_results['hor'][1]
        pcov_oop = self.fit_results['ver'][1]

        popt_ip = self.fit_results['hor'][0]
        popt_oop = self.fit_results['ver'][0]

        data['pcov_ip'].append(pcov_ip)
        data['pcov_oop'].append(pcov_oop)

        data['cen_ip'].append(popt_ip[0])
        data['FWHM_ip'].append(popt_ip[1])
        data['amp_ip'].append(popt_ip[2])
        data['lfrac_ip'].append(popt_ip[3])
        data['bg_slope_ip'].append(popt_ip[4])
        data['bg_offset_ip'].append(popt_ip[5])

        data['cen_oop'].append(popt_oop[0])
        data['FWHM_oop'].append(popt_oop[1])
        data['amp_oop'].append(popt_oop[2])
        data['lfrac_oop'].append(popt_oop[3])
        data['bg_slope_oop'].append(popt_oop[4])
        data['bg_offset_oop'].append(popt_oop[5])
        # print popt_oop[0],data['cen_oop']
        return data

    def update_data(self,data):
        # data = container
        pcov_ip = self.fit_results['hor'][1]
        pcov_oop = self.fit_results['ver'][1]

        popt_ip = self.fit_results['hor'][0]
        popt_oop = self.fit_results['ver'][0]

        data['pcov_ip'][-1]=pcov_ip
        data['pcov_oop'][-1]=pcov_oop

        data['cen_ip'][-1]=popt_ip[0]
        data['FWHM_ip'][-1]=popt_ip[1]
        data['amp_ip'][-1]=popt_ip[2]
        data['lfrac_ip'][-1]=popt_ip[3]
        data['bg_slope_ip'][-1]=popt_ip[4]
        data['bg_offset_ip'][-1]=popt_ip[5]

        data['cen_oop'][-1]=popt_oop[0]
        data['FWHM_oop'][-1]=popt_oop[1]
        data['amp_oop'][-1]=popt_oop[2]
        data['lfrac_oop'][-1]=popt_oop[3]
        data['bg_slope_oop'][-1]=popt_oop[4]
        data['bg_offset_oop'][-1]=popt_oop[5]
        # print popt_oop[0],data['cen_oop']
        return data

class Reciprocal_Space_Mapping():
    def __init__(self, img, E_keV=19.5, cen=(234,745), pixelsize=(0.055,0.055), sdd=714, UB=[],motor_angles=None):
        self.img = img
        self.intensity = None
        self.E_keV = E_keV
        self.wavelength = 12.39854*self.E_keV
        self.k0 = 2.*np.pi/self.wavelength
        self.cen = cen
        self.pixelsize = pixelsize
        self.sdd = sdd
        self.UB=UB
        self.motor_angles = motor_angles
        self.q=None
        self.grid_intensity = None
        # self.prepare_frame()
        # self.get_grid_q()

    def update_img(self, img,UB=None, motor_angles=None,update_q = True):
        self.img = img
        if UB!=None:
            self.UB =UB
        if motor_angles!=None:
            self.motor_angles = motor_angles
        self.prepare_frame()
        if update_q:
            self.get_grid_q()

    def prepare_frame(self, norm_mon=True, norm_transm=True,trans='attenpos',mon='avg_beamcurrent'):
        transm_= 1
        mon_= 1
        th_= self.motor_angles['mu']
        gam_= self.motor_angles['delta']
        del_= self.motor_angles['gamma']
        #the chi and phi values are arbitrary in the fio file, should be set to the same values as the ones that are usd to cal UB matrix(all 0 so far)
        phi_= self.motor_angles['phi']
        chi_= self.motor_angles['chi']
        mu_= self.motor_angles['omega_t']
        # print del_,gam_
        #first item is the incident angle (mu_ here)
        del_,gam_=np.rad2deg(vliegDiffracAngles(np.deg2rad([mu_,del_,gam_,mu_,0,0]))[1:3])
        # print del_,gam_
        intensity = self.img
        #detector dimension is (516,1556)
        #You may need to put a negative sign in front, check the rotation sense of delta and gamma motors at P23
        delta_range = np.arctan((np.arange(intensity.shape[1])-self.cen[0])*self.pixelsize[0]/self.sdd)*180/ np.pi + del_
        #the minus sign here because the column index increase towards bottom, then 0 index(top most) will give a negative gam offset
        #a minus sign in front correct this.
        gamma_range =-np.arctan((np.arange(intensity.shape[0])-self.cen[1])*self.pixelsize[1]/self.sdd)*180/ np.pi + gam_
        #polarisation correction
        # TODO: what is this doing?
        delta_grid , gamma_grid= np.meshgrid(delta_range,gamma_range)
        Pver = 1 - np.sin(delta_grid * np.pi / 180.)**2 * np.cos(gamma_grid * np.pi / 180.)**2
        intensity=np.divide(intensity,Pver)
        self.intensity = intensity
        self.vlieg_angles = {'gamma_range':gamma_range, 'delta_range':delta_range,'th_':th_, 'mu_':mu_, 'chi_':chi_, 'phi_':phi_}

    def get_HKL(self):
        for each in ['gamma_range','delta_range', 'th_', 'mu_', 'chi_', 'phi_']:
            locals()[each]=self.vlieg_angles[each]
        d = SixCircle.SixCircle()
        d.setEnergy(self.E_keV)
        d.setUB(self.UB)
        HKL = d.getHKL(delta=delta_range, theta=th_, chi=chi_, phi=phi_, mu=mu_, gamma=gamma_range, gamma_first=False)
        shape =  gamma_range.size,delta_range.size
        # shape =  delta_range.size,gamma_range.size
        H = HKL[0,:].reshape(shape)
        K = HKL[1,:].reshape(shape)
        L = HKL[2,:].reshape(shape)
        self.HKL = {'H':H, 'K':K, 'L':L}

    def get_grid_q(self):
        for each in ['gamma_range','delta_range', 'th_', 'mu_', 'chi_', 'phi_']:
            globals()[each]=self.vlieg_angles[each]
        d = SixCircle.SixCircle()
        d.setEnergy(self.E_keV)
        d.setUB(self.UB)
        Q = d.getQSurface(theta=th_, chi=chi_, phi=phi_, mu=mu_, delta=delta_range, gamma=gamma_range, gamma_first=False)
        shape =  gamma_range.size,delta_range.size
        qx = Q[0,:].reshape(shape)
        qy = Q[1,:].reshape(shape)
        qz = Q[2,:].reshape(shape)
        q_para = np.sqrt(qx**2 + qy**2)
        size =self.intensity.shape
        #shape=(vertical,horizontal),len(vertical)=size(1),len(horizontal)=size(0)
        grid_q_perp, grid_q_para = np.mgrid[np.max(qz):np.min(qz):(1.j*size[0]), np.min(q_para):np.max(q_para):(1.j*size[1])]
        grid_intensity = griddata((q_para.ravel(), qz.ravel()), self.intensity.ravel(), (grid_q_para, grid_q_perp), method='nearest')
        self.grid_intensity = grid_intensity
        self.q={'qx':qx,'qy':qy,'qz':qz,\
                'q_par':q_para, 'q_perp':qz,\
                'grid_q_par':grid_q_para,'grid_q_perp':grid_q_perp}

    def show_image(self):
        grid_q_para, grid_q_perp, grid_intensity = self.get_grid_q_in_out_plane(self.scan_no, self.frame_no,self.frame_prefix)
        plt.figure()
        # plt.imshow(grid_intensity, vmin=0, vmax=100.05)
        plt.imshow(grid_intensity,cmap='jet',vmin=0, vmax=100.05)
        plt.title("plt.imshow(grid_intensity")
        # plt.colorbar(extend='both',orientation='Vertical')
        plt.clim(0,90)
        plt.show()
