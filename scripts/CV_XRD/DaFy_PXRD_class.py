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
from DataFilterPool import create_mask, save_data_pxrd
from FitEnginePool import fit_pot_profile,backcor
from util.XRD_tools import reciprocal_space_v3 as rsp
from util.UtilityFunctions import image_generator_bkg
from util.UtilityFunctions import scan_generator
from util.UtilityFunctions import nexus_image_loader
from util.UtilityFunctions import extract_vars_from_config
from util.UtilityFunctions import show_status_bar_2
#make compatibility of py 2 and py 3#
if (sys.version_info > (3, 0)):
    raw_input = input

class run_app(object):
    def __init__(self):
        self.stop = True
        self.conf_file = None
        self.data = {}
        self.int_range = []
        self.int_range_bkg = []
        self.data_path = os.path.join(DaFy_path,'data')
        self.conf_path_temp = os.path.join(DaFy_path,'config','config_p23_pxrd_new.ini')
        self.run()

    def run(self, config = 'C:\\apps\\DaFy_P23\\config\\config_p23_pxrd_new.ini'):
        #extract global vars from config
        if config == None:
            #self.kwarg_global = extract_vars_from_config(self.conf_file, section_var ='Global')
            pass
        else:
            self.conf_file = config
            self.kwarg_global = extract_vars_from_config(self.conf_file, section_var ='Global')
        
        for each in self.kwarg_global:
            setattr(self,each,self.kwarg_global[each])

        if self.time_scan:
            self.plot_pxrd = plot_pxrd_profile_time_scan
        else:
            self.plot_pxrd = plot_pxrd_profile

        #pars lib for everything else
        self.kwarg_image = extract_vars_from_config(self.conf_file, section_var = 'Image_Loader')
        self.kwarg_mask = extract_vars_from_config(self.conf_file,section_var = 'Mask')
        self.kwarg_bkg = extract_vars_from_config(self.conf_file,section_var = 'Background_Subtraction')

        #recal clip_boundary and cen(you need to remove the edges)
        self.ver_offset = self.clip_width['ver']
        self.hor_offset = self.clip_width['hor']
        self.clip_boundary = {"ver":[self.ver_offset,self.dim_detector[0]-self.ver_offset],"hor":[self.hor_offset,self.dim_detector[1]-self.hor_offset]}
        self.cen_clip = [self.cen[0]-self.ver_offset,self.cen[1]-self.hor_offset]

        #init peak fit, bkg subtraction and reciprocal space and image loader instance
        self.img_loader = nexus_image_loader(clip_boundary = self.clip_boundary, kwarg = self.kwarg_image)
        self.create_mask_new = create_mask(kwarg = self.kwarg_mask)

        #build generator funcs
        self._scans = scan_generator(scans = self.scan_nos)
        self._images = image_generator_bkg(self._scans,self.img_loader,self.create_mask_new)

    def run_script(self):
        try:
            img = next(self._images)
            self.img = img
            self.merge_data_image_loader()
            self.fit_background()
            self._merge_data_bkg(tweak = False)
            return True
        except:
            return False

    def run_update(self):
        self.fit_background()
        self._merge_data_bkg(tweak = True)

    def merge_data_image_loader(self):
        if self.img_loader.scan_number not in self.data:
            if not self.time_scan:
                self.data[self.img_loader.scan_number] = {'2theta':[],'intensity':[],'2theta_previous':[],'intensity_previous':[]}
            else:
                self.data[self.img_loader.scan_number] = {'2theta':[],'intensity':[],'frame_number':[],'potential':[],'current':[]}
                for i in range(len(self.delta_segment_time_scan)):
                    self.data[img_loader.scan_number]['intensity_peak{}'.format(i+1)]=[]
        else:
            if self.time_scan:
                self.data[self.img_loader.scan_number]['frame_number'].append(self.img_loader.frame_number)
                self.data[self.img_loader.scan_number]['potential'].append(self.img_loader.potential)
                self.data[self.img_loader.scan_number]['current'].append(self.img_loader.current)
            else:
                pass

    def fit_background_old(self):
        int_range = np.sum(self.img[:,:],axis = 1)
        #bkg subtraction
        n=np.array(range(len(int_range)))
        #by default the contamination rate is 25%
        #the algorithm may fail if the peak cover >40% of the cut profile
        bkg_n = int(len(n)/4)
        y_sorted = list(copy.deepcopy(int_range))
        y_sorted.sort()
        std_bkg =np.array(y_sorted[0:bkg_n*3]).std()/(max(y_sorted)-min(y_sorted))
        t3=time.time()
        int_range[np.argmin(int_range)] =  int_range[np.argmin(int_range)+1]#???what for???
        int_range_bkg, *discard = backcor(range(len(int_range)),int_range,\
                                            ord_cus=self.kwarg_bkg['ord_cus_s'],\
                                            s=std_bkg*self.kwarg_bkg['ss_factor'],fct=self.kwarg_bkg['fct'])
        self.int_range = int_range-int_range_bkg
        self.int_range_bkg = list(int_range_bkg)


    def fit_background(self):
        self.int_range = np.sum(self.img[:,:],axis = 1)
        self.int_range_bkg=self.int_range*0

    def _merge_data_bkg(self, tweak = False):
        #run this after fit_background
        t0=time.time()
        #data_temp = pd.DataFrame(self.data[self.img_loader.scan_number])
        if not self.time_scan:
            delta = self.img_loader.motor_angles['delta']
            delta_range = list(np.round(delta + np.arctan((self.cen[0] - np.array(range(self.dim_detector[0]-self.ver_offset*2)))*self.ps/self.sd)/np.pi*180, 4))
            #overlap_start_index = np.argmin(abs(min(delta_range)-np.array(self.data[self.img_loader.scan_number]['2theta'])))
            #overlap_end_index = np.argmin(abs(np.array(delta_range)-max(self.data[self.img_loader.scan_number]['2theta'])))
            #append results

            time_step1=0
            time_step2=0
            time_step3=0
            if not tweak:
                '''
                index_vec=[]
                for each in delta_range:
                    if each in self.data[self.img_loader.scan_number]['2theta']:
                        index_vec.append(self.data[self.img_loader.scan_number]['2theta'].index(each))
                    else:
                        index_vec.append(None)
                index_vec_modify, index_vec_append=[],[]
                index_vec_modify_original_data = []
                for j in range(len(index_vec)):
                    if index_vec[j]==None:
                        index_vec_append.append(j)
                    else:
                        index_vec_modify.append(j)
                        index_vec_modify_original_data.append(index_vec[j])
                theta_temp = np.array(self.data[self.img_loader.scan_number]['2theta'])
                intensity_temp = np.array(self.data[self.img_loader.scan_number]['intensity'])
                intensity_temp[index_vec_modify_original_data]=0.5*intensity_temp[index_vec_modify_original_data]+0.5*self.int_range[index_vec_modify]
                intensity_temp=np.append(intensity_temp,self.int_range[index_vec_append])
                theta_temp=np.append(theta_temp,np.array(delta_range)[index_vec_append])
                '''
                time_0=time.time()-t0
                for j in delta_range:
                    t1=time.time()
                    if j in self.data[self.img_loader.scan_number]['2theta']:
                        t2=time.time()
                        jj = self.data[self.img_loader.scan_number]['2theta'].index(j)
                        time_step1=time_step1+t2-t1

                        self.data[self.img_loader.scan_number]['intensity'][jj] = 0.5*self.data[self.img_loader.scan_number]['intensity'][jj] + 0.5*self.int_range[delta_range.index(j)]
                        time_step2=time_step2+time.time()-t2
                    else:
                        self.data[self.img_loader.scan_number]['2theta'].append(j)
                        self.data[self.img_loader.scan_number]['intensity'].append(self.int_range[delta_range.index(j)])
                        time_step3=time_step3+time.time()-t1

                self.data[self.img_loader.scan_number]['2theta_previous'] = copy.deepcopy(self.data[self.img_loader.scan_number]['2theta'])
                self.data[self.img_loader.scan_number]['intensity_previous'] = copy.deepcopy(self.data[self.img_loader.scan_number]['intensity'])
                #print('test',time_0,time_step1,time_step2,time_step3)
            else:
                for j in delta_range:
                    if j in self.data[self.img_loader.scan_number]['2theta_previous']:
                        jj = self.data[self.img_loader.scan_number]['2theta_previous'].index(j)
                        self.data[self.img_loader.scan_number]['intensity'][jj] = 0.5*self.data[self.img_loader.scan_number]['intensity_previous'][jj] + 0.5*self.int_range[delta_range.index(j)]
                    else:
                        #self.data[self.img_loader.scan_number]['2theta'].append(j)
                        self.data[self.img_loader.scan_number]['intensity'].append(self.int_range[delta_range.index(j)])
            #self.data[self.img_loader.scan_number]=data_temp
            

        else:
            if len(self.data[self.img_loader.scan_number]['2theta']) == 0:
                delta = self.img_loader.motor_angles['delta']
                delta_range = np.round(delta + np.arctan((self.cen[0] - np.array(range(self.dim_detector[0]-self.ver_offset*2)))*self.ps/self.sd)/np.pi*180, 3)
                self.data[self.img_loader.scan_number]['2theta'] = delta_range
                self.data[self.img_loader.scan_number]['intensity'] = self.int_range
            else:#time scan: same 2theta value for each frame
                self.data[self.img_loader.scan_number]['intensity'] = np.array(self.int_range)
            k=0
            for each_segment in self.delta_segment_time_scan:
                k = k+1
                index_left = np.argmin(np.abs(np.array(self.data[self.img_loader.scan_number]['2theta'])-each_segment[0]))
                index_right = np.argmin(np.abs(np.array(self.data[self.img_loader.scan_number]['2theta'])-each_segment[1]))
                if not tweak:
                    self.data[self.img_loader.scan_number]['intensity_peak{}'.format(k)].append(np.array(self.int_range)[min([index_left,index_right]):max([index_left,index_right])].sum())
                else:
                    self.data[self.img_loader.scan_number]['intensity_peak{}'.format(k)][-1] = np.array(self.int_range)[min([index_left,index_right]):max([index_left,index_right])].sum()

    def _merge_data_bkg_old(self, tweak = False):
        #run this after fit_background
        if not self.time_scan:
            delta = self.img_loader.motor_angles['delta']
            delta_range = list(np.round(delta + np.arctan((self.cen[0] - np.array(range(self.dim_detector[0]-self.ver_offset*2)))*self.ps/self.sd)/np.pi*180, 3))
            #overlap_start_index = np.argmin(abs(min(delta_range)-np.array(self.data[self.img_loader.scan_number]['2theta'])))
            #overlap_end_index = np.argmin(abs(np.array(delta_range)-max(self.data[self.img_loader.scan_number]['2theta'])))
            #append results
            time_step1=0
            time_step2=0
            time_step3=0
            if not tweak:
                for j in delta_range:
                    t1=time.time()
                    if j in self.data[self.img_loader.scan_number]['2theta']:
                        t2=time.time()
                        time_step1=time_step1+t2-t1
                        #jj = np.where(np.array(self.data[self.img_loader.scan_number]['2theta'])==j)[0][0]
                        jj = self.data[self.img_loader.scan_number]['2theta'].index(j)
                        time_step2=time_step2+time.time()-t2
                        
                        self.data[self.img_loader.scan_number]['intensity'][jj] = 0.5*self.data[self.img_loader.scan_number]['intensity'][jj] + 0.5*self.int_range[delta_range.index(j)]
                        
                    else:
                        self.data[self.img_loader.scan_number]['2theta'].append(j)
                        self.data[self.img_loader.scan_number]['intensity'].append(self.int_range[delta_range.index(j)])
                        time_step3=time_step3+time.time()-t1
                self.data[self.img_loader.scan_number]['2theta_previous'] = copy.deepcopy(self.data[self.img_loader.scan_number]['2theta'])
                self.data[self.img_loader.scan_number]['intensity_previous'] = copy.deepcopy(self.data[self.img_loader.scan_number]['intensity'])
            else:
                for j in delta_range:
                    if j in self.data[self.img_loader.scan_number]['2theta_previous']:
                        jj = self.data[self.img_loader.scan_number]['2theta_previous'].index(j)
                        self.data[self.img_loader.scan_number]['intensity'][jj] = 0.5*self.data[self.img_loader.scan_number]['intensity_previous'][jj] + 0.5*self.int_range[delta_range.index(j)]
                    else:
                        #self.data[self.img_loader.scan_number]['2theta'].append(j)
                        self.data[self.img_loader.scan_number]['intensity'].append(self.int_range[delta_range.index(j)])
            print(time_step1,time_step2,time_step3)

        else:
            if len(self.data[self.img_loader.scan_number]['2theta']) == 0:
                delta = self.img_loader.motor_angles['delta']
                delta_range = np.round(delta + np.arctan((self.cen[0] - np.array(range(self.dim_detector[0]-self.ver_offset*2)))*self.ps/self.sd)/np.pi*180, 3)
                self.data[self.img_loader.scan_number]['2theta'] = delta_range
                self.data[self.img_loader.scan_number]['intensity'] = self.int_range
            else:#time scan: same 2theta value for each frame
                self.data[self.img_loader.scan_number]['intensity'] = np.array(self.int_range)
            k=0
            for each_segment in self.delta_segment_time_scan:
                k = k+1
                index_left = np.argmin(np.abs(np.array(self.data[self.img_loader.scan_number]['2theta'])-each_segment[0]))
                index_right = np.argmin(np.abs(np.array(self.data[self.img_loader.scan_number]['2theta'])-each_segment[1]))
                if not tweak:
                    self.data[self.img_loader.scan_number]['intensity_peak{}'.format(k)].append(np.array(self.int_range)[min([index_left,index_right]):max([index_left,index_right])].sum())
                else:
                    self.data[self.img_loader.scan_number]['intensity_peak{}'.format(k)][-1] = np.array(self.int_range)[min([index_left,index_right]):max([index_left,index_right])].sum()

    def save_data_file(self,path):
        #to be finished
        df = pd.DataFrame(self.data)
        if path.endswith('.xlsx'):          
            df.to_excel(path)
        else:
            df.to_excel(path+'.xlsx')
        #save data
        save_data_pxrd(data=int_intensity[scan_number], scan_number=scan_number, path=DaFy_path, time_scan = time_scan)



if __name__ == "__main__":
    run_app()
