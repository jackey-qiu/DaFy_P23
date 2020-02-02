import models.structure_tools.sxrd_dafy as model
from models.utils import UserVars
from datetime import datetime
import numpy as np
import sys,pickle,__main__,os
import batchfile.locate_path as batch_path
import dump_files.locate_path as output_path
import models.domain_creator as domain_creator
import accessory_functions.make_par_table.make_parameter_table_GenX_hematite_rcut as make_grid
import accessory_functions.data_formating.formate_xyz_to_vtk as xyz
from copy import deepcopy
import models.setup_domain_hematite_rcut as setup_domain_hematite_rcut
from models.structure_tools import tool_box
from models.structure_tools import sorbate_tool


COUNT_TIME=False
if COUNT_TIME:t_0=datetime.now()

#setting slabs##
wal=0.551#wavelength of x ray
unitcell = model.UnitCell(3.615, 3.615, 3.615, 90, 90, 90)
inst = model.Instrument(wavel = wal, alpha = 2.0)
bulk = model.Slab()
surface_1 =  model.Slab(c = 1.0)
surface_2 =  model.Slab(c = 1.0)
sorbate_1 = model.Slab(c = 1.0)
sorbate_2 = model.Slab(c = 1.0)
rgh=UserVars()
rgh.new_var('beta', 0.0)#roughness factor
rgh.new_var('mu',1)#liquid film thickness
scales=['scale_CTR']
for scale in scales:
    rgh.new_var(scale,1.)
rgh_co2_1=UserVars()
rgh_co2_1.new_var('r1',1.2)
rgh_co2_1.new_var('r2',1.2)
rgh_co2_1.new_var('r3',1.2)
rgh_co2_1.new_var('delta1',30)
rgh_co2_1.new_var('delta2',30)
rgh_co2_1.new_var('gamma',0)

rgh_co2_2=UserVars()
rgh_co2_2.new_var('r1',1.2)
rgh_co2_2.new_var('r2',1.2)
rgh_co2_2.new_var('r3',1.2)
rgh_co2_2.new_var('delta1',30)
rgh_co2_2.new_var('delta2',30)
rgh_co2_2.new_var('gamma',0)

rgh_wt = UserVars()
rgh_wt.new_var('wt_domain1',1)

#set up experimental constant#
re = 2.818e-5#electron radius
kvect=2*np.pi/wal#k vector
Egam = 6.626*(10**-34)*3*(10**8)/wal*10**10/1.602*10**19#energy in ev
LAM=1.5233e-22*Egam**6 - 1.2061e-17*Egam**5 + 2.5484e-13*Egam**4 + 1.6593e-10*Egam**3 + 1.9332e-06*Egam**2 + 1.1043e-02*Egam
exp_const = 4*kvect/LAM
auc=unitcell.a*unitcell.b*np.sin(unitcell.gamma)

##############################################end of main setup zone############################################
#                                                                                                              #
#                                                                                                              #
#                           You seldomly need to touch script lines hereafter!!!                               #
#                                                                                                              #
#                                                                                                              #
#                                                                                                              #
################################################################################################################
#depository path for output files(structure model files(.xyz,.cif), optimized values (CTR,RAXR,E_Density) for plotting
output_file_path=output_path.module_path_locator()

################################################build up ref domains############################################
#add atoms for bulk and two ref domains (ref_domain1<half layer> and ref_domain2<full layer>)                  #
#In those two reference domains, the atoms are ordered according to first hight (z values), then y values      #
#it is a super surface structure by stacking the surface slab on bulk slab, the repeat vector was counted      #
################################################################################################################
batch_path_head=batch_path.module_path_locator()
tool_box.add_atom_in_slab(bulk,os.path.join(batch_path_head,'Cu100_bulk.str'))
tool_box.add_atom_in_slab(surface_1,os.path.join(batch_path_head,'Cu100_surface_1.str'))
tool_box.add_atom_in_slab(sorbate_1,os.path.join(batch_path_head,'Cu100_sorbate_OCCO_1.str'))

tool_box.add_atom_in_slab(surface_2,os.path.join(batch_path_head,'Cu100_surface_2.str'))
tool_box.add_atom_in_slab(sorbate_2,os.path.join(batch_path_head,'Cu100_sorbate_OCCO_2.str'))

co2_1 = sorbate_tool.CarbonOxygenMotif(sorbate_1,['C2', 'O1', 'O2'], 'C1', 1.,0,0,flat_down_index=[])
co2_1.set_coordinate_all(gamma_list=[180,180,0],delta_list = [0,30,30], r_list=[1.5,1.2,1.2], new_anchor_list = [None,'C2',None])

co2_2 = sorbate_tool.CarbonOxygenMotif(sorbate_2,['C2','O1', 'O2'], 'C1', 1.,0,0,flat_down_index=[])
co2_2.set_coordinate_all(gamma_list=[180,180,0],delta_list = [0,30,30], r_list=[1.5,1.2,1.2], new_anchor_list = [None,'C2',None])

atm_gp_sorbate_1 = model.AtomGroup()
for id in ['O1', 'O2', 'C1', 'C2']:
    atm_gp_sorbate_1.add_atom(sorbate_1,id)

atm_gp_sorbate_2 = model.AtomGroup()
for id in ['O1', 'O2', 'C1', 'C2']:
    atm_gp_sorbate_2.add_atom(sorbate_2,id)

atm_gp_surface_1st_layer_1 = model.AtomGroup()
for id in ['Cu_top1_2','Cu_top2_2']:
    atm_gp_surface_1st_layer_1.add_atom(surface_1,id)

atm_gp_surface_1st_layer_2 = model.AtomGroup()
for id in ['Cu_top1_2','Cu_top2_2']:
    atm_gp_surface_1st_layer_2.add_atom(surface_2,id)

atm_gp_surface_2nd_layer_1 = model.AtomGroup()
for id in ['Cu1_2','Cu2_2']:
    atm_gp_surface_2nd_layer_1.add_atom(surface_1,id)

atm_gp_surface_2nd_layer_2 = model.AtomGroup()
for id in ['Cu3_2','Cu4_2']:
    atm_gp_surface_2nd_layer_2.add_atom(surface_2,id)
###################################fitting function part##########################################
VARS=vars()#pass local variables to sim function
if COUNT_TIME:t_1=datetime.now()

sorbate_syms_1 = [model.SymTrans([[1,0],[0,1]]), model.SymTrans([[0,1],[-1,0]]), model.SymTrans([[1,0],[0,1]]), model.SymTrans([[0,-1],[1,0]])][0:1]
sorbate_syms_2 = [model.SymTrans([[1,0],[0,1]]), model.SymTrans([[0,1],[-1,0]]), model.SymTrans([[1,0],[0,1]]), model.SymTrans([[0,-1],[1,0]])][0:1]
#sorbate_syms = []
sample = model.Sample(inst, bulk, {'domain1':{'slab':surface_1, 'sorbate':sorbate_1,'wt':rgh_wt.wt_domain1,'sorbate_sym':sorbate_syms_1},\
                                   'domain2':{'slab':surface_2, 'sorbate':sorbate_2,'wt':1-rgh_wt.wt_domain1,'sorbate_sym':sorbate_syms_2} }, unitcell)

def Sim(data,VARS=VARS):
    #co2.make_cif_file('D://test.cif')

    #for command in commands:eval(command)
    co2_1.set_coordinate_all(gamma_list=[180+rgh_co2_1.gamma,180+rgh_co2_1.gamma, rgh_co2_1.gamma],delta_list = [0,rgh_co2_1.delta1,rgh_co2_1.delta2], r_list=[rgh_co2_1.r1,rgh_co2_1.r2,rgh_co2_1.r3], new_anchor_list = [None,'C2',None])
    co2_2.set_coordinate_all(gamma_list=[180+rgh_co2_2.gamma,180+rgh_co2_2.gamma,rgh_co2_2.gamma],delta_list = [0,rgh_co2_2.delta1,rgh_co2_2.delta2], r_list=[rgh_co2_2.r1,rgh_co2_2.r2,rgh_co2_2.r3], new_anchor_list = [None,'C2',None])
    #co2_1.make_xyz_file('D://test.xyz')
    # sample.make_xyz_file('D://test.xyz')
    sample.domain['domain1']['wt']=rgh_wt.wt_domain1
    if rgh_wt.wt_domain1>1:
        sample.domain['domain2']['wt']=0
    else:
        sample.domain['domain2']['wt']=1-rgh_wt.wt_domain1

    for each in sorbate_syms_1:
        each.set_t([atm_gp_surface_1st_layer_1.getdx(),atm_gp_surface_1st_layer_1.getdy()])
    for each in sorbate_syms_2:
        each.set_t([atm_gp_surface_1st_layer_2.getdx(),atm_gp_surface_1st_layer_2.getdy()])

    VARS=VARS
    F =[]
    bv=0
    bv_container={}
    fom_scaler=[]
    beta=rgh.beta
    SCALES=[getattr(rgh,scale) for scale in scales]
    total_wt=0

    #for i in range(DOMAIN_NUMBER):
        #grap wt for each domain and cal the total wt
    #    vars()['wt_domain'+str(int(i+1))]=VARS['rgh_domain'+str(int(i+1))].wt
    #    total_wt=total_wt+vars()['wt_domain'+str(int(i+1))]



    if COUNT_TIME:t_2=datetime.now()

    #cal structure factor for each dataset in this for loop
    #fun to deal with the symmetrical shape of 10,30 and 20L rod at positive and negative sides
    def formate_hkl(h_,k_,x_):
        new_h,new_k,new_x=[],[],[]
        if np.around(h_[0],0) in [1,2,3] and np.around(k_[0],0)==0:
            for iii in range(len(x_)):
                if x_[iii]>0:
                    new_h.append(-h_[iii])
                    new_k.append(-k_[iii])
                    new_x.append(-x_[iii])
                else:
                    new_h.append(h_[iii])
                    new_k.append(k_[iii])
                    new_x.append(x_[iii])
            return  np.array(new_h),np.array(new_k),np.array(new_x)
        else:
            return np.array(h_),np.array(k_),np.array(x_)

    i=0
    for data_set in data:
        f=np.array([])
        h = data_set.extra_data['h']
        k = data_set.extra_data['k']
        x = data_set.x
        y = data_set.extra_data['Y']
        LB = data_set.extra_data['LB']
        dL = data_set.extra_data['dL']


        if data_set.use:
            if data_set.x[0]>100:#doing RAXR calculation(x is energy column typically in magnitude of 10000 ev)
                h_,k_,y_=formate_hkl(h,k,y)
                rough = (1-beta)/((1-beta)**2 + 4*beta*np.sin(np.pi*(y-LB)/dL)**2)**0.5#roughness model, double check LB and dL values are correctly set up in data file
                f=sample.cal_structure_factor_hematite_RAXR(i,VARS,RAXR_FIT_MODE,RESONANT_EL_LIST,RAXR_EL,h_, k_, y_, x, E0, F1F2,SCALES,rough)
                F.append(abs(f))
                fom_scaler.append(1)
                i+=1
            else:#doing CTR calculation (x is perpendicular momentum transfer L typically smaller than 15)
                h_,k_,x_=formate_hkl(h,k,x)
                rough = (1-beta)/((1-beta)**2 + 4*beta*np.sin(np.pi*(x-LB)/dL)**2)**0.5#roughness model, double check LB and dL values are correctly set up in data file
                if h[0]==0 and k[0]==0:#consider layered water only for specular rod if existent
                    q=np.pi*2*unitcell.abs_hkl(h,k,x)
                    pre_factor=(np.exp(-exp_const*rgh.mu/q))*(4*np.pi*re/auc)*3e6
                    f = pre_factor*SCALES[0]*rough*sample.calc_f_all(h, k, x)
                else:
                    f = rough*sample.calc_f_all(h, k, x)
                F.append(abs(f*f))
                fom_scaler.append(1)
        else:
            if x[0]>100:
                i+=1
            f=np.zeros(len(y))
            F.append(f)
            fom_scaler.append(1)

    #some ducumentation about using this script#
    print_help_info=False
    if print_help_info:
        setup_domain_hematite_rcut.print_help_doc()
    #do this in shell 'model.script_module.setup_domain_hematite_rcut.print_help_doc()' to get help info

    #output how fast the code is running#
    if COUNT_TIME:t_3=datetime.now()
    if COUNT_TIME:
        print("It took "+str(t_1-t_0)+" seconds to setup")
        print("It took "+str(t_2-t_1)+" seconds to calculate bv weighting")
        print("It took "+str(t_3-t_2)+" seconds to calculate structure factor")

    #you may play with the weighting rule by setting eg 2**bv, 5**bv for the wt factor, that way you are pushing the GenX to find a fit btween a good fit (low wt factor) and a reasonable fit (high wt factor)
    return F,1+1*bv,fom_scaler
