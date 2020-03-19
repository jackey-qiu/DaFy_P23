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
import accessory_functions.format_data.format_hkl as format_hkl
from copy import deepcopy
import models.setup_domain_hematite_rcut as setup_domain_hematite_rcut
from models.structure_tools import tool_box
from models.structure_tools import sorbate_tool

#--global settings--#
#/globalsetting/begin#
#/path/begin#
batch_path_head=batch_path.module_path_locator()
output_file_path=output_path.module_path_locator()
#/path/end#
#/wavelength/begin#
wal=0.551
#/wavelength/end#
#/slabnumber/begin#
num_sorbate_slabs = 2
num_surface_slabs = 2
#/slabnumber/end#
#/sorbatemotif/begin#
sorbate_motifs = ['OCCO','OCCO']
#/sorbatemotif/end#
#/expconstant/begin#
re = 2.818e-5#electron radius
kvect=2*np.pi/wal#k vector
Egam = 6.626*(10**-34)*3*(10**8)/wal*10**10/1.602*10**19#energy in ev
LAM=1.5233e-22*Egam**6 - 1.2061e-17*Egam**5 + 2.5484e-13*Egam**4 + 1.6593e-10*Egam**3 + 1.9332e-06*Egam**2 + 1.1043e-02*Egam
exp_const = 4*kvect/LAM
#/expconstant/end#
#globalsetting/end#

#--set unitcell--#
#/unitcell/begin#
unitcell = model.UnitCell(3.615, 3.615, 3.615, 90, 90, 90)
#/unitcell/end#

#--set instrument--#
#/instrument/begin#
inst = model.Instrument(wavel = wal, alpha = 2.0)
#/instrument/end#

#--set bulk slab--#
#/bulk/begin#
bulk = model.Slab()
bulk_file = 'Cu100_bulk.str'
tool_box.add_atom_in_slab(bulk,os.path.join(batch_path_head,bulk_file))
#/bulk/end#

#--set surface slabs--#
#/surfaceslab/begin#
surface_slab_head = 'Cu100_surface_'
for i in range(num_surface_slabs):
    globals()['surface_{}'.format(i+1)] = model.Slab(c = 1.0)
    tool_box.add_atom_in_slab(globals()['surface_{}'.format(i+1)],os.path.join(batch_path_head,'{}{}.str'.format(surface_slab_head, i+1)))
#/surfaceslab/end#

#--set sorbate slabs--#
#/sorbateslab/begin#
sorbate_slab_head = 'Cu100_sorbate_OCCO_'
for i in range(num_sorbate_slabs):
    globals()['sorbate_{}'.format(i+1)] = model.Slab(c = 1.0)
    tool_box.add_atom_in_slab(globals()['sorbate_{}'.format(i+1)],os.path.join(batch_path_head,'{}{}.str'.format(sorbate_slab_head,i+1)))
#/sorbateslab/end#

#--set sorbate coordinates--#
#/sorbatecoordinates/begin#
sorbate_instance_head = 'co2_'
atom_ids_sorbate = ['C1','C2', 'O1', 'O2']
gamma_list=[180,180,0]
delta_list = [0,30,30]
r_list=[1.5,1.2,1.2]
new_anchor_list = [None,'C2',None]
for i in range(num_sorbate_slabs):
    globals()['{}{}'.format(sorbate_instance_head, i+1)] = sorbate_tool.CarbonOxygenMotif(globals()['sorbate_{}'.format(i+1)],atom_ids_sorbate[1:], atom_ids_sorbate[0], 1.,0,0,flat_down_index=[])
    globals()['{}{}'.format(sorbate_instance_head, i+1)].set_coordinate_all(gamma_list = gamma_list,delta_list = delta_list, r_list=r_list, new_anchor_list = new_anchor_list)
#/sorbatecoordinates/end#

#--set rgh--#
#/rgh/begin#

#/rgh/global/begin#
rgh = UserVars()
rgh.new_var('beta', 0.0)#roughness factor
rgh.new_var('mu',1)#liquid film thickness
scales=['scale_CTR']
for scale in scales:
    rgh.new_var(scale,1.)
#/rgh/global/end#

#/rgh/sorbate/begin#
sorbate_rgh_head = 'rgh_co2_'
vars = ['r1', 'r2', 'r3', 'delta1', 'delta2', 'gamma']
ini_vals = [1.2, 1.2, 1.2, 30, 30, 0]
for i in range(num_sorbate_slabs):
    globals()['{}{}'.format(sorbate_rgh_head,i+1)] = UserVars()
    for each in vars:
        globals()['{}{}'.format(sorbate_rgh_head, i+1)].new_vars(each,ini_vals[vars.index(each)])
#/rgh/sorbate/end#

#/rgh/domain_weight/begin#
rgh_wt = UserVars()
for i in range(num_surface_slabs):
    rgh_wt.new_var('wt_domain{}'.format(i+1),1)
#/rgh/domain_weight/end#
#/rgh/end#
"C<-O->C<-O<-C
#/atmgroup/begin#
#/sorbate/begin#
for i in range(num_sorbate_slabs):
    globals()['atm_gp_sorbate_{}'.format(i+1)] = model.AtomGroup()
    for id in atom_ids_sorbate:
        globals()['atm_gp_sorbate_{}'.format(i+1)].add_atom(globals()['sorbate_{}'.format(i+1)], id)
#/sorbate/end#
#/substrate/start#
n_substrate_layers = 2
map_layer_index = {1:'1st',2:'2nd',3:'3rd',4:'4th',5:'5th',6:'6th',7:'7th',8:'8th',9:'9th',10:'10th'}
map_layer_atoms = {1:['Cu_top1_2','Cu_top2_2'],\
                   2:['Cu1_2','Cu2_2'],\
                   3:['Cu1_2','Cu2_2']}
for i in range(num_surface_slabs):
    for j in range(n_substrate_layers):
        globals()['atm_gp_surface_{}_layer_{}'.format(map_layer_index[j+1], i+1)] = model.AtomGroup()
        for id in map_layer_atoms[j+1]:
            globals()['atm_gp_surface_{}_layer_{}'.format(map_layer_index[j+1], i+1)].add_atom(globals()['surface_{}'.format(i+1)],id)
#/substrate/end#
#/atmgroup/end#

#/sorbatesym/begin#
n_sym = 1
for i in range(num_sorbate_slabs):
    globals()['sorbate_syms_{}'.format(i+1)] =  [model.SymTrans([[1,0],[0,1]]), model.SymTrans([[0,1],[-1,0]]), model.SymTrans([[1,0],[0,1]]), model.SymTrans([[0,-1],[1,0]])][0:n_sym]
#/sorbatesym/end#

#/sample/begin#
domains = {}
for i in range(num_surface_slabs):
    domains['domain{}'.format(i+1)] = {}
    domains['domain{}'.format(i+1)]['slab'] = globals()['surface_{}'.format(i+1)]
    domains['domain{}'.format(i+1)]['sorbate'] = globals()['sorbate_{}'.format(i+1)]
    domains['domain{}'.format(i+1)]['wt'] = globals()['rgh_wt']['wt_domain{}'.format(i+1)]
    domains['domain{}'.format(i+1)]['sorbate_sym'] = globals()['sorbate_syms_{}'.format(i+1)]
sample = model.Sample(inst, bulk, domains, unitcell)
#/sample/end#

def Sim(data,VARS=vars()):
    F =[]
    fom_scaler=[]
    beta=rgh.beta
    SCALES=[getattr(rgh,scale) for scale in scales]

    #/update_sorbate/begin#
    for i in range(num_sorbate_slabs):
        gamma_value = VARS['{}{}'.format(VARS['sorbate_rgh_head'],i+1)].gamma
        gamma_list = [180+gamma_value, 180+gamma_value, gamma_value]
        delta_list = [0, VARS['{}{}'.format(VARS['sorbate_rgh_head'],i+1)].delta1, VARS['{}{}'.format(VARS['sorbate_rgh_head'],i+1)].delta2]
        r_list = [VARS['{}{}'.format(VARS['sorbate_rgh_head'],i+1)].r1, VARS['{}{}'.format(VARS['sorbate_rgh_head'],i+1)].r2, VARS['{}{}'.format(VARS['sorbate_rgh_head'],i+1)].r3]
        new_anchor_list = [None,'C2',None]
        VARS['{}{}'.format(VARS['sorbate_instance_head'],i+1)].set_coordinate_all(gamma_list = gamma_list, delta_list = delta_list, r_list = r_list, new_anchor_list = new_anchor_list)
    #/update_sorbate/end#

    #normalize the domain weight to make total = 1
    wt_list = [getattr(rgh_wt, 'wt_domain{}'.format(i+1)) for i in range(num_surface_slabs)]
    total_wt = sum(wt_list)
    for i in range(num_surface_slabs):
        sample.domain['domain{}'.format(i+1)]['wt']=wt_list[i]/total_wt
    
    #update sorbate symmetry
    for i in range(num_sorbate_slabs):
        for each in VARS['sorbate_syms_{}'.format(i+1)]:
            each.set_t([VARS['atm_gp_surface_1st_layer_{}'.format(i+1)].getdx(),VARS['atm_gp_surface_1st_layer_{}'.format(i+1)].getdy()]) 

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
                #not yet implemented
                pass
            else:
                h,k,x=format_hkl(h,k,x)
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
            f=np.zeros(len(y))
            F.append(f)
            fom_scaler.append(1)

    #you may play with the weighting rule by setting eg 2**bv, 5**bv for the wt factor, that way you are pushing the GenX to find a fit btween a good fit (low wt factor) and a reasonable fit (high wt factor)
    panelty_factor = [sample.bond_distance_constraint(which_domain=i, max_distance =2.2) for i in num_surface_slabs]
    return F,1+panelty_factor,fom_scaler