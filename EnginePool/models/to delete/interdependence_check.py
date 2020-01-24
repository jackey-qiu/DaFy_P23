import os,sys
os.chdir("C:\\apps\\genx\\")
import model,diffev,time,fom_funcs
import filehandling as io
import glob
import numpy as np

mod = model.Model()
config = io.Config()
opt = diffev.DiffEv()
io.load_gx(sys.argv[1],mod,opt,config)
fom_containor=np.array(np.zeros([1,10])[0:0])

for i in range(mod.parameters.get_len_rows()):
    if (mod.parameters.get_value(i,0)!='')&(mod.parameters.get_value(i,2)==True):
        left,right=mod.parameters.get_value(i,3),mod.parameters.get_value(i,4)
        val_init=mod.parameters.get_value(i,1)
        temp_fom=[]
        for val in np.arange(left,right,(right-left)/5):
            mod.parameters.set_value(i,1,val)
            mod.simulate()
            temp_fom.append(mod.fom)
        fom_containor=np.append(fom_containor,[list(np.arange(left,right,(right-left)/5))+temp_fom],axis=0)
        mod.parameters.set_value(i,1,val_init)
np.savetxt(sys.argv[2],fom_containor,"%5.3f")

