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
slopes=[]
names,values,switch,low,high=[],[],[],[],[]

for i in range(mod.parameters.get_len_rows()):
    if (mod.parameters.get_value(i,0)!='')&(mod.parameters.get_value(i,2)==True):
        left,right=mod.parameters.get_value(i,3),mod.parameters.get_value(i,4)
        val_init=mod.parameters.get_value(i,1)
        temp_fom=[]
        for val in np.arange(left,right,(right-left)/5):
            mod.parameters.set_value(i,1,val)
            mod.simulate()
            temp_fom.append(mod.fom)
        fom_max,fom_min=max(temp_fom),min(temp_fom)
        fom_max_index,fom_min_index=temp_fom.index(fom_max),temp_fom.index(fom_min)
        par_max,par_min=np.arange(left,right,(right-left)/5)[fom_max_index],np.arange(left,right,(right-left)/5)[fom_min_index]
        slopes.append(abs((fom_max-fom_min)/(par_max-par_min)))
        print mod.parameters.get_value(i,0),"par_max,par_min=",par_max,par_min,"sensitivity=",slopes[-1]
        mod.parameters.set_value(i,1,val_init)
    else:
        slopes.append(-1000+i)
indexs=np.array(slopes).argsort()[::-1]
print indexs
for index in indexs:
    names.append(mod.parameters.get_value(index,0))
    values.append(mod.parameters.get_value(index,1))
    switch.append(mod.parameters.get_value(index,2))
    low.append(mod.parameters.get_value(index,3))
    high.append(mod.parameters.get_value(index,4))
for i in range(len(indexs)):
    mod.parameters.set_value(i,0,names[i])
    mod.parameters.set_value(i,1,values[i])
    mod.parameters.set_value(i,2,switch[i])
    mod.parameters.set_value(i,3,low[i])
    mod.parameters.set_value(i,4,high[i])
io.save_gx(sys.argv[2],mod,opt,config)


        
