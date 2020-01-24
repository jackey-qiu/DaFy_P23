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

rows=mod.parameters.get_len_rows()
i_false=0
for row in range(rows):
	if (mod.parameters.get_value(row,2)==True):
		mod.parameters.set_value(row,2,False)
		if mod.parameters.get_value(row+1,2)==False:
			i_false=row+1
			break
for row in range(i_false,i_false+10):
	try:
		mod.parameters.set_value(row,2,True)
	except:
		break
io.save_gx(sys.argv[1],mod,opt,config)