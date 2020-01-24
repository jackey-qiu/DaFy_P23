import numpy as np
ref_id_list_1=["O1_1_0","O1_2_0","O1_3_0","O1_4_0","Fe1_4_0","Fe1_6_0","O1_5_0","O1_6_0","O1_7_0","O1_8_0","Fe1_8_0","Fe1_9_0","O1_9_0","O1_10_0","Fe1_10_0","Fe1_12_0","O1_11_0","O1_12_0",\
'O1_1_1','O1_2_1','Fe1_2_1','Fe1_3_1','O1_3_1','O1_4_1','Fe1_4_1','Fe1_6_1','O1_5_1','O1_6_1','O1_7_1','O1_8_1','Fe1_8_1','Fe1_9_1','O1_9_1','O1_10_1','Fe1_10_1','Fe1_12_1','O1_11_1','O1_12_1']
#generate library for using in next function
#arguments for this function:
#ref_id_list:reference id list for a half or full domain (here is half domain)
#N_dxdy: how many atoms will be considered for inplane movement
#N_dz, N_oc, N_u: how many surface atoms will be considered for dz, occupancy and thermol factor
#N_metals, N_OH: number of metals and oxygen ligands

def lib_creator(ref_id_list=ref_id_list_1,N_dxdy=4,N_dz=14,N_oc=6,N_u=14,N_metal=2,N_OH=0,el='Pb'):
    displace={}
    debye={}
    occupation={}
    #inplane movement for surface atoms
    for i in range(N_dxdy):
        displace[1+2*i]=(ref_id_list[i],'dx')
        displace[2+2*i]=(ref_id_list[i],'dy')
    #dz movement for surface atoms
    for i in range(N_dz/2):
        displace[N_dxdy*2+i+1]=(ref_id_list[i*2],'dz')
    #debye factor for surface atoms
    for i in range(N_u/2):
        debye[i+1]=(ref_id_list[i*2],'u')
    #occupation for surface atoms
    for i in range(N_oc/2):
        occupation[i+1]=(ref_id_list[i*2],'oc')
    #movement for sorbates
    for i in range(N_metal):
        N=max(displace.keys())
        displace[N+1]=(el+str(i+1),'dx')
        displace[N+2]=(el+str(i+1),'dy')
    for i in range(N_OH):
        N=max(displace.keys())
        displace[N+1]=('HO'+str(i+1),'dx')
        displace[N+2]=('HO'+str(i+1),'dy')
    for i in range(N_metal):
        N=max(displace.keys())
        displace[N+1]=(el+str(i+1),'dz')
    for i in range(N_OH):
        N=max(displace.keys())
        displace[N+1]=('HO'+str(i+1),'dz')
    #debye for sorbates
    for i in range(N_metal):
        N=max(debye.keys())
        debye[N+1]=(el+str(i+1),'u')
    for i in range(N_OH):
        N=max(debye.keys())
        debye[N+1]=('HO'+str(i+1),'u')
    #occupation for sorbates
    for i in range(N_metal):
        N=max(occupation.keys())
        occupation[N+1]=(el+str(i+1),'oc')
    for i in range(N_OH):
        N=max(occupation.keys())
        occupation[N+1]=('HO'+str(i+1),'oc')
    return {'displace':displace,'debye':debye,'occupation':occupation}
#here the returned lib looks like: displace, debye, occupation={serial number:(id,var)}, var can be any in list of ['dx','dy','dz','u','oc']

#this function will update current par file using the parameters from genx (so should be called inside sim function)
#You must manually fix the debye term in the returned par file (insert two items right before sorbates to represent the debye factors for freezed atoms (one for iron one for oxygen), and make sure the serial number are continued after this step.
def from_tab_to_par(domain,par_file_old,lib):
    f_par_old=open(par_file_old)
    f_par_new=open(par_file_old+'new','w')
    def _find_index(id_list,id):
        index=None
        for each_id in id_list:
            if id in each_id:
                index=id_list.index(each_id)
                break
        return index
    def _cal_dx(domain,index):
        return domain.dx1[index]+domain.dx2[index]+domain.dx3[index]
    def _cal_dy(domain,index):
        return domain.dy1[index]+domain.dy2[index]+domain.dy3[index]
    def _cal_dz(domain,index):
        return domain.dz1[index]+domain.dz2[index]+domain.dz3[index]        
    for line in f_par_old.readlines():
        items=line.rstrip().rsplit()
        if items[0] not in lib.keys():
            f_par_new.write(line)
        else:
            id,var=lib[items[0]][int(items[1])]
            index=_find_index(list(domain.id),id)
            if var=='dx':
                items[2]=str(_cal_dx(domain,index))
            elif var=='dy':
                items[2]=str(_cal_dy(domain,index))
            elif var=='dz':
                items[2]=str(_cal_dz(domain,index))
            elif var=='u':
                items[2]=str(domain.u[index])
            elif var=='oc':
                items[2]=str(domain.oc[index])
            f_par_new.write(' '.join(items)+'\n')
    return "ran successfully"
