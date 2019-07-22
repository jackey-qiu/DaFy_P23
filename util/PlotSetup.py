import numpy as np

class ScanInfo():

    def __init__(self, scan_id, data, structure_lattice, HKL_position, scan_direction_ranges, color, analog_data_filename=None, analog_data_plot_range=[0,-1], ids_filename=None, analog_data_interval=None, scan_label=None, pH=13):
        self.scan_id = scan_id
        self.data = data
        self.structure_lattice = structure_lattice
        self.HKL_position = HKL_position
        self.scan_direction_ranges = scan_direction_ranges
        self.color = color
        self.analog_data_filename = analog_data_filename
        self.analog_data_plot_range = analog_data_plot_range
        self.ids_filename = ids_filename
        self.analog_data_interval = analog_data_interval
        self.scan_label = scan_label
        self.pH = pH

class ScanInfoContainer(dict):
    def __init__(self):
        return
    def add(self, scan_id, *args, **kwargs):
        self[scan_id] = ScanInfo(scan_id, *args, **kwargs)

#####################################################################
def RHE(E_AgAgCl, pH=13):
    # electrode is 3.4 M KCl
    return 0.205 + E_AgAgCl + 0.059*pH

def POT(E_AgAgCl, plot_vs_RHE=False, pH=13):
    if(plot_vs_RHE):
        return RHE(E_AgAgCl, pH)
    else:
        return E_AgAgCl

def find_tick_ticklables(x, num_ticks = 5, endpoint = True, dec_place = 1):
    x_ticks = np.linspace(min(x),max(x),num = num_ticks, endpoint = endpoint)
    steps_x = np.round(x_ticks[1]-x_ticks[0], dec_place)
    first_x = np.round(x_ticks[0],dec_place)
    x_ticks =[first_x+i*steps_x for i in range(num_ticks)]
    x_tick_labels = [str(np.round(each_x,dec_place)) for each_x in x_ticks]
    x_ticks = [float(each) for each in x_tick_labels]
    return x_ticks, x_tick_labels

def select_cycle_new(pot_result_frame_triple, bin_mode = 'select', return_cycle = 0, n_pots_rebin = 120, tailing_cut = 50):
    #bin_mode = 'average' or 'select'
    #   'average':average every bin_level data points
    #   'select': select every bin_level data point
    #plot_mode = 'CV' or 'pot_step'
    #   'CV': plot CV data
    #   'pot_step': plot potential step data
    pot_bin, result_bin = [], []
    pot, result, frame_numbers = pot_result_frame_triple
    pot, result, frame_numbers = pot[0:len(pot)-tailing_cut], result[0:len(result)-tailing_cut], frame_numbers[0:len(frame_numbers)-tailing_cut]
    pot = np.around(pot,decimals=6)
    pot_offset = pot.max()-pot.min()
    pot_keys_number = n_pots_rebin
    pot_keys = np.array([i*abs(pot[0]-pot[1])+pot.min() for i in range(pot_keys_number)])
    result_lib ={}
    for key in pot_keys:
        result_lib[(key,0)] = []
        result_lib[(key,1)] = []
    for i in range(len(pot)):
        current_pot = pot[i]
        up_or_down = 0
        if i ==0:
            up_or_down = int(pot[i+1]>pot[i])
        elif i == len(pot)-1:
            up_or_down = int(pot[i-1]<pot[i])
        else:
            if pot[i]>pot[i+1] and (frame_numbers[i+1]-frame_numbers[i])==1:
                up_or_down = 0
            elif pot[i]>pot[i-1] and  (frame_numbers[i]-frame_numbers[i-1])==1:
                up_or_down = 1
            else:
                up_or_down = 0
        current_key = (pot_keys[np.argmin(np.abs(pot_keys - current_pot))],up_or_down)
        result_lib[current_key].append(result[i])

    if bin_mode == 'select':
        for key in pot_keys:
            if result_lib[(key,0)]!=[]:
                try:
                    result_bin.append(result_lib[(key,0)][return_cycle])
                except:
                    result_bin.append(result_lib[(key,0)][0])
                pot_bin.append(key)
            else:
                pass
        for key in pot_keys:
            if result_lib[(key,1)]!=[]:
                try:
                    result_bin.append(result_lib[(key,1)][return_cycle])
                except:
                    result_bin.append(result_lib[(key,1)][0])
                pot_bin.append(key)
            else:
                pass
    elif bin_mode == 'average':
        for key in pot_keys:
            if result_lib[(key,0)]!=[]:
                result_bin.append(np.mean(result_lib[(key,0)]))
                pot_bin.append(key)
            else:
                pass
        for key in pot_keys:
            if result_lib[(key,1)]!=[]:
                result_bin.append(np.mean(result_lib[(key,1)]))
                pot_bin.append(key)
            else:
                pass
    return [0,int(len(pot_bin)/2), -2], np.array(pot_bin), np.array(result_bin)


def select_cycle(pot_result_couple, bin_level = 1, bin_mode = 'select', return_cycle = 2, pot_step = 0.10, plot_mode='CV'):
    #bin_mode = 'average' or 'select'
    #   'average':average every bin_level data points
    #   'select': select every bin_level data point
    #plot_mode = 'CV' or 'pot_step'
    #   'CV': plot CV data
    #   'pot_step': plot potential step data
    pot, result = pot_result_couple
    pot = np.around(pot,decimals=6)
    index_of_valley_pot = []
    min_pot_value = np.min(pot)
    temp_index = []
    #get valley points
    for i in range(1,len(pot)):
        if abs(pot[i]-min_pot_value)<pot_step:
            if temp_index!=[]:
                if (i-temp_index[-1])<10:
                    temp_index.append(i)
                else:
                    index_of_valley_pot.append(temp_index[int(len(temp_index)/2)])
                    temp_index = [i]
            else:
                temp_index.append(i)
    if temp_index!=[]:
        index_of_valley_pot.append(temp_index[int(len(temp_index)/2)])

    if plot_mode == 'pot_step':
        index_of_valley_pot = [0, len(pot)]
        bin_level = 1# you dont want to bin for pot_step dataset
        return [0, len(pot), -2], pot, result
    if return_cycle > len(index_of_valley_pot)-2:
        return_cycle = len(index_of_valley_pot)-2

    total_cycles = len(index_of_valley_pot)
    length_of_each_cycle = int(len(pot)/total_cycles)

    pot_partial, result_partial = pot[length_of_each_cycle*return_cycle:length_of_each_cycle*(return_cycle+1)],\
                                 result[length_of_each_cycle*return_cycle:length_of_each_cycle*(return_cycle+1)]
    # pot_partial, result_partial = pot[index_of_valley_pot[return_cycle]:index_of_valley_pot[return_cycle+1]],\
                                 # result[index_of_valley_pot[return_cycle]:index_of_valley_pot[return_cycle+1]]
    # pot_partial, result_partial = pot[0:index_of_valley_pot[return_cycle]],\
                                  # result[0:index_of_valley_pot[return_cycle]]
    def _bin_array(data, bin_level, bin_mode):
        if bin_mode =='average':
            data_bin =np.array(data[0::bin_level])
            for i in range(1,bin_level):
                temp = np.array(data[i::bin_level])
                if len(temp)<len(data_bin):
                    temp = np.append(temp, np.zeros(len(data_bin) - len(temp))+temp[-1])
                    data_bin = data_bin + temp
            return data_bin/bin_level
        elif bin_mode =='select':
            return np.array(data[0::bin_level])


    if bin_level>1:
        pot_bin = _bin_array(pot_partial, bin_level,bin_mode)
        result_bin = _bin_array(result_partial, bin_level,bin_mode)
        return [0, int((index_of_valley_pot[return_cycle+1]-index_of_valley_pot[return_cycle])/(bin_level*2)),-2],pot_bin, result_bin
    else:
        return [0, int((index_of_valley_pot[return_cycle+1]-index_of_valley_pot[return_cycle])/2),-2],pot_partial, result_partial

def set_max_to_0(data_list,slice_index):
    return np.array(data_list)[slice_index[0]:slice_index[1]]-np.array(data_list).max()

def extract_ids_file(file_path,which_cycle=3):
    data = []
    data_lines =[]
    current_cycle = 0
    with open(file_path,encoding="ISO-8859-1") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            line = lines[i]
            if line.startswith('primary_data'):
                print(current_cycle)
                current_cycle=current_cycle+1
                if current_cycle == which_cycle:
                    for j in range(i+3,i+3+int(lines[i+2].rstrip())):
                        data.append([float(each) for each in lines[j].rstrip().rsplit()])
                    break
                else:
                    pass
            else:
                pass
    return np.array(data)[:,0], np.array(data)[:,1]

def strain_ip(q_ip, HKL, lattice):
    q_bulk = lattice.q(HKL)
    q_ip_bulk = np.sqrt(q_bulk[0]**2 + q_bulk[1]**2)
    return (q_ip_bulk/q_ip - 1.0)*100.

def strain_ip_with_uncertainty(q_ip, HKL, lattice, uncertainty_q_ip):
    q_bulk = lattice.q(HKL)
    q_ip_bulk = np.sqrt(q_bulk[0]**2 + q_bulk[1]**2)
    _strain_ip = (q_ip_bulk/q_ip - 1.0)*100.
    uncedrtainty_strain_ip = np.abs(q_ip_bulk/q_ip**2*100*uncertainty_q_ip)
    return (_strain_ip, uncedrtainty_strain_ip)

def strain_oop(q_oop, HKL, lattice):
    return (lattice.q(HKL)[2]/q_oop - 1.0)*100.

def strain_oop_with_uncertainty(q_oop, HKL, lattice, uncertainty_q_oop):
    q_oop_bulk = lattice.q(HKL)[2]
    _strain_oop = (q_oop_bulk/q_oop - 1.0)*100.
    uncertainty_strain_oop = np.abs(q_oop_bulk/q_oop**2*100*uncertainty_q_oop)
    return  (_strain_oop, uncertainty_strain_oop)
#####################################################################################################################################

