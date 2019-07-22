import sys
import numpy as np

def cut_profile_from_2D_img(img, cut_range, cut_direction, sum_result=True):
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

def extract_subset_of_zap_scan(n_data=0, l_boundary=[0, 1], bragg_ls = [], skip_bragg_l_offset = 0.05, delta_l_norm = 0.05, min_delta_l =0.01):
    def _point_density_calculator(l, l_bragg_adjacent,delta_l_bt_adj_bragg_peak):
        density_factor = 1/np.sin((l-l_bragg_adjacent)/delta_l_bt_adj_bragg_peak*np.pi)**2
        delta_l_current_point = delta_l_norm/density_factor
        if delta_l_current_point < min_delta_l:
            delta_l_current_point = min_delta_l
        # print(density_factor,delta_l_current_point)
        return delta_l_current_point
    l_all = np.linspace(l_boundary[0],l_boundary[1],n_data)
    l_new = []
    current_l = l_boundary[0]
    l_new.append(current_l)
    partial_index = []
    while current_l<l_boundary[1]:
        adjacent_bragg_l_index = np.argmin(np.abs(np.array(bragg_ls)-current_l))
        adjacent_bragg_l = bragg_ls[adjacent_bragg_l_index]
        if current_l > adjacent_bragg_l:
            right = True
        else:
            right = False
        if right:
            next_bragg_l = bragg_ls[adjacent_bragg_l_index +1]
        else:
            next_bragg_l = bragg_ls[adjacent_bragg_l_index -1]
        delta_l_current_point = _point_density_calculator(current_l, adjacent_bragg_l, abs(next_bragg_l-adjacent_bragg_l))
        if abs(current_l+delta_l_current_point-adjacent_bragg_l)>skip_bragg_l_offset:
            l_new.append(current_l+delta_l_current_point)
        else:
            pass
        current_l+= delta_l_current_point
    for each_l in l_new:
        partial_index.append(np.argmin(np.abs(np.array(l_all)-each_l)))
    # print(len(l_new),l_new)
    return partial_index


def cut_profile_from_2D_img_around_center(img, cut_offset = {'hor':10, 'ver':20}, data_range_offset = {'hor':50, 'ver':50}, center_index = None, sum_result = True):
    size = img.shape
    f = lambda x, y: [x-y, x+y]
    if center_index == None:
        center_index = [int(each/2) for each in size]
    data_range = {'hor':f(center_index[1],data_range_offset['hor']),'ver':f(center_index[0], data_range_offset['ver'])}
    cut_range = {'hor': f(center_index[0],cut_offset['hor']),'ver': f(center_index[1], cut_offset['ver'])}
    cut = {'hor': cut_profile_from_2D_img(img, cut_range['hor'], cut_direction ='horizontal', sum_result = sum_result)[data_range['hor'][0]:data_range['hor'][1]],
            'ver': cut_profile_from_2D_img(img, cut_range['ver'], cut_direction ='vertical', sum_result = sum_result)[data_range['ver'][0]:data_range['ver'][1]]}
    return cut

def create_mask(img, img_q_par, img_q_ver, threshold = 10000, compare_method ='larger',remove_columns = [], \
                remove_rows = [], remove_pix = None, remove_q_range = {'par':[], 'ver':[]}, \
                remove_partial_range = {'point_couple':[{'p1':[2.4,3.2],'p2':[2.5,3.0]},{'p1':[2.55,3.3],'p2':[2.4,3.0]}],'pixel_width':[]}):
    #remove_partial_range, each point is of form [q_hor, q_ver]; two point couple will be used to calculate a line equation, which will be
    #coupled with pixel width to get those pixel index points to be excluded.
    #point group item must be a list of lib with only two items with keys of p1 and p2.
    mask = np.ones(np.array(img).shape)
    mask_ip_q = np.ones(np.array(img_q_par).shape)
    mask_oop_q = np.ones(np.array(img_q_ver).shape)

    def _find_pixel_index_from_q(grid_q_par, grid_q_ver, point):
        q_par_one_row = grid_q_par[0,:]
        q_ver_one_col = grid_q_ver[:,0]
        qx,qy = point
        index_point = [np.argmin(abs(q_ver_one_col - qy)),np.argmin(abs(q_par_one_row - qx))]
        return index_point

    if remove_q_range['par']!=[]:
        for each_range in remove_q_range['par']:
            mask_ip_q[np.where(np.logical_and(img_q_par>each_range[0], img_q_par<each_range[1]))] = 0
            # remove_rows = remove_rows + range(np.argmin(abs(img_q_par[:,0]- each_range[1])), np.argmin(abs(img_q_par[:,0]- each_range[0])))
            # print(np.argmin(abs(img_q_ver[:,0]- each_range[1])), np.argmin(abs(img_q_ver[:,0]- each_range[0])))
            # print img_q_ver.max(), img_q_par.max()
    if remove_q_range['ver']!=[]:
        for each_range in remove_q_range['ver']:
            mask_oop_q[np.where(np.logical_and(img_q_ver>each_range[0], img_q_ver<each_range[1]))] = 0
            # mask_oop_q[img_q_ver>each_range[0] && img_q_ver<each_range[1]] = 0
            # remove_columns = remove_columns + range(np.argmin(abs(img_q_ver[:,0]- each_range[0])), np.argmin(abs(img_q_ver[:,0]- each_range[1])))
            # print(np.argmin(abs(img_q_par[:,0]- each_range[0])), np.argmin(abs(img_q_par[:,0]- each_range[1])))
    if remove_partial_range['pixel_width']!=[]:
        for i in range(len(remove_partial_range['pixel_width'])):
            p1, p2 = remove_partial_range['point_couple'][i]['p1'], remove_partial_range['point_couple'][i]['p2']
            p1_index = _find_pixel_index_from_q(img_q_ver, img_q_par, p1)
            p2_index = _find_pixel_index_from_q(img_q_ver, img_q_par, p2)
            slope = float(p2_index[1]-p1_index[1])/float(p2_index[0]-p1_index[0])
            offset = p1_index[1] - slope*p1_index[0]
            h, w = img.shape[:2]
            j, k = np.ogrid[:h,:w]
            temp_mask = abs(j*slope+offset-k) < remove_partial_range['pixel_width'][i]
            remove_pix.extend([each for each in np.argwhere(temp_mask == True) if (each[0]>min(p1_index[0],p2_index[0])) and (each[0]<max(p1_index[0],p2_index[0]))])
    if compare_method =='larger':
        mask[img>threshold]=0
    elif compare_method =='smaller':
        mask[img<threshold]=0
    elif compare_method =='equal':
        maks[img == threshold] =0
    mask[:,remove_columns] = 0
    mask[remove_rows,:] = 0
    if remove_pix!=None:
        for each in remove_pix:
            mask[tuple(each)] = 0
    return mask*mask_ip_q*mask_oop_q

def create_mask_bkg(img, threshold = 10000, compare_method ='larger',remove_columns = [], \
                remove_rows = [], remove_pix = None, remove_xy_range = {'par':[], 'ver':[]}, \
                remove_partial_range = {'point_couple':[{'p1':[10,20],'p2':[20,30]},{'p1':[20,30],'p2':[24,35]}],'pixel_width':[]}):
    #remove_partial_range, each point is of form [q_hor, q_ver]; two point couple will be used to calculate a line equation, which will be
    #coupled with pixel width to get those pixel index points to be excluded.
    #point group item must be a list of lib with only two items with keys of p1 and p2.For each point,(x_hor_index, y_ver_index)
    mask = np.ones(np.array(img).shape)

    if remove_xy_range['par']!=[]:
        for each_range in remove_xy_range['par']:
            mask[each_range[0]:each_range[1],:] = 0
            # remove_rows = remove_rows + range(np.argmin(abs(img_q_par[:,0]- each_range[1])), np.argmin(abs(img_q_par[:,0]- each_range[0])))
            # print(np.argmin(abs(img_q_ver[:,0]- each_range[1])), np.argmin(abs(img_q_ver[:,0]- each_range[0])))
            # print img_q_ver.max(), img_q_par.max()
    if remove_xy_range['ver']!=[]:
        for each_range in remove_xy_range['ver']:
            mask[:,each_range[0]:each_range[1]] = 0
            # mask_oop_q[img_q_ver>each_range[0] && img_q_ver<each_range[1]] = 0
            # remove_columns = remove_columns + range(np.argmin(abs(img_q_ver[:,0]- each_range[0])), np.argmin(abs(img_q_ver[:,0]- each_range[1])))
            # print(np.argmin(abs(img_q_par[:,0]- each_range[0])), np.argmin(abs(img_q_par[:,0]- each_range[1])))
    if remove_partial_range['pixel_width']!=[]:
        for i in range(len(remove_partial_range['pixel_width'])):
            p1_index, p2_index = remove_partial_range['point_couple'][i]['p1'], remove_partial_range['point_couple'][i]['p2']
            slope = float(p2_index[1]-p1_index[1])/float(p2_index[0]-p1_index[0])
            offset = p1_index[1] - slope*p1_index[0]
            h, w = img.shape[:2]
            j, k = np.ogrid[:h,:w]
            temp_mask = abs(j*slope+offset-k) < remove_partial_range['pixel_width'][i]
            remove_pix.extend([each for each in np.argwhere(temp_mask == True) if (each[0]>min(p1_index[0],p2_index[0])) and (each[0]<max(p1_index[0],p2_index[0]))])
    if compare_method =='larger':
        mask[img>threshold]=0
    elif compare_method =='smaller':
        mask[img<threshold]=0
    elif compare_method =='equal':
        maks[img == threshold] =0
    mask[:,remove_columns] = 0
    mask[remove_rows,:] = 0
    if remove_pix!=None:
        for each in remove_pix:
            mask[tuple(each)] = 0
    return mask



