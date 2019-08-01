import sys
import matplotlib
# matplotlib.use("tkAgg")
from numpy import dtype
sys.path.append('./XRD_tools/')
import numpy as np
from scipy.interpolate import griddata
import scipy.optimize as opt
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
# from pyspec import spec
import subprocess
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.animation import FFMpegWriter
import matplotlib.patches as patches

def draw_lines_on_image(ax_handle,x_y_grid,variable_list,direction = 'horizontal',\
                        color='gray',line_style='-',marker = None,\
                        xlabel=r'$q_\parallel$ / $\AA^{-1}$',\
                        ylabel='$q_\perp$ / $\AA^{-1}$',\
                        fontsize=20,
                        debug = False):
    line_ax_container = []
    x_couples, y_couples = [], []
    if direction == 'horizontal':
        # print(x_y_grid[0,:][0])
        x_couples = [x_y_grid[0,:][[0,-1]]]*len(variable_list)
        y_couples = [[each,each] for each in variable_list]
    elif direction == 'vertical':
        y_couples = [x_y_grid[:,0][[0,-1]]]*len(variable_list)
        x_couples = [[each,each] for each in variable_list]
    for i in range(len(x_couples)):
        temp_line_ax = ax_handle.plot(x_couples[i],y_couples[i],line_style, color = color, marker = marker)
        # line_ax_container.append(temp_line_ax)
    if debug:
        print(x_couples,y_couples)
    ax_handle.set_xlabel(xlabel,fontsize=fontsize)
    ax_handle.set_ylabel(ylabel,fontsize=fontsize)
    return ax_handle

def show_all_plots(fig,grid_intensity,grid_q_ip,grid_q_oop, vmin, vmax, cmap, is_zap_scan, fit_data, model, fit_results, processed_data_container={},bkg_int=None, cut_offset=None,peak_center=None,title = None):
    fig.clear()
    if processed_data_container == {}:
        ax_im=fig.add_subplot(131)
        ax_ip=fig.add_subplot(132)
        ax_oop=fig.add_subplot(133)
    else:
        ax_im=fig.add_subplot(341)
        ax_ip=fig.add_subplot(342)
        ax_oop=fig.add_subplot(343)

        ax_cv = fig.add_subplot(345)
        ax_v = fig.add_subplot(349)
        ax_strain_ip=fig.add_subplot(346)
        ax_strain_oop=fig.add_subplot(347)
        ax_width_ip=fig.add_subplot(3,4,10)
        ax_width_oop=fig.add_subplot(3,4,11)

        ax_bkg_profile = fig.add_subplot(344)
        ax_ctr_profile_pot = fig.add_subplot(348)
        ax_ctr_profile_time = fig.add_subplot(3,4,12)

    ax_im.set_title(title)
    # ax_im.pcolormesh(grid_q_ip, grid_q_oop, grid_intensity, vmin = vmin, vmax = vmax, cmap = cmap)
    ax_im.pcolormesh(grid_q_ip, grid_q_oop, grid_intensity, vmin = vmin, vmax = grid_intensity.max()*.5, cmap = cmap)
    ax_im = draw_lines_on_image(ax_im, grid_q_ip, variable_list = [fit_data['ver']['x'][0].min(),fit_data['ver']['x'][0].max()],direction = 'horizontal',color = 'gray')
    ax_im = draw_lines_on_image(ax_im, grid_q_ip, variable_list = [fit_data['ver']['x'][1].min(),fit_data['ver']['x'][1].max()],direction = 'horizontal',color = 'red')
    ax_im = draw_lines_on_image(ax_im, grid_q_oop, variable_list = [fit_data['hor']['x'][0].min(),fit_data['hor']['x'][0].max()],direction = 'vertical',color = 'gray')
    ax_im = draw_lines_on_image(ax_im, grid_q_oop, variable_list = [fit_data['hor']['x'][1].min(),fit_data['hor']['x'][1].max()],direction = 'vertical',color = 'red')
    # print([fit_data['ver']['x'][0].min(),fit_data['ver']['x'][0].max()])
    cut_values_oop=[peak_center[0]-cut_offset['hor'][-1],peak_center[0]+cut_offset['hor'][-1]]
    cut_values_ip = [peak_center[1]-cut_offset['ver'][-1],peak_center[1]+cut_offset['ver'][-1]]
    # print(cut_values_oop, cut_values_ip)
    # print('sensor',peak_center,cut_values_ip, cut_values_oop)
    ax_im = draw_lines_on_image(ax_im, grid_q_ip, variable_list = [grid_q_oop[each,0] for each in cut_values_oop], direction = 'horizontal',color = 'm')
    ax_im = draw_lines_on_image(ax_im, grid_q_oop, variable_list = [grid_q_ip[0,each] for each in cut_values_ip], direction = 'vertical',color = 'm')

    x_span =  abs(grid_q_ip[0,bkg_int.x_min+bkg_int.x_span]-grid_q_ip[0,bkg_int.x_min])
    y_span =  abs(grid_q_oop[-bkg_int.y_min-bkg_int.y_span,0]-grid_q_oop[-bkg_int.y_min,0])
    rect = patches.Rectangle((grid_q_ip[0,bkg_int.x_min],grid_q_oop[-bkg_int.y_min,0]),x_span, y_span,linewidth=1,edgecolor='g',ls='-',facecolor='none')
    ax_im.add_patch(rect)

    ax_ip.plot(fit_data['hor']['x'][-1],fit_data['hor']['y'][-1])
    ax_ip.plot(fit_data['hor']['x'][0],model(fit_data['hor']['x'][0],*fit_results['hor'][0]))
    q_boundary_ip = [fit_data['hor']['x'][0].min(), fit_data['hor']['x'][0].max(),fit_data['hor']['x'][-1].min(), fit_data['hor']['x'][-1].max()]
    intensity_boundary_ip = [fit_data['hor']['y'][-1].min(), fit_data['hor']['y'][-1].max()]
    for i in range(len(q_boundary_ip)):
        if i in [0,1]:
            color = 'gray'
        else:
            color = 'red'
        ax_ip.plot([q_boundary_ip[i],q_boundary_ip[i]],intensity_boundary_ip,color = color)
    ax_ip.set_ylim((0,fit_data['hor']['y'][-1].max()*1.2))

    ax_ip.set_xlabel(r'$q_\parallel$ / $\AA^{-1}$', fontsize=20)
    ax_ip.set_ylabel(r'Intensity / a.u.', fontsize=20)
    ax_oop.plot(fit_data['ver']['x'][-1],fit_data['ver']['y'][-1])
    ax_oop.plot(fit_data['ver']['x'][0],model(fit_data['ver']['x'][0],*fit_results['ver'][0]))
    q_boundary_oop = [fit_data['ver']['x'][0].min(), fit_data['ver']['x'][0].max(),fit_data['ver']['x'][-1].min(), fit_data['ver']['x'][-1].max()]
    intensity_boundary_oop = [fit_data['ver']['y'][-1].min(), fit_data['ver']['y'][-1].max()]
    for i in range(len(q_boundary_oop)):
        if i in [0,1]:
            color = 'gray'
        else:
            color = 'red'
        ax_oop.plot([q_boundary_oop[i],q_boundary_oop[i]],intensity_boundary_oop,color = color)
    ax_oop.set_ylim((0,fit_data['ver']['y'][-1].max()*1.2))
    ax_oop.set_xlabel(r'$q_\perp$ / $\AA^{-1}$', fontsize=20)
    ax_oop.set_ylabel(r'Intensity / a.u.', fontsize=20)

    ax_bkg_profile.plot(bkg_int.fit_data['x'],bkg_int.fit_data['y_total'],color = 'blue', label ='data')
    ax_bkg_profile.plot(bkg_int.fit_data['x'],bkg_int.fit_data['y_bkg'],color = 'red', label ='background')
    ax_bkg_profile.plot(bkg_int.fit_data['x'],bkg_int.fit_data['y_total']-bkg_int.fit_data['y_bkg'],color = 'm', label ='data-background')
    ax_bkg_profile.plot(bkg_int.fit_data['x'],[0]*len(bkg_int.fit_data['x']))

    ax_bkg_profile.set_xlabel('pixel',fontsize=20)
    ax_bkg_profile.set_ylabel('I',fontsize=20)
    ax_bkg_profile.set_title('bkg_profile')

    ax_ctr_profile_pot.set_xlabel('E(V)',fontsize=20)
    ax_ctr_profile_pot.set_ylabel('I',fontsize=20)
    ax_ctr_profile_pot.set_title('ctr_profile_L')

    ax_ctr_profile_time.set_xlabel('t',fontsize=20)
    ax_ctr_profile_time.set_ylabel('I',fontsize=20)
    ax_ctr_profile_time.set_title('CTR_profile_t')
    



    index_ctr = np.where(np.array(processed_data_container['mask_ctr'])==1)
    if len(index_ctr[0])==0:
        index_ctr = [0]
    ax_ctr_profile_pot.errorbar(np.array(processed_data_container['potential'])[index_ctr], np.array(processed_data_container['peak_intensity'])[index_ctr],xerr=None,yerr=np.array(processed_data_container['peak_intensity_error'])[index_ctr],fmt='ro:', markersize=4, label='CTR profile')

    ax_ctr_profile_time.errorbar(np.array(range(len(processed_data_container['L'])))[index_ctr], np.array(processed_data_container['peak_intensity'])[index_ctr],xerr=None,yerr=np.array(processed_data_container['peak_intensity_error'])[index_ctr],fmt='ro:', markersize=4, label='CTR profile')

    if is_zap_scan: 
        ax_strain_oop.plot(processed_data_container['oop_strain'],'r-.')
        ax_strain_ip.plot(processed_data_container['ip_strain'],'g-.')
        ax_width_oop.plot(processed_data_container['oop_grain_size'],'r-.')
        ax_width_ip.plot(processed_data_container['ip_grain_size'],'g-.')
        ax_cv.plot(processed_data_container['potential'], processed_data_container['current_density'],'k--.')
        ax_v.plot(processed_data_container['potential'],'m:d')
        ax_strain_oop.plot(len(processed_data_container['potential'])-1,processed_data_container['oop_strain'][-1],'g-o')
        ax_strain_ip.plot(len(processed_data_container['potential'])-1,processed_data_container['ip_strain'][-1],'r-o')
        ax_width_oop.plot(len(processed_data_container['potential'])-1,processed_data_container['oop_grain_size'][-1],'g-o')
        ax_width_ip.plot(len(processed_data_container['potential'])-1,processed_data_container['ip_grain_size'][-1],'r-o')
        ax_cv.plot(processed_data_container['potential'][-1], processed_data_container['current_density'][-1],'r--o')
        ax_v.plot(len(processed_data_container['potential'])-1,processed_data_container['potential'][-1],'k:d')
    else:
        # ax_strain_oop.plot(processed_data_container['potential'],processed_data_container['oop_strain'],'r-.')
        # ax_strain_ip.plot(processed_data_container['potential'],processed_data_container['ip_strain'],'g-.')
        # ax_width_oop.plot(processed_data_container['potential'],processed_data_container['oop_grain_size'],'r-.')
        # ax_width_ip.plot(processed_data_container['potential'],processed_data_container['ip_grain_size'],'g-.')
        # ax_cv.plot(processed_data_container['potential'], processed_data_container['current_density'],'k--.')
        # ax_v.plot(processed_data_container['potential'],'m:d')
        # ax_strain_oop.plot(processed_data_container['potential'][-1],processed_data_container['oop_strain'][-1],'g-o')
        # ax_strain_ip.plot(processed_data_container['potential'][-1],processed_data_container['ip_strain'][-1],'r-o')
        # ax_width_oop.plot(processed_data_container['potential'][-1],processed_data_container['oop_grain_size'][-1],'g-o')
        # ax_width_ip.plot(processed_data_container['potential'][-1],processed_data_container['ip_grain_size'][-1],'r-o')
        # ax_cv.plot(processed_data_container['potential'][-1], processed_data_container['current_density'][-1],'r--o')
        index_cv = np.where(np.array(processed_data_container['mask_cv_xrd'])==1)
        if len(index_cv[0])==0:
            index_cv = [0]
        ax_strain_oop.plot(np.array(processed_data_container['oop_strain'])[index_cv],'r-.')
        ax_strain_ip.plot(np.array(processed_data_container['ip_strain'])[index_cv],'g-.')
        ax_width_oop.plot(np.array(processed_data_container['oop_grain_size'])[index_cv],'r-.')
        ax_width_ip.plot(np.array(processed_data_container['ip_grain_size'])[index_cv],'g-.')
        ax_cv.plot(processed_data_container['potential'], processed_data_container['current_density'],'k--.')
        ax_v.plot(processed_data_container['potential'],'m:d')
        ax_strain_oop.plot(processed_data_container['oop_strain'][-1],'g-o')
        ax_strain_ip.plot(processed_data_container['ip_strain'][-1],'r-o')
        ax_width_oop.plot(processed_data_container['oop_grain_size'][-1],'g-o')
        ax_width_ip.plot(processed_data_container['ip_grain_size'][-1],'r-o')
        ax_cv.plot(processed_data_container['potential'][-1], processed_data_container['current_density'][-1],'r--o')

    ax_v.set_xlabel('Time')
    ax_v.set_title('potential')
    ax_strain_oop.set_title('out-of-plane strain')
    ax_strain_ip.set_title('inplane strain')
    ax_width_oop.set_title('out-of-plane width')
    ax_width_ip.set_title('inplane width')
    ax_cv.set_title('CV')
    if is_zap_scan:
        ax_cv.set_xlabel('E(V)')
        ax_strain_ip.set_xlabel('Time')
        ax_strain_oop.set_xlabel('Time')
        ax_width_ip.set_xlabel('Time')
        ax_width_oop.set_xlabel('Time')
    else:
        ax_cv.set_xlabel('E(V)')
        ax_strain_ip.set_xlabel('E(V)')
        ax_strain_oop.set_xlabel('E(V)')
        ax_width_ip.set_xlabel('E(V)')
        ax_width_oop.set_xlabel('E(V)')

    return fig

def show_3_plots(fig,grid_intensity,grid_q_ip,grid_q_oop, vmin, vmax, cmap, is_zap_scan, fit_data, model, fit_results, processed_data_container={}, cut_offset=None,peak_center=None,title = None):
    fig.clear()
    if processed_data_container == {}:
        ax_im=fig.add_subplot(131)
        ax_ip=fig.add_subplot(132)
        ax_oop=fig.add_subplot(133)
    else:
        ax_im=fig.add_subplot(331)
        ax_ip=fig.add_subplot(332)
        ax_oop=fig.add_subplot(333)

        ax_cv = fig.add_subplot(334)
        ax_v = fig.add_subplot(337)
        ax_strain_ip=fig.add_subplot(335)
        ax_strain_oop=fig.add_subplot(336)
        ax_width_ip=fig.add_subplot(338)
        ax_width_oop=fig.add_subplot(339)
    ax_im.set_title(title)
    ax_im.pcolormesh(grid_q_ip, grid_q_oop, grid_intensity, vmin = vmin, vmax = vmax, cmap = cmap)
    ax_im = draw_lines_on_image(ax_im, grid_q_ip, variable_list = [fit_data['ver']['x'][0].min(),fit_data['ver']['x'][0].max()],direction = 'horizontal',color = 'gray')
    ax_im = draw_lines_on_image(ax_im, grid_q_ip, variable_list = [fit_data['ver']['x'][1].min(),fit_data['ver']['x'][1].max()],direction = 'horizontal',color = 'red')
    ax_im = draw_lines_on_image(ax_im, grid_q_oop, variable_list = [fit_data['hor']['x'][0].min(),fit_data['hor']['x'][0].max()],direction = 'vertical',color = 'gray')
    ax_im = draw_lines_on_image(ax_im, grid_q_oop, variable_list = [fit_data['hor']['x'][1].min(),fit_data['hor']['x'][1].max()],direction = 'vertical',color = 'red')
    # print([fit_data['ver']['x'][0].min(),fit_data['ver']['x'][0].max()])
    cut_values_oop=[peak_center[0]-cut_offset['hor'][-1],peak_center[0]+cut_offset['hor'][-1]]
    cut_values_ip = [peak_center[1]-cut_offset['ver'][-1],peak_center[1]+cut_offset['ver'][-1]]
    # print(cut_values_oop, cut_values_ip)
    # print('sensor',peak_center,cut_values_ip, cut_values_oop)
    ax_im = draw_lines_on_image(ax_im, grid_q_ip, variable_list = [grid_q_oop[each,0] for each in cut_values_oop], direction = 'horizontal',color = 'm')
    ax_im = draw_lines_on_image(ax_im, grid_q_oop, variable_list = [grid_q_ip[0,each] for each in cut_values_ip], direction = 'vertical',color = 'm')

    ax_ip.plot(fit_data['hor']['x'][-1],fit_data['hor']['y'][-1])
    ax_ip.plot(fit_data['hor']['x'][0],model(fit_data['hor']['x'][0],*fit_results['hor'][0]))
    q_boundary_ip = [fit_data['hor']['x'][0].min(), fit_data['hor']['x'][0].max(),fit_data['hor']['x'][-1].min(), fit_data['hor']['x'][-1].max()]
    intensity_boundary_ip = [fit_data['hor']['y'][-1].min(), fit_data['hor']['y'][-1].max()]
    for i in range(len(q_boundary_ip)):
        if i in [0,1]:
            color = 'gray'
        else:
            color = 'red'
        ax_ip.plot([q_boundary_ip[i],q_boundary_ip[i]],intensity_boundary_ip,color = color)
    ax_ip.set_ylim((0,fit_data['hor']['y'][-1].max()*1.2))

    ax_ip.set_xlabel(r'$q_\parallel$ / $\AA^{-1}$', fontsize=20)
    ax_ip.set_ylabel(r'Intensity / a.u.', fontsize=20)
    ax_oop.plot(fit_data['ver']['x'][-1],fit_data['ver']['y'][-1])
    ax_oop.plot(fit_data['ver']['x'][0],model(fit_data['ver']['x'][0],*fit_results['ver'][0]))
    q_boundary_oop = [fit_data['ver']['x'][0].min(), fit_data['ver']['x'][0].max(),fit_data['ver']['x'][-1].min(), fit_data['ver']['x'][-1].max()]
    intensity_boundary_oop = [fit_data['ver']['y'][-1].min(), fit_data['ver']['y'][-1].max()]
    for i in range(len(q_boundary_oop)):
        if i in [0,1]:
            color = 'gray'
        else:
            color = 'red'
        ax_oop.plot([q_boundary_oop[i],q_boundary_oop[i]],intensity_boundary_oop,color = color)
    ax_oop.set_ylim((0,fit_data['ver']['y'][-1].max()*1.2))
    ax_oop.set_xlabel(r'$q_\perp$ / $\AA^{-1}$', fontsize=20)
    ax_oop.set_ylabel(r'Intensity / a.u.', fontsize=20)
    if is_zap_scan: 
        ax_strain_oop.plot(processed_data_container['oop_strain'],'r-.')
        ax_strain_ip.plot(processed_data_container['ip_strain'],'g-.')
        ax_width_oop.plot(processed_data_container['oop_grain_size'],'r-.')
        ax_width_ip.plot(processed_data_container['ip_grain_size'],'g-.')
        ax_cv.plot(processed_data_container['potential'], processed_data_container['current_density'],'k--.')
        ax_v.plot(processed_data_container['potential'],'m:d')
        ax_strain_oop.plot(len(processed_data_container['potential'])-1,processed_data_container['oop_strain'][-1],'g-o')
        ax_strain_ip.plot(len(processed_data_container['potential'])-1,processed_data_container['ip_strain'][-1],'r-o')
        ax_width_oop.plot(len(processed_data_container['potential'])-1,processed_data_container['oop_grain_size'][-1],'g-o')
        ax_width_ip.plot(len(processed_data_container['potential'])-1,processed_data_container['ip_grain_size'][-1],'r-o')
        ax_cv.plot(processed_data_container['potential'][-1], processed_data_container['current_density'][-1],'r--o')
        ax_v.plot(len(processed_data_container['potential'])-1,processed_data_container['potential'][-1],'k:d')
    else:
        # ax_strain_oop.plot(processed_data_container['potential'],processed_data_container['oop_strain'],'r-.')
        # ax_strain_ip.plot(processed_data_container['potential'],processed_data_container['ip_strain'],'g-.')
        # ax_width_oop.plot(processed_data_container['potential'],processed_data_container['oop_grain_size'],'r-.')
        # ax_width_ip.plot(processed_data_container['potential'],processed_data_container['ip_grain_size'],'g-.')
        # ax_cv.plot(processed_data_container['potential'], processed_data_container['current_density'],'k--.')
        # ax_v.plot(processed_data_container['potential'],'m:d')
        # ax_strain_oop.plot(processed_data_container['potential'][-1],processed_data_container['oop_strain'][-1],'g-o')
        # ax_strain_ip.plot(processed_data_container['potential'][-1],processed_data_container['ip_strain'][-1],'r-o')
        # ax_width_oop.plot(processed_data_container['potential'][-1],processed_data_container['oop_grain_size'][-1],'g-o')
        # ax_width_ip.plot(processed_data_container['potential'][-1],processed_data_container['ip_grain_size'][-1],'r-o')
        # ax_cv.plot(processed_data_container['potential'][-1], processed_data_container['current_density'][-1],'r--o')
        ax_strain_oop.plot(processed_data_container['oop_strain'],'r-.')
        ax_strain_ip.plot(processed_data_container['ip_strain'],'g-.')
        ax_width_oop.plot(processed_data_container['oop_grain_size'],'r-.')
        ax_width_ip.plot(processed_data_container['ip_grain_size'],'g-.')
        ax_cv.plot(processed_data_container['potential'], processed_data_container['current_density'],'k--.')
        ax_v.plot(processed_data_container['potential'],'m:d')
        ax_strain_oop.plot(processed_data_container['oop_strain'][-1],'g-o')
        ax_strain_ip.plot(processed_data_container['ip_strain'][-1],'r-o')
        ax_width_oop.plot(processed_data_container['oop_grain_size'][-1],'g-o')
        ax_width_ip.plot(processed_data_container['ip_grain_size'][-1],'r-o')
        ax_cv.plot(processed_data_container['potential'][-1], processed_data_container['current_density'][-1],'r--o')
    ax_v.set_xlabel('Time')
    ax_v.set_title('potential')
    ax_strain_oop.set_title('out-of-plane strain')
    ax_strain_ip.set_title('inplane strain')
    ax_width_oop.set_title('out-of-plane width')
    ax_width_ip.set_title('inplane width')
    ax_cv.set_title('CV')
    if is_zap_scan:
        ax_cv.set_xlabel('E(V)')
        ax_strain_ip.set_xlabel('Time')
        ax_strain_oop.set_xlabel('Time')
        ax_width_ip.set_xlabel('Time')
        ax_width_oop.set_xlabel('Time')
    else:
        ax_cv.set_xlabel('E(V)')
        ax_strain_ip.set_xlabel('E(V)')
        ax_strain_oop.set_xlabel('E(V)')
        ax_width_ip.set_xlabel('E(V)')
        ax_width_oop.set_xlabel('E(V)')

    return fig

def plot_after_fit(data,is_zap_scan):
    if not is_zap_scan:
        plt.figure(111)
        plt.plot(data['potential'],data['ip_strain'],'o-',markersize=5)
        plt.title('Inplane strain(%)')
        plt.xlabel('Potential(v)')
        plt.figure(112)
        plt.plot(data['potential'],data['ip_grain_size'])
        plt.title('Inplane crystallite size (nm)')
        plt.xlabel('Potential(v)')
        plt.figure(222)
        plt.plot(data['potential'],data['oop_strain'],'o-',markersize=5)
        plt.title('Out of plane strain(%)')
        plt.xlabel('Potential(v)')
        plt.figure(223)
        plt.plot(data['potential'],data['oop_grain_size'])
        plt.xlabel('Potential(v)')
        plt.title('Out of plane crystallite size(nm)')
        plt.figure(224)
        plt.plot(data['potential'],data['peak_intensity'])
        plt.xlabel('Potential(v)')
        plt.title('Peak intensity (counts)')
        plt.show()
    else:
        plt.figure(111)
        plt.plot(data['cen_ip'],'o-',markersize=5)
        plt.title('Inplane peak position')
        plt.figure(112)
        plt.plot(data['FWHM_ip'])
        plt.title('Inplane peak width')
        plt.figure(222)
        plt.plot(data['cen_oop'],'o-',markersize=5)
        plt.title('Out of plane peak position')
        plt.figure(223)
        plt.plot(data['FWHM_oop'])
        plt.title('Out of plane peak width')
        plt.show()


def movie_creator(fig_handle, movie_name,fps = 5):
    canvas_width, canvas_height = fig_handle.canvas.get_width_height()
    # Open an ffmpeg process
    outf = movie_name
    cmdstring = ('ffmpeg',
            '-y', '-r', str(fps), # overwrite, 30fps
            '-s', '%dx%d' % (canvas_width, canvas_height), # size of image string
            '-pix_fmt', 'argb', # format
            '-f', 'rawvideo',  '-i', '-', # tell ffmpeg to expect raw video from the pipe
            '-vcodec', 'mpeg4', outf) # output encoding
    p = subprocess.Popen(cmdstring, stdin=subprocess.PIPE)
    agg=fig_handle.canvas.switch_backends(FigureCanvasAgg)
    return p, agg

def update_line(line_ax_handles, atr_value_dic={'set_data':None,'set_cmap':'gnuplot2'}):
    num_ax = len(line_ax_handles)
    updated_lines = []
    for i in range(num_ax):
        line_ax_handle = line_ax_handles[i]
        for key in atr_value_dic.keys():
            getattr(line_ax_handle, key)(atr_value_dic[key])
        updated_lines.append(line_ax_handle)
    return line_ax_handles

