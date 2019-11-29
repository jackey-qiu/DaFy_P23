import sys,os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
from PyQt5 import uic
import random
import numpy as np
import matplotlib.pyplot as plt
try:
    from . import locate_path
except:
    import locate_path
script_path = locate_path.module_path_locator()
DaFy_path = os.path.dirname(os.path.dirname(script_path))
sys.path.append(DaFy_path)
sys.path.append(os.path.join(DaFy_path,'EnginePool'))
sys.path.append(os.path.join(DaFy_path,'FilterPool'))
sys.path.append(os.path.join(DaFy_path,'util'))
import pandas as pd
import time
import matplotlib
matplotlib.use("Qt5Agg")
#from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

class MyMainWindow(QMainWindow):
    def __init__(self, parent = None):
        super(MyMainWindow, self).__init__(parent)
        uic.loadUi(os.path.join(DaFy_path,'scripts','Viewer','data_viewer.ui'),self)
        # self.setupUi(self)
        # plt.style.use('ggplot')
        self.setWindowTitle('XRV and CTR data Viewer')
        self.set_plot_channels()
        matplotlib.rc('xtick', labelsize=10)
        matplotlib.rc('ytick', labelsize=10)
        plt.rcParams.update({'axes.labelsize': 10})
        plt.rc('font',size = 10)
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['axes.grid'] = True
        plt.rcParams['xtick.major.size'] = 6
        plt.rcParams['xtick.major.width'] = 2
        plt.rcParams['xtick.minor.size'] = 4
        plt.rcParams['xtick.minor.width'] = 1
        plt.rcParams['ytick.major.size'] = 6
        plt.rcParams['ytick.major.width'] = 2
        plt.rcParams['ytick.minor.size'] = 4
        plt.rcParams['ytick.minor.width'] = 1
        plt.rcParams['mathtext.default']='regular'
        self.open.clicked.connect(self.load_file)
        self.plot.clicked.connect(self.plot_figure)
        # self.apply.clicked.connect(self.replot_figure)
        self.PushButton_append_scans.clicked.connect(self.append_scans)
        self.PushButton_fold_or_unfold.clicked.connect(self.fold_or_unfold)
        self.checkBox_time_scan.clicked.connect(self.set_plot_channels)
        self.radioButton_ctr.clicked.connect(self.set_plot_channels)
        self.radioButton_xrv.clicked.connect(self.set_plot_channels)
        self.checkBox_mask.clicked.connect(self.append_scans)
        self.data = None
       
    def set_plot_channels(self):
        xrv = self.radioButton_xrv.isChecked()
        time_scan = self.checkBox_time_scan.isChecked()
        if xrv:
            if time_scan:
                self.lineEdit_x.setText('image_no')
                self.lineEdit_y.setText('current,strain_ip,strain_oop,grain_size_ip,grain_size_oop')
            else:
                self.lineEdit_x.setText('potential')
                self.lineEdit_y.setText('current,strain_ip,strain_oop,grain_size_ip,grain_size_oop') 
        else:
            if time_scan:
                self.lineEdit_x.setText('image_no')
                self.lineEdit_y.setText('peak_intensity,potential')
            else:
                self.lineEdit_x.setText('L')
                self.lineEdit_y.setText('peak_intensity') 

    def load_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","Data Files (*.xlsx);;All Files (*.csv)", options=options)
        if fileName:
            self.lineEdit_data_file.setText(fileName)
            self.data = pd.read_excel(fileName)
        col_labels = 'col_labels\n'+str(list(self.data.columns))+'\n'
        scans = list(set(list(self.data['scan_no'])))
        self.scans_all = scans
        scans.sort()
        scan_numbers = 'scan_nos\n'+str(scans)+'\n'
        # print(list(self.data[self.data['scan_no']==scans[0]]['phs'])[0])
        if self.radioButton_xrv.isChecked():
            self.phs_all = [list(self.data[self.data['scan_no']==scan]['phs'])[0] for scan in scans]
            phs = 'pHs\n'+str(self.phs_all)+'\n'
        else:
            phs = ''
        
        self.textEdit_summary_data.setText('\n'.join([col_labels,scan_numbers,phs]))

    #to fold or unfold the config file editor
    def fold_or_unfold(self):
        text = self.PushButton_fold_or_unfold.text()
        if text == "<":
            self.frame.setVisible(False)
            self.PushButton_fold_or_unfold.setText(">")
        elif text == ">":
            self.frame.setVisible(True)
            self.PushButton_fold_or_unfold.setText("<")

    def plot_figure(self):
        if self.radioButton_xrv.isChecked():
            self.plot_figure_xrv()
        else:
            self.plot_figure_ctr()


    def plot_figure_xrv(self):
        self.mplwidget.fig.clear()
        plot_dim = [len(self.plot_labels_y), len(self.scans)]
        for scan in self.scans:
            setattr(self,'plot_axis_scan{}'.format(scan),[])
            j = self.scans.index(scan) + 1
            for i in range(plot_dim[0]):
                getattr(self,'plot_axis_scan{}'.format(scan)).append(self.mplwidget.canvas.figure.add_subplot(plot_dim[0], plot_dim[1],j+plot_dim[1]*i))
        y_max_values,y_min_values = [-100000000]*len(self.plot_labels_y),[100000000]*len(self.plot_labels_y)
        for scan in self.scans:
            for each in self.plot_labels_y:
                i = self.plot_labels_y.index(each)
                try:
                    fmt = self.lineEdit_fmt.text().rsplit(',')[self.scans.index(scan)]
                except:
                    fmt = 'b-'
                getattr(self,'plot_axis_scan{}'.format(scan))[i].plot(self.data_to_plot[scan][self.plot_label_x],self.data_to_plot[scan][self.plot_labels_y[i]],fmt,markersize = 3)
                if i==0:
                    getattr(self,'plot_axis_scan{}'.format(scan))[i].set_title(r'pH {}_scan{}'.format(self.phs[self.scans.index(scan)],scan),fontsize=10)
                temp_max, temp_min = max(list(self.data_to_plot[scan][self.plot_labels_y[i]])),min(list(self.data_to_plot[scan][self.plot_labels_y[i]]))
                if y_max_values[i]<temp_max:
                    y_max_values[i] = temp_max
                if y_min_values[i]>temp_min:
                    y_min_values[i] = temp_min
                if i!=(len(self.plot_labels_y)-1):
                    ax = getattr(self,'plot_axis_scan{}'.format(scan))[i]
                    ax.set_xticklabels([])
                else:
                    ax = getattr(self,'plot_axis_scan{}'.format(scan))[i]
                    x_label = [r'Image_no','E / V$_{RHE}$'][self.plot_label_x=='potential']
                    ax.set_xlabel(x_label)
                if scan!=self.scans[0]:
                    getattr(self,'plot_axis_scan{}'.format(scan))[i].set_yticklabels([])
                else:
                    y_label_map = {'potential':'E / V$_{RHE}$',
                                   'current':r'j / mAcm$^{-2}$',
                                   'strain_ip':r'$\Delta\varepsilon_\parallel$  (%)',
                                   'strain_oop':r'$\Delta\varepsilon_\perp$  (%)',
                                   'grain_size_oop':r'$\Delta d_\perp$ / nm',
                                   'grain_size_ip':r'$\Delta d_\parallel$ / nm',
                                   'peak_intensity':r'Intensity / a.u.'}
                    if each in y_label_map:
                        getattr(self,'plot_axis_scan{}'.format(scan))[i].set_ylabel(y_label_map[each])
                    else:
                        pass
        for scan in self.scans:
            for each in self.plot_labels_y:
                i = self.plot_labels_y.index(each)
                getattr(self,'plot_axis_scan{}'.format(scan))[i].set_ylim(y_min_values[i],y_max_values[i])
        # self.mplwidget.fig.tight_layout()
        self.mplwidget.fig.subplots_adjust(wspace=0.04,hspace=0.04)
        self.mplwidget.canvas.draw()

    def plot_figure_ctr(self):
        self.mplwidget.fig.clear()
        col_num=2#two columns only
        plot_dim = [int(len(self.scan_numbers_all.text().rsplit('+'))/col_num)+int(len(self.scan_numbers_all.text().rsplit('+'))%col_num != 0), col_num]
        for i in range(len(self.scan_numbers_all.text().rsplit('+'))):
            setattr(self,'plot_axis_plot_set{}'.format(i+1),self.mplwidget.canvas.figure.add_subplot(plot_dim[0], plot_dim[1],i+1))
            each = self.scan_numbers_all.text().rsplit('+')[i]
            each = each[1:-1]#remove []
            scan_list_temp = each.rsplit(';')
            for each_set in scan_list_temp:
                j = scan_list_temp.index(each_set)
                sub_scan_list = each_set.rsplit(',')
                scans_temp= [int(i) for i in sub_scan_list]
                for scan in scans_temp:
                    fmt = self.lineEdit_fmt.text().rsplit('+')[i].rsplit(';')[j]
                    if scan == scans_temp[0]:
                        label = self.lineEdit_labels.text().rsplit('+')[i].rsplit(';')[j]
                    else:
                        label = None
                    #remove [ or ] in the fmt and label
                    if ('[' in fmt) and (']' in fmt):
                        fmt = fmt[1:-1]
                    elif '[' in fmt:
                        fmt = fmt[1:]
                    elif ']' in fmt:
                        fmt = fmt[0:-1]
                    else:
                        pass
                    if label != None:
                        if ('[' in label) and (']' in label):
                            label = label[1:-1]
                        elif '[' in label:
                            label = label[1:]
                        elif ']' in label:
                            label = label[0:-1]
                        else:
                            pass
                    if self.checkBox_time_scan.isChecked():
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).plot(self.data_to_plot[scan][self.plot_label_x],self.data_to_plot[scan][self.plot_labels_y[0]],fmt,label =label)
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).set_xlabel(self.plot_label_x)
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).set_ylabel('Intensity')
                        pot_ax = getattr(self,'plot_axis_plot_set{}'.format(i+1)).twinx()
                        pot_ax.plot(self.data_to_plot[scan][self.plot_label_x],self.data_to_plot[scan][self.plot_labels_y[1]],'b-',label = None)
                        pot_ax.set_ylabel(self.plot_labels_y[1],color = 'b')
                        pot_ax.tick_params(axis = 'y', labelcolor = 'b')
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).legend()
                    else:
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).plot(self.data_to_plot[scan][self.plot_label_x],self.data_to_plot[scan][self.plot_labels_y[0]],fmt,label =label)
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).set_ylabel('Intensity')
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).set_xlabel('L')
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).set_title('{}{}L'.format(int(round(list(self.data[self.data['scan_no']==scan]['H'])[0],0)),int(round(list(self.data[self.data['scan_no']==scan]['K'])[0],0))))
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).set_yscale('log')
                        getattr(self,'plot_axis_plot_set{}'.format(i+1)).legend()
        self.mplwidget.fig.tight_layout()
        self.mplwidget.canvas.draw()

    def prepare_data_to_plot_xrv(self,plot_label_list, scan_number):
        scans_ = [int(each) for each in self.lineEdit_scan_numbers.text().rsplit(',')]
        img_ranges_ = [[int(each_) for each_ in each.rsplit('-')] for each in self.lineEdit_point_ranges.text().rsplit(',')]
        
        if scan_number in scans_:
            which = scans_.index(scan_number)
            l , r = img_ranges_[which]
        else:
            l, r = 0, 100000000
        if hasattr(self,'data_to_plot'):
            self.data_to_plot[scan_number] = {}
        else:
            self.data_to_plot = {}
            self.data_to_plot[scan_number] = {}
        if self.checkBox_mask.isChecked():
            condition = (self.data['mask_cv_xrd'] == True)&(self.data['scan_no'] == scan_number)
        else:
            condition = self.data['scan_no'] == scan_number
        for each in plot_label_list:
            if each=='potential':#RHE potential
                self.data_to_plot[scan_number][each] = 0.205+np.array(self.data[condition][each])[l:r]+0.059*np.array(self.data[self.data['scan_no'] == scan_number]['phs'])[0]
            elif each=='current':#RHE potential
                self.data_to_plot[scan_number][each] = -np.array(self.data[condition][each])[l:r]
            else:
                if each in ['peak_intensity','peak_intensity_error','strain_ip','strain_oop','grain_size_ip','grain_size_oop']:
                    temp_data = np.array(self.data[condition][each])[l:r]
                    self.data_to_plot[scan_number][each] = list(temp_data-max(temp_data))
                else:
                    self.data_to_plot[scan_number][each] = list(self.data[condition][each])[l:r]

    def prepare_data_to_plot_ctr(self,plot_label_list, scan_number):
        scans_ = [int(each) for each in self.lineEdit_scan_numbers.text().rsplit(',')]
        img_ranges_ = [[int(each_) for each_ in each.rsplit('-')] for each in self.lineEdit_point_ranges.text().rsplit(',')]
        
        if scan_number in scans_:
            which = scans_.index(scan_number)
            l , r = img_ranges_[which]
        else:
            l, r = 0, 100000000
        if hasattr(self,'data_to_plot'):
            self.data_to_plot[scan_number] = {}
        else:
            self.data_to_plot = {}
            self.data_to_plot[scan_number] = {}
        if self.checkBox_mask.isChecked():
            condition = (self.data['mask_ctr'] == True)&(self.data['scan_no'] == scan_number)
        else:
            condition = self.data['scan_no'] == scan_number
        for each in plot_label_list:
            if each=='current':#RHE potential
                self.data_to_plot[scan_number][each] = -np.array(self.data[condition][each])[l:r]
                #self.data_to_plot[scan_number][each] = 0.205+np.array(self.data[condition][each])[l:r]+0.059*np.array(self.data[self.data['scan_no'] == scan_number]['phs'])[0]               
            else:
                self.data_to_plot[scan_number][each] = np.array(self.data[condition][each])[l:r]

    def append_scans(self):
        if self.radioButton_xrv.isChecked():
            self.append_scans_xrv()
        else:
            self.append_scans_ctr()

    def append_scans_xrv(self):
        text = self.scan_numbers_append.text()
        text_original = self.scan_numbers_all.text()
        if text_original!='':
            text_new = ','.join([text_original, text])
        else:
            text_new = text
        scans = list(set([int(each) for each in text_new.rstrip().rsplit(',')]))
        scans.sort()
        self.scan_numbers_all.setText(','.join([str(scan) for scan in scans]))
        assert (self.lineEdit_x.text()!='' and self.lineEdit_y.text()!=''), 'No channels for plotting have been selected!'
        assert self.scan_numbers_all.text()!='', 'No scans have been selected!'
        plot_labels = self.lineEdit_x.text() + ',' + self.lineEdit_y.text()
        plot_labels = plot_labels.rstrip().rsplit(',')
        for scan in scans:
            self.prepare_data_to_plot_xrv(plot_labels,scan)
            print('Prepare data for scan {} now!'.format(scan))
        self.scans = scans
        self.phs = [self.phs_all[self.scans_all.index(each)] for each in scans]
        self.plot_label_x = self.lineEdit_x.text()
        self.plot_labels_y = self.lineEdit_y.text().rstrip().rsplit(',')

    #the text looks like: [1,2,3;4,5,6]+[7,8,9;10,11,12]+...
    def append_scans_ctr(self):
        text = self.scan_numbers_append.text()
        text_original = self.scan_numbers_all.text()
        if (text_original!='') and (text not in text_original):
            text_new = '+'.join([text_original, text])
        else:
            text_new = text
        scans = []
        for each in text_new.rstrip().rsplit('+'):
            each = each[1:-1]#remove []
            scan_list_temp = each.rsplit(';')
            for each_set in scan_list_temp:
                sub_scan_list = each_set.rsplit(',')
                scans+= [int(i) for i in sub_scan_list]

        # scans = list(set([int(each) for each in text_new.rstrip().rsplit('+')]))
        #scans.sort()
        self.scan_numbers_all.setText(text_new)
        assert (self.lineEdit_x.text()!='' and self.lineEdit_y.text()!=''), 'No channels for plotting have been selected!'
        assert self.scan_numbers_all.text()!='', 'No scans have been selected!'
        plot_labels = self.lineEdit_x.text() + ',' + self.lineEdit_y.text()
        plot_labels = plot_labels.rstrip().rsplit(',')
        for scan in scans:
            self.prepare_data_to_plot_ctr(plot_labels,scan)
            print('Prepare data for scan {} now!'.format(scan))
        self.scans = scans
        #self.phs = [self.phs_all[self.scans_all.index(each)] for each in scans]
        self.plot_label_x = self.lineEdit_x.text()
        self.plot_labels_y = self.lineEdit_y.text().rstrip().rsplit(',')


if __name__ == "__main__":
    app = QApplication(sys.argv)
    myWin = MyMainWindow()
    myWin.show()
    sys.exit(app.exec_())