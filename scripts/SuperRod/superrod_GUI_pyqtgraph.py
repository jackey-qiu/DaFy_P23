import sys,os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from PyQt5 import uic
import random
import numpy as np
import pandas as pd
import types
import matplotlib.pyplot as plt
try:
    from . import locate_path
except:
    import locate_path
script_path = locate_path.module_path_locator()
DaFy_path = os.path.dirname(os.path.dirname(script_path))
sys.path.append(DaFy_path)
sys.path.append(os.path.join(DaFy_path,'dump_files'))
sys.path.append(os.path.join(DaFy_path,'EnginePool'))
sys.path.append(os.path.join(DaFy_path,'FilterPool'))
sys.path.append(os.path.join(DaFy_path,'util'))
from fom_funcs import *
import parameters
from VisualizationEnginePool import plot_bkg_fit_gui_pyqtgraph,replot_bkg_profile
import data_superrod as data
import time
import matplotlib
matplotlib.use("Qt5Agg")
import pyqtgraph as pg
import pyqtgraph.exporters
from PyQt5 import QtCore
from PyQt5.QtWidgets import QCheckBox, QRadioButton, QTableWidgetItem, QHeaderView, QAbstractItemView
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QTransform, QFont, QBrush, QColor
from pyqtgraph.Qt import QtGui
import syntax_pars

#from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

class MyMainWindow(QMainWindow):
    def __init__(self, parent = None):
        super(MyMainWindow, self).__init__(parent)
        pg.setConfigOptions(imageAxisOrder='row-major')
        pg.mkQApp()
        uic.loadUi(os.path.join(DaFy_path,'scripts','SuperRod','superrod.ui'),self)
        self.setWindowTitle('Data analysis factory: CTR data modeling')
        self.stop = False
        self.show_checkBox_list = []
        #set fom_func
        self.fom_func = chi2bars_2
        #parameters
        self.parameters = parameters.Parameters()
        #scripts
        self.script = ''
        #script module
        self.script_module = types.ModuleType('genx_script_module')
        #tool bar buttons to operate modeling
        self.actionNew.triggered.connect(self.init_new_model)
        self.actionOpen.triggered.connect(self.open_model)
        self.actionSave.triggered.connect(self.save_model)
        self.actionSimulate.triggered.connect(self.simulate_model)
        self.actionRun.triggered.connect(self.run_model)
        self.actionStop.triggered.connect(self.stop_model)
        #pushbuttons for data handeling
        self.pushButton_load_data.clicked.connect(self.load_data_ctr)
        self.pushButton_append_data.clicked.connect(self.append_data)
        self.pushButton_delete_data.clicked.connect(self.delete_data)
        self.pushButton_save_data.clicked.connect(self.save_data)
        self.pushButton_calculate.clicked.connect(self.calculate)
        #pushbutton for changing plotting style
        self.pushButton_plot_style.clicked.connect(self.change_plot_style)
        #pushbutton to load/save script
        self.pushButton_load_script.clicked.connect(self.load_script)
        self.pushButton_save_script.clicked.connect(self.save_script)
        #pushbutton to load/save parameter file
        self.pushButton_load_table.clicked.connect(self.load_par)
        self.pushButton_save_table.clicked.connect(self.save_par)
        #select dataset in the viewer
        self.comboBox_dataset.activated.connect(self.update_data_view)

        #syntax highlight
        self.plainTextEdit_script.setStyleSheet("""QPlainTextEdit{
	                            font-family:'Consolas'; 
                                font-size:11pt;
	                            color: #ccc; 
	                            background-color: #2b2b2b;}""")
        self.plainTextEdit_script.setTabStopWidth(self.plainTextEdit_script.fontMetrics().width(' ')*4)
        self.data = data.DataList()

        #table view for parameters set to selecting row basis
        self.tableWidget_pars.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.timer_save_data = QtCore.QTimer(self)
        self.timer_save_data.timeout.connect(self.save_model)
        self.setup_plot()
        
    def setup_plot(self):
        self.selected_data_profile = self.widget_data.addPlot()
        self.fom_evolution_profile = self.widget_fom.addPlot()
        self.par_profile = self.widget_pars.addPlot()
        self.fom_scan_profile = self.widget_fom_scan.addPlot()
        #self.widget_edp.setup_view()

    def update_plot_data_view(self):
        plot_data_index = []
        for i in range(len(self.data)):
            if self.tableWidget_data.cellWidget(i,1).isChecked():
                # self.selected_data_profile.plot(self.data[i].x, self.data[i].y, clear = True)
                self.selected_data_profile.plot(self.data[i].x, self.data[i].y,pen={'color': 'y', 'width': 1},  symbolBrush=(255,0,0), symbolSize=5,symbolPen='w', clear = (len(plot_data_index) == 0))
                plot_data_index.append(i)
        self.selected_data_profile.setLogMode(x=False,y=True)
        self.selected_data_profile.autoRange()

    def update_plot_data_view_upon_simulation(self):
        plot_data_index = []
        for i in range(len(self.data)):
            if self.tableWidget_data.cellWidget(i,1).isChecked():
                # self.selected_data_profile.plot(self.data[i].x, self.data[i].y, clear = True)
                self.selected_data_profile.plot(self.data[i].x, self.data[i].y,pen={'color': 'y', 'width': 0},  symbolBrush=(255,0,0), symbolSize=5,symbolPen='w', clear = (len(plot_data_index) == 0))
                self.selected_data_profile.plot(self.data[i].x, self.data[i].y_sim,pen={'color': 'r', 'width': 2},  clear = False)
                
                plot_data_index.append(i)
        self.selected_data_profile.setLogMode(x=False,y=True)
        self.selected_data_profile.autoRange()

    def test_table_widget(self):
        self.tableWidget_data.setRowCount(3)
        self.tableWidget_data.setColumnCount(4)
        self.tableWidget_data.setHorizontalHeaderLabels(['Name','Show','Use','Errors'])
        self.tableWidget_data.setItem(0,0,QTableWidgetItem('(00L)'))
        # self.tableWidget_data.setItem(0,1,QTableWidgetItem('No'))
        self.tableWidget_data.setCellWidget(0,1,QCheckBox())
        self.tableWidget_data.setItem(0,2,QTableWidgetItem('Yes'))
        self.tableWidget_data.setItem(0,3,QTableWidgetItem('No'))

        self.tableWidget_pars.setRowCount(3)
        self.tableWidget_pars.setColumnCount(6)
        self.tableWidget_pars.setHorizontalHeaderLabels(['Parameter','Value','Fit','Min','Max','Error'])
        self.tableWidget_pars.setItem(0,0,QTableWidgetItem('rgh.setBeta'))
        self.tableWidget_pars.setItem(0,1,QTableWidgetItem('0.099726'))
        self.tableWidget_pars.setCellWidget(0,2,QCheckBox())
        self.tableWidget_pars.setItem(0,3,QTableWidgetItem('0'))
        self.tableWidget_pars.setItem(0,4,QTableWidgetItem('1'))
        self.tableWidget_pars.setItem(0,5,QTableWidgetItem('0.001,0.003'))

    def update_plot(self):
        pass

    def init_new_model(self):
        pass

    def open_model(self):
        pass

    def save_model(self):
        pass

    def simulate_model(self):
        self.compile_script()
        (funcs, vals) = self.get_sim_pars()
        # Set the parameter values in the model
        #[func(val) for func,val in zip(funcs, vals)]
        i = 0
        for func, val in zip(funcs,vals):
            try:
                func(val)
            except Exception as e:
                (sfuncs_tmp, vals_tmp) = self.parameters.get_sim_pars()
                raise ParameterError(sfuncs_tmp[i], i, str(e), 1)
            i += 1

        self.evaluate_sim_func()
        self.update_plot_data_view_upon_simulation()
        self.update_structure_view()

    def calc_fom(self, simulated_data,wt=1,wt_list=[]):
        '''calc_fom(self, fomlist) -> fom_raw (list of arrays),
                                      fom_indiv(list of floats),
                                      fom(float)

        Sums up the evaluation of the fom values calculated for each
         data point to form the overall fom function for all data sets.
        '''
        fom_raw = self.fom_func(simulated_data, self.data)
        # Sum up a unique fom for each data set in use
        fom_indiv=[]
        if wt_list==[]:
            fom_indiv = [np.sum(np.abs(fom_set)) for fom_set in fom_raw]
        else:
            fom_indiv = [np.sum(np.abs(fom_set))*scale for (fom_set,scale) in zip(fom_raw,wt_list)]
        fom = np.sum([f for f, d in zip(fom_indiv, self.data) if d.use])
        # Lets extract the number of datapoints as well:
        N = np.sum([len(fom_set) for fom_set, d in zip(fom_raw, self.data) if d.use])
        # And the number of fit parameters
        p = self.parameters.get_len_fit_pars()
        #self.fom_dof = fom/((N-p)*1.0)
        try:
            use_dif = self.fom_func.__div_dof__
        except Exception:
            use_dif = False
        if use_dif:
            if N>p:
                fom = fom/((N-p)*1.0)
            else:
                fom=fom/(N*1.0)

        return fom_raw, fom_indiv, fom*wt
        #return fom_raw,fom_indiv,fom_raw[0][3]

    def evaluate_sim_func(self):
        '''evaluate_sim_func(self) --> None

        Evalute the Simulation function and updates the data simulated data
        as well as the fom of the model. Use this one for calculating data to
        update plots, simulations and such.
        '''
        try:
            simulated_data,wt,wt_list = self.script_module.Sim(self.data)
        except:
            try:
                simulated_data,wt = self.script_module.Sim(self.data)
            except Exception as e:
                outp = StringIO()
                traceback.print_exc(200, outp)
                val = outp.getvalue()
                outp.close()
                raise ModelError(str(val), 1)

        # check so that the Sim function returns anything
        if not simulated_data:
            text = 'The Sim function does not return anything, it should' +\
            ' return a list of the same length as the number of data sets.'
            raise ModelError(text, 1)
        # Check so the number of data sets is correct
        if len(simulated_data) != len(self.data):
            text = 'The number of simulated data sets returned by the Sim function'\
             + ' has to be same as the number of loaded data sets.\n' +\
             'Number of loaded data sets: ' + str(len(self.data)) +\
             '\nNumber of simulated data sets: ' + str(len(simulated_data))
            raise ModelError(text, 1)

        self.data.set_simulated_data(simulated_data)

        try:
            #self.fom = self.fom_func(simulated_data, self.data)
            fom_raw, fom_inidv, fom = self.calc_fom(simulated_data,wt,wt_list)
            self.fom = fom
        except:
            try:
                fom_raw, fom_inidv, fom = self.calc_fom(simulated_data,wt)
                self.fom = fom
            except Exception as e:
                outp = StringIO()
                traceback.print_exc(200, outp)
                val = outp.getvalue()
                outp.close()
                raise FomError(str(val))
        #print len(fom_raw)
        self.data.set_fom_data(fom_raw)

    def get_sim_pars(self):
        ''' get_sim_pars(self) --> (funcs, values)

        Returns the parameters used with simulations. i.e. the function to
        set the parameters, the guess value (values). Used for simulation,
        for fitting see get_fit_pars(self).s
        '''
        (sfuncs, vals) = self.parameters.get_sim_pars()
        # Compile the strings to create the functions..
        funcs = []
        for func in sfuncs:
            # funcs.append(self.create_fit_func(func))
            try:
                funcs.append(self.create_fit_func(func))
            except Exception as e:
                raise ParameterError(func, len(funcs), str(e),0)
        return (funcs, vals)

    def create_fit_func(self, str):
        '''create_fit_func(self, str) --> function

        Creates a function from the string expression in string.
        If the string is a function in the model this function will be
        returned if string represents anything else a function that sets that
        object will be returned.
        '''
        object = self.eval_in_model(str)
        #print type(object)
        # Is it a function or a method!
        name = type(object).__name__
        if name == 'instancemethod' or name == 'function' or name == 'method':
            return object
        # Nope lets make a function of it
        else:
            print(name)
            #print 'def __tempfunc__(val):\n\t%s = val'%str
            #The function must be created in the module in order to acess
            # the different variables
            exec('def __tempfunc__(val):\n\t%s = val'%str\
                in self.script_module.__dict__)

            #print self.script_module.__tempfunc__
            return self.script_module.__tempfunc__

    def eval_in_model(self, codestring):
        '''
        Excecute the code in codestring in the namespace of
        model module
        '''
        #exec codestring in self.script_module.__dict__
        # result = eval(codestring, self.script_module.__dict__)
        result = eval(codestring, self.script_module.__dict__)
        # print('Sucessfully evaluted: ', codestring)
        return result

    def run_model(self):
        pass

    def stop_model(self):
        pass

    def load_data(self, loader = 'ctr'):
        exec('self.load_data_{}()'.format(loader))

    def load_data_ctr(self):
        #8 columns in total
        #X, H, K, Y, I, eI, LB, dL 
        #for CTR data, X column is L column, Y column all 0
        #for RAXR data, X column is energy column, Y column is L column
        self.data = data.DataList()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","csv Files (*.csv);;data Files (*.dat);txt Files (*.txt)", options=options)
        if fileName:
            with open(fileName,'r') as f:
                data_loaded = np.loadtxt(f,comments = '#',delimiter=None)
                data_loaded_pd = pd.DataFrame(data_loaded, columns = ['X','h','k','Y','I','eI','LB','dL'])
                data_loaded_pd['h'] = data_loaded_pd['h'].apply(lambda x:int(np.round(x)))
                data_loaded_pd['k'] = data_loaded_pd['k'].apply(lambda x:int(np.round(x)))
                data_loaded_pd.sort_values(by = ['h','k'], inplace = True)
                # print(data_loaded_pd)
                hk_unique = list(set(zip(list(data_loaded_pd['h']), list(data_loaded_pd['k']))))
                hk_unique.sort()
                h_unique = [each[0] for each in hk_unique]
                k_unique = [each[1] for each in hk_unique]
                for i in range(len(h_unique)):
                    h_temp, k_temp = h_unique[i], k_unique[i]
                    name = 'Data-{}{}L'.format(h_temp, k_temp)
                    self.data.add_new(name = name)
                    self.data.items[-1].x = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['X'].to_numpy()
                    self.data.items[-1].y = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['I'].to_numpy()
                    self.data.items[-1].error = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['eI'].to_numpy()
                    self.data.items[-1].x_raw = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['X'].to_numpy()
                    self.data.items[-1].y_raw = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['I'].to_numpy()
                    self.data.items[-1].error_raw = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['eI'].to_numpy()
                    self.data.items[-1].set_extra_data(name = 'h', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['h'].to_numpy())
                    self.data.items[-1].set_extra_data(name = 'k', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['k'].to_numpy())
                    self.data.items[-1].set_extra_data(name = 'Y', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['Y'].to_numpy())
                    self.data.items[-1].set_extra_data(name = 'LB', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['LB'].to_numpy())
                    self.data.items[-1].set_extra_data(name = 'dL', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['dL'].to_numpy())
        
        #now remove the empty datasets
        empty_data_index = []
        i=0
        for each in self.data.items:
            if len(each.x_raw) == 0:
                empty_data_index.append(i)
            i += 1
        for i in range(len(empty_data_index)):
            self.data.delete_item(empty_data_index[i])
            for ii in range(len(empty_data_index)):
                if empty_data_index[ii]>empty_data_index[i]:
                    empty_data_index[ii] = empty_data_index[ii]-1
                else:
                    pass
        #update script_module
        self.script_module.__dict__['data'] = self.data
        #update the view
        self.update_table_widget_data()
        self.update_combo_box_dataset()
        self.update_plot_data_view()

    def compile_script(self):
        assert self.script!='', 'script file not loaded!'
        exec(self.script, self.script_module.__dict__)
        # print(self.script_module.__dict__)

    

    def update_table_widget_data(self):
        self.tableWidget_data.setRowCount(len(self.data))
        self.tableWidget_data.setColumnCount(4)
        self.tableWidget_data.setHorizontalHeaderLabels(['DataID','Show','Use','Errors'])
        # self.tableWidget_pars.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        for i in range(len(self.data)):
            current_data = self.data[i]
            name = current_data.name
            for j in range(4):
                if j == 0:
                    qtablewidget = QTableWidgetItem(name)
                    self.tableWidget_data.setItem(i,j,qtablewidget)
                else:
                    check_box = QCheckBox()
                    #self.show_checkBox_list.append(check_box)
                    check_box.setChecked(True)
                    check_box.stateChanged.connect(self.update_plot_data_view)
                    self.tableWidget_data.setCellWidget(i,j,check_box)

    def update_combo_box_dataset(self):
        new_items = [each.name for each in self.data]
        self.comboBox_dataset.clear()
        self.comboBox_dataset.addItems(new_items)

    def update_data_view(self):
        dataset_name = self.comboBox_dataset.currentText()
        dataset = None
        for each in self.data:
            if each.name == dataset_name:
                dataset = each
                break
            else:
                pass
        column_labels_main = ['x','y','error']
        extra_labels = ['h', 'k', 'dL', 'LB']
        all_labels = ['x','y','error','h','k','dL','LB','mask']
        self.tableWidget_data_view.setRowCount(len(dataset.x))
        self.tableWidget_data_view.setColumnCount(len(all_labels))
        self.tableWidget_data_view.setHorizontalHeaderLabels(all_labels)
        for i in range(len(dataset.x)):
            for j in range(len(all_labels)):
                if all_labels[j] in column_labels_main:
                    # print(getattr(dataset,'x')[i])
                    qtablewidget = QTableWidgetItem(str(getattr(dataset,all_labels[j])[i]))
                elif all_labels[j] in extra_labels:
                    qtablewidget = QTableWidgetItem(str(dataset.get_extra_data(all_labels[j])[i]))
                else:
                    qtablewidget = QTableWidgetItem('True')
                self.tableWidget_data_view.setItem(i,j,qtablewidget)
    
    def update_structure_view(self):
        xyz = self.script_module.sample.extract_xyz(1)
        self.widget_edp.show_structure(xyz)



    def append_data(self):
        pass

    def delete_data(self):
        pass

    def save_data(self):
        pass

    def calculate(self):
        pass

    def change_plot_style(self):
        pass

    def load_script(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","script Files (*.py);;text Files (*.txt)", options=options)
        if fileName:
            with open(fileName,'r') as f:
                self.plainTextEdit_script.setPlainText(f.read())
        self.script = (self.plainTextEdit_script.toPlainText())
        self.compile_script()

    def save_script(self):
        pass

    def load_par(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","Table Files (*.tab);;text Files (*.txt)", options=options)
        vertical_labels = []
        if fileName:
            with open(fileName,'r') as f:
                lines = f.readlines()
                # self.parameters.set_ascii_input(f)
                lines = [each for each in lines if not each.startswith('#')]
                how_many_pars = len(lines)
                self.tableWidget_pars.setRowCount(how_many_pars)
                self.tableWidget_pars.setColumnCount(6)
                self.tableWidget_pars.setHorizontalHeaderLabels(['Parameter','Value','Fit','Min','Max','Error'])
                # self.tableWidget_pars.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
                for i in range(len(lines)):
                    line = lines[i]
                    items = line.rstrip().rsplit('\t')
                    j = 0
                    if items[0] == '':
                        vertical_labels.append('')
                        j += 1
                    else:
                        #add items to parameter attr
                        self.parameters.data.append([items[0],float(items[1]),items[2]=='True',float(items[3]), float(items[4]),items[5]])
                        #add items to table view
                        if len(vertical_labels)==0:
                            vertical_labels.append('1')
                        else:
                            if vertical_labels[-1] != '':
                                vertical_labels.append('{}'.format(int(vertical_labels[-1])+1))
                            else:
                                vertical_labels.append('{}'.format(int(vertical_labels[-2])+1))
                        for item in items:
                            if j == 2:
                                check_box = QCheckBox()
                                check_box.setChecked(item=='True')
                                self.tableWidget_pars.setCellWidget(i,2,check_box)
                            else:
                                qtablewidget = QTableWidgetItem(item)
                                # qtablewidget.setTextAlignment(Qt.AlignCenter)
                                if j == 0:
                                    qtablewidget.setFont(QFont('Times',10,QFont.Bold))
                                elif j == 1:
                                    qtablewidget.setForeground(QBrush(QColor(255,0,255)))
                                self.tableWidget_pars.setItem(i,j,qtablewidget)
                            j += 1 
        self.tableWidget_pars.resizeColumnsToContents()
        self.tableWidget_pars.resizeRowsToContents()
        self.tableWidget_pars.setShowGrid(False)
        self.tableWidget_pars.setVerticalHeaderLabels(vertical_labels)
        print(self.parameters.data)
        

    def save_par(self):
        pass

#Some Exception definition for errorpassing
class GenericError(Exception):
    ''' Just a empty class used for inheritance. Only useful
    to check if the errors are originating from the model library.
    All these errors are controllable. If they not originate from
    this class something has passed trough and that should be impossible '''
    pass

class ParameterError(GenericError):
    ''' Class for yielding Parameter errors
    '''
    def __init__(self, parameter, parameter_number, error_message, what = -1):
        '''__init__(self, parameter, parameter_number, error_message) --> None

        parameter: the name of the parameter [string]
        parameter_number: the position of the parameter in the list [int]
        error_mesage: pythons error message from the original exception
        set: int to show where the error lies.
            -1 : undefined
             0 : an not find the parameter
             1 : can not evaluate i.e. set the parameter
             2 : value are larger than max
             3 : value are smaller than min
             4 : No parameters to fit
        '''
        self.parameter = parameter
        self.parameter_number = parameter_number
        self.error_message = error_message
        self.what = what

    def __str__(self):
        ''' __str__(self) --> text [string]
        Yields a human readable description of the problem
        '''
        text = ''
        text += 'Parameter number %i, %s, '%(self.parameter_number+1,\
            self.parameter)

        # Take care of the different cases
        if self.what == 0:
            text += 'could not be found. Check the spelling.\n'
        elif self.what == 1:
            text += 'could not be evaluated. Check the code of the function.\n'
        elif self.what == 2:
            text += 'is larger than the value in the max column.\n'
        elif self.what == 3:
            text += 'is smaller than the value in the min column\n'
        elif self.what == 4:
            text = 'There are no parameter selcted to be fitted.\n' + \
                    'Select the parameters you want to fit by checking the ' +\
                    'boxes in the fit column, folder grid'
        else:
            text += 'yielded an undefined error. Check the Python output\n'

        if self.error_message != 'None':
            text += '\nPython error output:\n' + self.error_message

        return text

class ModelError(GenericError):
    ''' Class for yielding compile or evaluation errors in the model text
    '''
    def __init__(self, error_message, where):
        '''__init__(self, error_message, where = -1) --> None

        error_mesage: pythons error message from the original exception
        where: integer describing where the error was raised.
                -1: undef
                 0: compile error
                 1: evaulation error
        '''
        self.error_message = error_message
        self.where = where

    def __str__(self):
        ''' __str__(self) --> text [string]
        Yields a human readable description of the problem
        '''
        text = ''
        if self.where == 0:
            text += 'It was not possible to compile the model script.\n'
        elif self.where == 1:
            text += 'It was not possible to evaluate the model script.\n'\
                    + 'Check the Sim function.\n'
        elif self.where == -1:
            text += 'Undefined error from the Model. See below.\n'

        text += '\n' + self.error_message

        return text

class FomError(GenericError):
    '''Error class for the fom evaluation'''
    def __init__(self, error_message):
        ''' __init__(self, error_message) --> None'''
        self.error_message = error_message

    def __str__(self):
        text = 'Could not evaluate the FOM function. See python output.\n'\
            + '\n' + self.error_message
        return text

class IOError(GenericError):
    ''' Error class for input output, mostly concerning files'''

    def __init__(self, error_message, file = ''):
        '''__init__(self, error_message)'''
        self.error_message = error_message
        self.file = file

    def __str__(self):
        text = 'Input/Output error for file:\n' + self.file +\
                '\n\n Python error:\n ' + self.error_message
        return text


# Some small default function that are needed for initilization

def default_fom_func(simulated_data, data):
    '''
    The default fom function. Its just a dummy so far dont use it!
    '''
    return sum([abs(d.y-sim_d).sum() for sim_d, d \
                in zip(simulated_data,data)])

if __name__ == "__main__":
    QApplication.setStyle("windows")
    app = QApplication(sys.argv)
    myWin = MyMainWindow()
    myWin.setWindowIcon(QtGui.QIcon('dafy.PNG'))
    hightlight = syntax_pars.PythonHighlighter(myWin.plainTextEdit_script.document())
    myWin.plainTextEdit_script.show()
    myWin.plainTextEdit_script.setPlainText(myWin.plainTextEdit_script.toPlainText())
    myWin.show()
    sys.exit(app.exec_())