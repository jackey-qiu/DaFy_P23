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
import data_superrod as data
import model
import solvergui
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

class RunFit(QtCore.QObject):

    updateplot = QtCore.pyqtSignal(str,object)
    def __init__(self,solver):
        super(RunFit, self).__init__()
        self.solver = solver
        self.running = True

    def run(self):
        if self.running:
            self.solver.optimizer.stop = False
            self.solver.StartFit(self.updateplot)
    
    def stop(self):
        self.solver.optimizer.stop = True

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
        #self.fom_func = chi2bars_2
        #parameters
        #self.parameters = parameters.Parameters()
        #scripts
        #self.script = ''
        #script module
        #self.script_module = types.ModuleType('genx_script_module')
        self.model = model.Model()
        # self.solver = solvergui.SolverController(self)
        self.run_fit = RunFit(solvergui.SolverController(self.model))
        self.fit_thread = QtCore.QThread()
        self.run_fit.moveToThread(self.fit_thread)
        self.run_fit.updateplot.connect(self.update_plot_data_view_upon_simulation)
        self.run_fit.updateplot.connect(self.update_par_during_fit)
        self.run_fit.updateplot.connect(self.update_status)
        # self.run_fit.updateplot.connect(self.update_structure_view)

        self.fit_thread.started.connect(self.run_fit.run)

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
        #self.data = data.DataList()

        #table view for parameters set to selecting row basis
        #self.tableWidget_pars.itemChanged.connect(self.update_par_upon_change)
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
        for i in range(len(self.model.data)):
            if self.tableWidget_data.cellWidget(i,1).isChecked():
                # self.selected_data_profile.plot(self.data[i].x, self.data[i].y, clear = True)
                self.selected_data_profile.plot(self.model.data[i].x, self.model.data[i].y,pen={'color': 'y', 'width': 1},  symbolBrush=(255,0,0), symbolSize=5,symbolPen='w', clear = (len(plot_data_index) == 0))
                plot_data_index.append(i)
        self.selected_data_profile.setLogMode(x=False,y=True)
        self.selected_data_profile.autoRange()

    def update_plot_data_view_upon_simulation(self):
        plot_data_index = []
        for i in range(len(self.model.data)):
            if self.tableWidget_data.cellWidget(i,1).isChecked():
                # self.selected_data_profile.plot(self.data[i].x, self.data[i].y, clear = True)
                self.selected_data_profile.plot(self.model.data[i].x, self.model.data[i].y,pen={'color': 'y', 'width': 0},  symbolBrush=(255,0,0), symbolSize=5,symbolPen='w', clear = (len(plot_data_index) == 0))
                self.selected_data_profile.plot(self.model.data[i].x, self.model.data[i].y_sim,pen={'color': 'r', 'width': 2},  clear = False)
                plot_data_index.append(i)
        self.selected_data_profile.setLogMode(x=False,y=True)
        self.selected_data_profile.autoRange()

    def update_plot(self):
        pass

    def init_new_model(self):
        pass

    def open_model(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","rod file (*.rod);;zip Files (*.rar)", options=options)
        if fileName:
            self.model.load(fileName)
        self.update_table_widget_data()
        self.update_combo_box_dataset()
        self.update_plot_data_view()
        self.update_par_upon_load()
        self.update_script_upon_load()

    def save_model(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save file", "", "rod file (*.rod);zip files (*.rar)")
        if path:
            self.model.save(path)

    def simulate_model(self):
        # self.update_par_upon_change()
        self.model.simulate()
        '''
        self.compile_script()
        # self.update_pars()
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
        '''
        self.update_plot_data_view_upon_simulation()
        self.update_structure_view()

    def run_model(self):
        # self.solver.StartFit()
        
        self.fit_thread.start()

    def stop_model(self):
        self.run_fit.stop()
        self.fit_thread.terminate()

    def load_data(self, loader = 'ctr'):
        exec('self.load_data_{}()'.format(loader))

    def load_data_ctr(self):
        #8 columns in total
        #X, H, K, Y, I, eI, LB, dL 
        #for CTR data, X column is L column, Y column all 0
        #for RAXR data, X column is energy column, Y column is L column
        # self.data = data.DataList()
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
                    self.model.data.add_new(name = name)
                    self.model.data.items[-1].x = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['X'].to_numpy()
                    self.model.data.items[-1].y = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['I'].to_numpy()
                    self.model.data.items[-1].error = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['eI'].to_numpy()
                    self.model.data.items[-1].x_raw = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['X'].to_numpy()
                    self.model.data.items[-1].y_raw = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['I'].to_numpy()
                    self.model.data.items[-1].error_raw = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['eI'].to_numpy()
                    self.model.data.items[-1].set_extra_data(name = 'h', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['h'].to_numpy())
                    self.model.data.items[-1].set_extra_data(name = 'k', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['k'].to_numpy())
                    self.model.data.items[-1].set_extra_data(name = 'Y', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['Y'].to_numpy())
                    self.model.data.items[-1].set_extra_data(name = 'LB', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['LB'].to_numpy())
                    self.model.data.items[-1].set_extra_data(name = 'dL', value = data_loaded_pd[(data_loaded_pd['h']==h_temp) & (data_loaded_pd['k']==k_temp)]['dL'].to_numpy())
        
        #now remove the empty datasets
        empty_data_index = []
        i=0
        for each in self.model.data.items:
            if len(each.x_raw) == 0:
                empty_data_index.append(i)
            i += 1
        for i in range(len(empty_data_index)):
            self.model.data.delete_item(empty_data_index[i])
            for ii in range(len(empty_data_index)):
                if empty_data_index[ii]>empty_data_index[i]:
                    empty_data_index[ii] = empty_data_index[ii]-1
                else:
                    pass
        #update script_module
        #self.model.script_module.__dict__['data'] = self.data
        #update the view
        self.update_table_widget_data()
        self.update_combo_box_dataset()
        self.update_plot_data_view()

    def update_table_widget_data(self):
        self.tableWidget_data.clear()
        self.tableWidget_data.setRowCount(len(self.model.data))
        self.tableWidget_data.setColumnCount(4)
        self.tableWidget_data.setHorizontalHeaderLabels(['DataID','Show','Use','Errors'])
        # self.tableWidget_pars.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        for i in range(len(self.model.data)):
            current_data = self.model.data[i]
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
        new_items = [each.name for each in self.model.data]
        self.comboBox_dataset.clear()
        self.comboBox_dataset.addItems(new_items)

    def update_data_view(self):
        dataset_name = self.comboBox_dataset.currentText()
        dataset = None
        for each in self.model.data:
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
        xyz = self.model.script_module.sample.extract_xyz(1)
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
        self.model.script = (self.plainTextEdit_script.toPlainText())
        #self.compile_script()

    def update_script_upon_load(self):
        self.plainTextEdit_script.setPlainText(self.model.script)

    def save_script(self):
        pass

    def update_par_upon_load(self):

        vertical_labels = []
        lines = self.model.parameters.data
        how_many_pars = len(lines)
        self.tableWidget_pars.clear()
        self.tableWidget_pars.setRowCount(how_many_pars)
        self.tableWidget_pars.setColumnCount(6)
        self.tableWidget_pars.setHorizontalHeaderLabels(['Parameter','Value','Fit','Min','Max','Error'])
        # self.tableWidget_pars.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        for i in range(len(lines)):
            items = lines[i]
            #items = line.rstrip().rsplit('\t')
            j = 0
            if items[0] == '':
                #self.model.parameters.data.append([items[0],0,False,0, 0,'-'])
                vertical_labels.append('')
                j += 1
            else:
                #add items to parameter attr
                #self.model.parameters.data.append([items[0],float(items[1]),items[2]=='True',float(items[3]), float(items[4]),items[5]])
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
                        check_box.setChecked(item==True)
                        self.tableWidget_pars.setCellWidget(i,2,check_box)
                    else:
                        qtablewidget = QTableWidgetItem(str(item))
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
                        self.model.parameters.data.append([items[0],0,False,0, 0,'-'])
                        vertical_labels.append('')
                        j += 1
                    else:
                        #add items to parameter attr
                        self.model.parameters.data.append([items[0],float(items[1]),items[2]=='True',float(items[3]), float(items[4]),items[5]])
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

    @QtCore.pyqtSlot(str,object)
    def update_par_during_fit(self,string,model):
        #labels = [data[0] for each in self.model.parameters.data]
        for i in range(len(model.parameters.data)):
            if model.parameters.data[i][0]!='':
                # print(self.model.parameters.data[i][0])
                #print(len(self.model.parameters.data))
                print(model.parameters.data[i][0])
                item_temp = self.tableWidget_pars.item(i,1)
                #print(type(item_temp))
                item_temp.setText(str(model.parameters.data[i][1]))
        self.tableWidget_pars.resizeColumnsToContents()
        self.tableWidget_pars.resizeRowsToContents()
        self.tableWidget_pars.setShowGrid(False)
        # self.update_structure_view()

    def update_par_upon_change(self):
        self.model.parameters.data = []
        for each_row in range(self.tableWidget_pars.rowCount()):
            if self.tableWidget_pars.item(each_row,0)==None:
                items = ['',0,False,0,0,'-']
            elif self.tableWidget_pars.item(each_row,0).text()=='':
                items = ['',0,False,0,0,'-']
            else:
                # print(each_row,type(self.tableWidget_pars.item(each_row,0)))
                items = [self.tableWidget_pars.item(each_row,0).text()] + [float(self.tableWidget_pars.item(each_row,i).text()) for i in [1,3,4]] + [self.tableWidget_pars.item(each_row,5).text()]
                items.insert(2, self.tableWidget_pars.cellWidget(each_row,2).isChecked())
                self.model.parameters.data.append(items)

    @QtCore.pyqtSlot(str,object)
    def update_status(self,string,model):
        self.statusbar.clearMessage()
        self.statusbar.showMessage(string)

    def save_par(self):
        pass

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