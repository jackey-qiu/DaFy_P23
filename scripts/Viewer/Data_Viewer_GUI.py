import sys,os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
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
        self.open.clicked.connect(self.load_file)
        self.plot.clicked.connect(self.plot_figure)
        self.apply.clicked.connect(self.replot_figure)
        self.PushButton_append_scans.clicked.connect(self.append_scans)
        self.PushButton_fold_or_unfold.clicked.connect(self.fold_or_unfold)
        self.data = None
       
    def load_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","Data Files (*.xlsx);;All Files (*.csv)", options=options)
        if fileName:
            self.lineEdit_data_file.setText(fileName)
            self.data = pd.read_excel(fileName)
        col_labels = 'col_labels\n'+str(list(self.data.columns))+'\n'
        scans = list(set(list(self.data['scan_no'])))
        self.scans = scans
        scans.sort()
        scan_numbers = 'scan_nos\n'+str(scans)+'\n'
        # print(list(self.data[self.data['scan_no']==scans[0]]['phs'])[0])
        self.phs = [list(self.data[self.data['scan_no']==scan]['phs'])[0] for scan in scans]
        phs = 'pHs\n'+str(self.phs)+'\n'
        
        self.textEdit_summary_data.setText('\n'.join([col_labels,scan_numbers,phs]))

    #to fold or unfold the config file editor
    def fold_or_unfold(self):
        text = self.PushButton_fold_or_unfold.text()
        if text == "<":
            self.frame.setVisible(False)
            self.fold_or_unfold.setText(">")
        elif text == ">":
            self.frame.setVisible(True)
            self.fold_or_unfold.setText("<")

    def plot_figure(self):
        # print(dir(self))
        self.mplwidget.canvas.ax_img = self.mplwidget.canvas.figure.add_subplot(121)
        self.mplwidget.canvas.ax_profile = self.mplwidget.canvas.figure.add_subplot(322)
        self.mplwidget.canvas.ax_ctr = self.mplwidget.canvas.figure.add_subplot(324)
        self.mplwidget.canvas.ax_pot = self.mplwidget.canvas.figure.add_subplot(326)
        # self.mplwidget.canvas.axes = self.mplwidget.canvas.figure.add_subplot(111)
        self.mplwidget.canvas.draw()
    def replot_figure(self):
        pass

    def append_scans(self):
        text = self.scan_numbers_append.text()
        text_original = self.scan_numbers_all.text()
        if text_original!='':
            text_new = ','.join([text_original, text])
        else:
            text_new = text
        scans = list(set([int(each) for each in text_new.rstrip().rsplit(',')]))
        scans.sort()
        self.scan_numbers_all.setText(','.join([str(scan) for scan in scans]))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    myWin = MyMainWindow()
    myWin.show()
    sys.exit(app.exec_())