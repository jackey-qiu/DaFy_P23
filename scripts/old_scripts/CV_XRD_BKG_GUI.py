import sys,os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from open_file import *
from mplwidget import MplWidget
import random
import numpy as np
import matplotlib.pyplot as plt
from DaFy_CTR_BKG_class import run_app
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
from VisualizationEnginePool import plot_bkg_fit_gui
import time
import matplotlib
matplotlib.use("Qt5Agg")
#from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

class MyMainWindow(QMainWindow, Ui_openfile):
    def __init__(self, parent = None):
        super(MyMainWindow, self).__init__(parent)
        self.app_ctr=run_app()
        self.current_image_no = 0
        self.current_scan_number = None
        self.setupUi(self)
        self.stop = False
        self.open.clicked.connect(self.load_file)
        self.stopBtn.clicked.connect(self.stop_func)
        self.save.clicked.connect(self.save_file)
        self.plot.clicked.connect(self.plot_figure)
        self.runstepwise.clicked.connect(self.plot_)
        self.spinBox_hor.valueChanged.connect(self.peak_cen_shift_hor)
        self.spinBox_ver.valueChanged.connect(self.peak_cen_shift_ver)
        self.horizontalSlider.valueChanged.connect(self.row_width_shift)
        self.verticalSlider.valueChanged.connect(self.col_width_shift)

    def stop_func(self):
        if not self.stop:
            self.stop = True
            self.stopBtn.setText('Resume')
        else:
            self.stop = False
            self.stopBtn.setText('Stop')
        
    def load_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            self.lineEdit.setText(fileName)
            with open(fileName,'r') as f:
                self.textEdit.setText(f.read())

    def save_file(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save file", "", "Text documents (*.txt);All files (*.*)")
        text = self.textEdit.toPlainText()
        with open(path, 'w') as f:
            f.write(text)

    def plot_figure(self):
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.plot_)
        timer.start(5)

    def plot_(self):
        #self.app_ctr.set_fig(self.MplWidget.canvas.figure)
        if self.stop:
            pass
        else:
            return_value = self.app_ctr.run_script()
            if return_value:
                if self.app_ctr.img_loader.scan_number!=self.current_scan_number:
                    self.current_scan_number = self.app_ctr.img_loader.scan_number 
                    self.current_image_no = 0
                else:
                    self.current_image_no += 1
                #print(dir(self.MplWidget.canvas))
                plot_bkg_fit_gui(self.MplWidget.canvas.ax_img, self.MplWidget.canvas.ax_profile, self.MplWidget.canvas.ax_ctr, self.MplWidget.canvas.ax_pot,self.app_ctr.data, self.app_ctr.bkg_sub)
                self.MplWidget.canvas.figure.tight_layout()
                self.MplWidget.canvas.draw()
                self.statusbar.clearMessage()
                self.statusbar.showMessage('Working on scan{}: we are now at frame{} of {} frames in total!'.format(self.current_scan_number,self.current_image_no,self.app_ctr.img_loader.total_frame_number))
                self.progressBar.setValue(self.current_image_no/float(self.app_ctr.img_loader.total_frame_number)*100)
            else:
                self.stop_func()

    def update_plot(self):
        #self.app_ctr.set_fig(self.MplWidget.canvas.figure)
        img = self.app_ctr.run_update()
        #print(dir(self.MplWidget.canvas))
        plot_bkg_fit_gui(self.MplWidget.canvas.ax_img, self.MplWidget.canvas.ax_profile, self.MplWidget.canvas.ax_ctr, self.MplWidget.canvas.ax_pot,self.app_ctr.data, self.app_ctr.bkg_sub)
        self.MplWidget.canvas.figure.tight_layout()
        self.MplWidget.canvas.draw()

    def peak_cen_shift_hor(self):
        offset = int(self.spinBox_hor.value())
        self.app_ctr.bkg_sub.update_center_pix_left_and_right(offset)
        #print(self.app_ctr.bkg_sub.center_pix)
        self.update_plot()

    def peak_cen_shift_ver(self):
        offset = int(self.spinBox_ver.value())
        self.app_ctr.bkg_sub.update_center_pix_up_and_down(offset)
        #print(self.app_ctr.bkg_sub.center_pix)
        self.update_plot()

    def row_width_shift(self):
        offset = int(self.horizontalSlider.value())
        self.app_ctr.bkg_sub.update_integration_window_row_width(offset)
        #print(self.app_ctr.bkg_sub.center_pix)
        self.update_plot()

    def col_width_shift(self):
        offset = int(self.verticalSlider.value())
        self.app_ctr.bkg_sub.update_integration_window_column_width(offset)
        #print(self.app_ctr.bkg_sub.center_pix)
        self.update_plot()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    myWin = MyMainWindow()
    myWin.show()
    sys.exit(app.exec_())