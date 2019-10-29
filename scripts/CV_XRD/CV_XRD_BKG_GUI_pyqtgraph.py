import sys,os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
#from open_file import *
from PyQt5 import uic
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
from VisualizationEnginePool import plot_bkg_fit_gui_pyqtgraph,replot_bkg_profile
import time
import matplotlib
matplotlib.use("Qt5Agg")
import pyqtgraph as pg
from PyQt5 import QtCore
from PyQt5.QtWidgets import QCheckBox, QRadioButton

#from matplotlib.backends.backend_qt5agg import (NavigationToolbar2QT as NavigationToolbar)

class MyMainWindow(QMainWindow):
    def __init__(self, parent = None):
        super(MyMainWindow, self).__init__(parent)
        pg.setConfigOptions(imageAxisOrder='row-major')
        pg.mkQApp()
        uic.loadUi('C:\\apps\\DaFy_P23\\scripts\\CV_XRD\\ctr_bkg_pyqtgraph3.ui',self)

        self.app_ctr=run_app()
        #self.app_ctr.run()
        self.current_image_no = 0
        self.current_scan_number = None
         
        #self.setupUi(self)
        self.stop = False
        self.open.clicked.connect(self.load_file)
        self.stopBtn.clicked.connect(self.stop_func)
        self.save.clicked.connect(self.save_file)
        self.plot.clicked.connect(self.plot_figure)
        self.runstepwise.clicked.connect(self.plot_)
        #self.update_poly_order(init_step = True)
        for each in self.groupBox_2.findChildren(QCheckBox):
            each.released.connect(self.update_poly_order)
        for each in self.groupBox_cost_func.findChildren(QRadioButton):
            each.toggled.connect(self.update_cost_func)
        self.doubleSpinBox_ss_factor.valueChanged.connect(self.update_ss_factor)
        self.setup_image()

    def update_poly_order(self, init_step = False):
        ord_total = 0
        i=1
        for each in self.groupBox_2.findChildren(QCheckBox):
            ord_total += int(bool(each.checkState()))*int(each.text())
            i+=i
        self.app_ctr.bkg_sub.update_integration_order(ord_total)
        #print(self.app_ctr.bkg_sub.ord_cus_s)
        
        if not init_step:
            self.updatePlot()

    def update_cost_func(self, init_step = False):
        for each in self.groupBox_cost_func.findChildren(QRadioButton):
            if each.isChecked():
                self.app_ctr.bkg_sub.update_integration_function(each.text())
                break
        try:
            self.updatePlot()
        except:
            pass

    def update_ss_factor(self, init_step = False):
        self.app_ctr.bkg_sub.update_ss_factor(self.doubleSpinBox_ss_factor.value())
        #print(self.app_ctr.bkg_sub.ss_factor)
        try:
            self.updatePlot()
        except:
            pass

    def setup_image(self):
        # Interpret image data as row-major instead of col-major
        global img, roi, data, p2, isoLine, iso

        win = self.widget_image
        win.setWindowTitle('pyqtgraph example: Image Analysis')

        # A plot area (ViewBox + axes) for displaying the image
        p1 = win.addPlot()

        # Item for displaying image data
        img = pg.ImageItem()
        self.img_pyqtgraph = img
        p1.addItem(img)

        # Custom ROI for selecting an image region
        roi = pg.ROI([100, 100], [100, 100])
        self.roi = roi
        roi.addScaleHandle([0.5, 1], [0.5, 0.5])
        roi.addScaleHandle([0, 0.5], [0.5, 0.5])
        p1.addItem(roi)
        #roi.setZValue(10)  # make sure ROI is drawn above image

        # Isocurve drawing
        iso = pg.IsocurveItem(level=0.8, pen='g')
        iso.setParentItem(img)
        self.iso = iso
        
        #iso.setZValue(5)

        # Contrast/color control
        hist = pg.HistogramLUTItem()
        self.hist = hist
        hist.setImageItem(img)
        win.addItem(hist)

        # Draggable line for setting isocurve level
        isoLine = pg.InfiniteLine(angle=0, movable=True, pen='g')
        hist.vb.addItem(isoLine)
        hist.vb.setMouseEnabled(y=False) # makes user interaction a little easier
        isoLine.setValue(0.8)
        isoLine.setZValue(100000) # bring iso line above contrast controls

        # Another plot area for displaying ROI data
        win.nextRow()
        p2 = win.addPlot(colspan=2)
        p2.setMaximumHeight(200)
        #p2.setLogMode(y = True)


        # plot to show intensity over time
        win.nextRow()
        p3 = win.addPlot(colspan=2)
        p3.setMaximumHeight(200)


        # plot to show intensity over time
        win.nextRow()
        p4 = win.addPlot(colspan=2)
        p4.setMaximumHeight(200)

        # Generate image data
        #data = np.random.normal(size=(500, 600))
        #data[20:80, 20:80] += 2.
        #data = pg.gaussianFilter(data, (3, 3))
        #data += np.random.normal(size=(500, 600)) * 0.1
        #img.setImage(data)
        ##hist.setLevels(data.min(), data.max())

        # build isocurves from smoothed data
        ##iso.setData(pg.gaussianFilter(data, (2, 2)))

        # set position and scale of image
        #img.scale(0.2, 0.2)
        #img.translate(-50, 0)

        # zoom to fit imageo
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        p1.autoRange()  


        # Callbacks for handling user interaction
        def updatePlot():
            #global data
            try:
                selected = roi.getArrayRegion(self.app_ctr.bkg_sub.img, self.img_pyqtgraph)
            except:
                #selected = roi.getArrayRegion(data, self.img_pyqtgraph)
                pass
            p2.plot(selected.sum(axis=0), clear=True)
            self.reset_peak_center_and_width()
            self.app_ctr.run_update()
            ##update iso curves
            x, y = [int(each) for each in self.roi.pos()]
            w, h = [int(each) for each in self.roi.size()]
            self.iso.setData(pg.gaussianFilter(self.app_ctr.bkg_sub.img[y:(y+h),x:(x+w)], (2, 2)))
            self.iso.setPos(x,y)
            isoLine.setValue(self.app_ctr.bkg_sub.img[y:(y+h),x:(x+w)].mean())
            #plot others
            plot_bkg_fit_gui_pyqtgraph(self.p2, self.p3, self.p4,self.app_ctr.data, self.app_ctr.bkg_sub)


        roi.sigRegionChanged.connect(updatePlot)
        self.updatePlot = updatePlot

        def updateIsocurve():
            global isoLine, iso
            iso.setLevel(isoLine.value())

        self.updateIsocurve = updateIsocurve

        isoLine.sigDragged.connect(updateIsocurve)

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
            self.app_ctr.run(self.lineEdit.text())
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
            if self.app_ctr.bkg_sub.img is not None:
                self.img_pyqtgraph.setImage(self.app_ctr.bkg_sub.img)
                self.p1.autoRange() 
                self.hist.setLevels(self.app_ctr.bkg_sub.img.min(), self.app_ctr.bkg_sub.img.mean()*10)
                #self.iso.setData(pg.gaussianFilter(self.app_ctr.bkg_sub.img[200:400,150:250], (2, 2)))

                self.updatePlot()

            if return_value:
                if self.app_ctr.img_loader.scan_number!=self.current_scan_number:
                    self.current_scan_number = self.app_ctr.img_loader.scan_number 
                    self.current_image_no = 0
                else:
                    self.current_image_no += 1
                #plot_bkg_fit_gui(self.MplWidget.canvas.ax_img, self.MplWidget.canvas.ax_profile, self.MplWidget.canvas.ax_ctr, self.MplWidget.canvas.ax_pot,self.app_ctr.data, self.app_ctr.bkg_sub)
                #self.MplWidget.canvas.figure.tight_layout()
                #self.MplWidget.canvas.draw()
                self.statusbar.clearMessage()
                self.statusbar.showMessage('Working on scan{}: we are now at frame{} of {} frames in total!'.format(self.current_scan_number,self.current_image_no,self.app_ctr.img_loader.total_frame_number))
                self.progressBar.setValue(self.current_image_no/float(self.app_ctr.img_loader.total_frame_number)*100)

            else:
                self.stop_func()

    def update_plot(self):
        img = self.app_ctr.run_update()
        plot_bkg_fit_gui_pyqtgraph(self.p2, self.p3, self.p4,self.app_ctr.data, self.app_ctr.bkg_sub)
        #self.MplWidget.canvas.figure.tight_layout()
        #self.MplWidget.canvas.draw()

    def reset_peak_center_and_width(self):
        roi_size = [int(each/2) for each in self.roi.size()][::-1]
        roi_pos = [int(each) for each in self.roi.pos()][::-1]
        #roi_pos[0] = self.app_ctr.cen_clip[0]*2-roi_pos[0]
        #new_center = [roi_pos[0]-roi_size[0],roi_pos[1]+roi_size[1]]
        #roi_pos[0] = self.app_ctr.cen_clip[0]*2-roi_pos[0]
        new_center = [roi_pos[0]+roi_size[0],roi_pos[1]+roi_size[1]]
        self.app_ctr.bkg_sub.center_pix = new_center
        self.app_ctr.bkg_sub.row_width = roi_size[1]
        self.app_ctr.bkg_sub.col_width = roi_size[0]

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