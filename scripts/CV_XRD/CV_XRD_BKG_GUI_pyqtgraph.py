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
        uic.loadUi(os.path.join(DaFy_path,'scripts','CV_XRD','ctr_bkg_pyqtgraph4_new.ui'),self)
        self.setWindowTitle('Data analysis factory: CTR data analasis')
        self.app_ctr=run_app()
        #self.app_ctr.run()
        self.current_image_no = 0
        self.current_scan_number = None
         
        #self.setupUi(self)
        self.stop = False
        self.open.clicked.connect(self.load_file)
        self.reload.clicked.connect(self.rload_file)
        self.stopBtn.clicked.connect(self.stop_func)
        self.saveas.clicked.connect(self.save_file_as)
        self.save.clicked.connect(self.save_file)
        self.plot.clicked.connect(self.plot_figure)
        self.runstepwise.clicked.connect(self.plot_)
        self.pushButton_filePath.clicked.connect(self.locate_data_folder)
        self.lineEdit_data_file_name.setText('temp_data.xlsx')
        self.lineEdit_data_file_path.setText(self.app_ctr.data_path)
        #self.lineEdit.setText(self.app_ctr.conf_path_temp)
        self.actionOpenConfig.triggered.connect(self.load_file)
        self.actionSaveConfig.triggered.connect(self.save_file)
        self.actionRun.triggered.connect(self.plot_)
        self.actionStop.triggered.connect(self.stop_func)
        self.actionSaveData.triggered.connect(self.save_data)
        #self.update_poly_order(init_step = True)
        for each in self.groupBox_2.findChildren(QCheckBox):
            each.released.connect(self.update_poly_order)
        for each in self.groupBox_cost_func.findChildren(QRadioButton):
            each.toggled.connect(self.update_cost_func)
        self.pushButton_remove_current_point.clicked.connect(self.remove_data_point)
        self.doubleSpinBox_ss_factor.valueChanged.connect(self.update_ss_factor)
        self.setup_image()
        self.timer_save_data = QtCore.QTimer(self)
        self.timer_save_data.timeout.connect(self.save_data)
        

    def save_data(self):
        data_file = os.path.join(self.lineEdit_data_file_path.text(),self.lineEdit_data_file_name.text())
        try:
            self.app_ctr.save_data_file(data_file)
            self.statusbar.showMessage('Data file is saved as {}!'.format(data_file))
        except:
            self.statusbar.showMessage('Failure to save data file!')

    def remove_data_point(self):
        self.app_ctr.data['mask_ctr'][-1]=False
        self.statusbar.showMessage('Current data point is masked!')
        self.updatePlot2()

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
        self.isoLine = isoLine
        hist.vb.addItem(isoLine)
        hist.vb.setMouseEnabled(y=True) # makes user interaction a little easier
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
            if self.app_ctr.img_loader.frame_number ==0:
                isoLine.setValue(self.app_ctr.bkg_sub.img[y:(y+h),x:(x+w)].mean())
            else:
                pass
            #print(isoLine.value(),self.current_image_no)
            #plot others
            plot_bkg_fit_gui_pyqtgraph(self.p2, self.p3, self.p4,self.app_ctr.data, self.app_ctr.bkg_sub)
            self.lcdNumber_potential.display(self.app_ctr.data['potential'][-1])
            self.lcdNumber_current.display(self.app_ctr.data['current'][-1])
            self.lcdNumber_intensity.display(self.app_ctr.data['peak_intensity'][-1])
            self.lcdNumber_iso.display(isoLine.value())

        def updatePlot_after_remove_point():
            #global data
            try:
                selected = roi.getArrayRegion(self.app_ctr.bkg_sub.img, self.img_pyqtgraph)
            except:
                #selected = roi.getArrayRegion(data, self.img_pyqtgraph)
                pass
            p2.plot(selected.sum(axis=0), clear=True)
            self.reset_peak_center_and_width()
            #self.app_ctr.run_update()
            ##update iso curves
            x, y = [int(each) for each in self.roi.pos()]
            w, h = [int(each) for each in self.roi.size()]
            self.iso.setData(pg.gaussianFilter(self.app_ctr.bkg_sub.img[y:(y+h),x:(x+w)], (2, 2)))
            self.iso.setPos(x,y)
            if self.app_ctr.img_loader.frame_number ==0:
                isoLine.setValue(self.app_ctr.bkg_sub.img[y:(y+h),x:(x+w)].mean())
            else:
                pass
            #print(isoLine.value(),self.current_image_no)
            #plot others
            plot_bkg_fit_gui_pyqtgraph(self.p2, self.p3, self.p4,self.app_ctr.data, self.app_ctr.bkg_sub)
            self.lcdNumber_potential.display(self.app_ctr.data['potential'][-2])
            self.lcdNumber_current.display(self.app_ctr.data['current'][-2])
            self.lcdNumber_intensity.display(self.app_ctr.data['peak_intensity'][-2])
            self.lcdNumber_iso.display(isoLine.value())

        roi.sigRegionChanged.connect(updatePlot)
        self.updatePlot = updatePlot
        self.updatePlot2 = updatePlot_after_remove_point

        def updateIsocurve():
            global isoLine, iso
            iso.setLevel(isoLine.value())
            self.lcdNumber_iso.display(isoLine.value())

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
            #self.app_ctr.run(self.lineEdit.text())
            self.timer_save_data.start(self.spinBox_save_frequency.value()*1000)
            #self.current_image_no = 0
            #self.current_scan_number = self.app_ctr.img_loader.scan_number
            #self.plot_()
            with open(fileName,'r') as f:
                self.textEdit.setText(f.read())

    def locate_data_folder(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            self.lineEdit_data_file_path.setText(os.path.dirname(fileName))

    def rload_file(self):
        self.save_file()
        #self.region_bounds = [0,1]
        try:
            self.app_ctr.run(self.lineEdit.text())
            self.timer_save_data.stop()
            self.timer_save_data.start(self.spinBox_save_frequency.value()*1000)
            self.plot_()
            self.statusbar.showMessage('Initialization succeed!')
        except:
            self.statusbar.showMessage('Initialization failed!')

        #self.app_ctr.run(self.lineEdit.text())
        #self.timer_save_data.stop()
        #self.timer_save_data.start(self.spinBox_save_frequency.value()*1000)
        #self.plot_()

    def save_file_as(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save file", "", "Text documents (*.txt);All files (*.*)")
        text = self.textEdit.toPlainText()
        with open(path, 'w') as f:
            f.write(text)
        self.statusbar.showMessage('Config file is saved as {}!'.format(path))

    def save_file(self):
        text = self.textEdit.toPlainText()
        if text=='':
            self.statusbar.showMessage('Text editor is empty. Config file is not saved!')
        else:
            with open(self.lineEdit.text(), 'w') as f:
                f.write(text)
            self.statusbar.showMessage('Config file is saved with the same file name!')

    def plot_figure(self):
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.plot_)
        self.timer.start(5)

    def plot_(self):
        #self.app_ctr.set_fig(self.MplWidget.canvas.figure)
        if self.stop:
            self.timer.stop()
        else:
            return_value = self.app_ctr.run_script()
            if self.app_ctr.bkg_sub.img is not None:
                #if self.current_scan_number == None:
                #    self.current_scan_number = self.app_ctr.img_loader.scan_number
                self.lcdNumber_scan_number.display(self.app_ctr.img_loader.scan_number)
                self.img_pyqtgraph.setImage(self.app_ctr.bkg_sub.img)
                if self.app_ctr.img_loader.frame_number == 0:
                    self.p1.autoRange() 
                self.hist.setLevels(self.app_ctr.bkg_sub.img.min(), self.app_ctr.bkg_sub.img.mean()*10)
                self.updatePlot()

            if return_value:
                self.statusbar.clearMessage()
                self.statusbar.showMessage('Working on scan{}: we are now at frame{} of {} frames in total!'.format(self.app_ctr.img_loader.scan_number,self.app_ctr.img_loader.frame_number+1,self.app_ctr.img_loader.total_frame_number))
                self.progressBar.setValue((self.app_ctr.img_loader.frame_number+1)/float(self.app_ctr.img_loader.total_frame_number)*100)
                self.lcdNumber_frame_number.display(self.app_ctr.img_loader.frame_number+1)
                #self.app_ctr.img_loader.frame_number
                #self.current_image_no += 1
            else:
                self.timer.stop()
                self.stop = False
                self.stopBtn.setText('Stop')
                self.statusbar.clearMessage()
                self.statusbar.showMessage('Run for scan{} is finished, {} frames in total have been processed!'.format(self.app_ctr.img_loader.scan_number,self.app_ctr.img_loader.total_frame_number))

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