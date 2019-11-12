import sys,os
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
#from open_file import *
from PyQt5 import uic
from mplwidget import MplWidget
import random
import numpy as np
import matplotlib.pyplot as plt
from DaFy_PXRD_class import run_app
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
from VisualizationEnginePool import plot_bkg_fit_gui_pyqtgraph,replot_bkg_profile, plot_pxrd_fit_gui_pyqtgraph
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
        #pg.mkQApp()
        uic.loadUi(os.path.join(DaFy_path,'scripts','CV_XRD','ctr_bkg_pyqtgraph6_new.ui'),self)
        self.setWindowTitle('Data analysis factory: PXRD data analasis')
        self.app_ctr=run_app()
        #self.app_ctr.run()
        self.current_image_no = 0
        self.current_scan_number = None
         
        #self.setupUi(self)
        self.stop = False
        self.open.clicked.connect(self.load_file)
        self.relaunch.clicked.connect(self.relaunch_file)
        self.launch.clicked.connect(self.launch_file)
        self.stopBtn.clicked.connect(self.stop_func)
        self.saveas.clicked.connect(self.save_file_as)
        self.save.clicked.connect(self.save_file)
        self.plot.clicked.connect(self.plot_figure)
        self.runstepwise.clicked.connect(self.plot_)
        self.pushButton_filePath.clicked.connect(self.locate_data_folder)
        self.lineEdit_data_file_name.setText('temp_data.xlsx')
        self.lineEdit_data_file_path.setText(self.app_ctr.data_path)
        #self.lineEdit.setText(self.app_ctr.conf_path_temp)
        #self.update_poly_order(init_step = True)
        for each in self.groupBox_2.findChildren(QCheckBox):
            each.released.connect(self.update_poly_order)
        for each in self.groupBox_cost_func.findChildren(QRadioButton):
            each.toggled.connect(self.update_cost_func)
        self.pushButton_remove_current_point.clicked.connect(self.remove_data_point)
        self.doubleSpinBox_ss_factor.valueChanged.connect(self.update_ss_factor)
        self.comboBox_p4.activated.connect(self.select_source_for_plot_p4)
        self.comboBox_p5.activated.connect(self.select_source_for_plot_p5)
        self.p4_data_source = self.comboBox_p4.currentText()
        self.p5_data_source = self.comboBox_p5.currentText()
        setattr(self.app_ctr,'p4_data_source',self.comboBox_p4.currentText())
        setattr(self.app_ctr,'p5_data_source',self.comboBox_p5.currentText())
        #self.setup_image()
        self.timer_save_data = QtCore.QTimer(self)
        self.timer_save_data.timeout.connect(self.save_data)


    def select_source_for_plot_p4(self):
        if self.app_ctr.time_scan:
            self.app_ctr.p4_data_source = self.comboBox_p4.currentText()
            self.updatePlot()
        else:
            pass

    def select_source_for_plot_p5(self):
        if self.app_ctr.time_scan:
            self.app_ctr.p5_data_source = self.comboBox_p5.currentText()
            self.updatePlot()
        else:
            pass

    def save_data(self):
        data_file = os.path.join(self.lineEdit_data_file_path.text(),self.lineEdit_data_file_name.text())
        #self.app_ctr.save_data_file(data_file)
        try:
            self.app_ctr.save_data_file(data_file)
            self.statusbar.showMessage('Data file is saved as {}!'.format(data_file))
        except:
            self.statusbar.showMessage('Failure to save data file!')

    def remove_data_point(self):
        pass
        '''
        self.app_ctr.data['mask_ctr'][-1]=False
        self.statusbar.showMessage('Current data point is masked!')
        self.updatePlot2()
        '''

    def update_poly_order(self, init_step = False):
        ord_total = 0
        i=1
        for each in self.groupBox_2.findChildren(QCheckBox):
            ord_total += int(bool(each.checkState()))*int(each.text())
            i+=i
        self.app_ctr.kwarg_bkg['ord_cus_s']=ord_total
        #print(self.app_ctr.bkg_sub.ord_cus_s)
        
        if not init_step:
            self.updatePlot()

    def update_cost_func(self, init_step = False):
        for each in self.groupBox_cost_func.findChildren(QRadioButton):
            if each.isChecked():
                self.app_ctr.kwarg_bkg['fct']=each.text()
                break
        try:
            self.updatePlot()
        except:
            pass

    def update_ss_factor(self, init_step = False):
        self.app_ctr.kwarg_bkg['ss_factor']=self.doubleSpinBox_ss_factor.value()
        #print(self.app_ctr.bkg_sub.ss_factor)
        try:
            self.updatePlot()
        except:
            pass

    def setup_image(self):
        # Interpret image data as row-major instead of col-major
        global img, roi, data, isoLine, iso,region
        win = self.widget_image
        #win.setWindowTitle('pyqtgraph example: Image Analysis')
        #iso.setZValue(5)
        img = pg.ImageItem()
        # Contrast/color control
        hist = pg.HistogramLUTItem()
        self.hist = hist
        hist.setImageItem(img)
        win.addItem(hist)

        # A plot area (ViewBox + axes) for displaying the image
        if self.app_ctr.time_scan:
            p1 = win.addPlot(row=1,col=0,rowspan=4,colspan=1)
        else:
            p1 = win.addPlot(row=1,col=0,rowspan=3,colspan=1)
        p1.setMaximumWidth(300)

        # Item for displaying image data
        self.img_pyqtgraph = img
        p1.addItem(img)

        # Custom ROI for selecting an image region
        
        roi = pg.ROI([0, 0], [self.app_ctr.dim_detector[1]-2*self.app_ctr.hor_offset, self.app_ctr.dim_detector[0]-2*self.app_ctr.ver_offset])
        self.roi = roi
        roi.addScaleHandle([0.5, 1], [0.5, 0.5])
        roi.addScaleHandle([0, 0.5], [0.5, 0.5])
        p1.addItem(roi)
        #roi.setZValue(10)  # make sure ROI is drawn above image

        # Isocurve drawing
        iso = pg.IsocurveItem(level=0.8, pen='g')
        iso.setParentItem(img)
        self.iso = iso
        # Draggable line for setting isocurve level
        isoLine = pg.InfiniteLine(angle=0, movable=True, pen='g')
        self.isoLine = isoLine
        hist.vb.addItem(isoLine)
        hist.vb.setMouseEnabled(y=True) # makes user interaction a little easier
        isoLine.setValue(0.)
        isoLine.setZValue(100000) # bring iso line above contrast controls

        # Another plot area for displaying ROI data
        #win.nextColumn()
        p2 = win.addPlot(0,1,rowspan=1,colspan=1, title = 'image profile (vertical direction)')
        #p2.setMaximumWidth(400)
        #p2.setMaximumHeight(200)
        #p2.setLogMode(y = True)

        region = pg.LinearRegionItem()
        region.setZValue(10)
        region.setRegion([10, 15])
        # Add the LinearRegionItem to the ViewBox, but tell the ViewBox to exclude this 
        # item when doing auto-range calculations.


        #monitor window
        p2_mon = win.addPlot(0,2,rowspan=1, colspan=2,title='Monitor window')

        def update():
            region.setZValue(10)
            region.show()
            minX, maxX = region.getRegion()
            #save the bounds of shape area for plotting in monitor window
            self.region_bounds = [minX, maxX]

            x, y = p2_mon.listDataItems()[0].xData, p2_mon.listDataItems()[0].yData 
            data_range =[np.argmin(abs(x-minX)),np.argmin(abs(x-maxX))]
            if data_range[1]>len(y):
                data_range[1] = len(y)
            minY = min(y[data_range[0]:data_range[1]])
            maxY = max(y[data_range[0]:data_range[1]])*1.1
            p2_mon.setYRange(minY, maxY, padding=0) 

            p2_mon.setXRange(minX, maxX, padding=0)  
            #print(p2_mon.listDataItems()[0].xData)
            #print(dir(p2_mon))
        region.sigRegionChanged.connect(update)

        # plot to show intensity over time
        #win.nextRow()
        p3 = win.addPlot(1,1,rowspan=1,colspan=3)
        if self.app_ctr.time_scan:
            p2.addItem(region, ignoreBounds=True)
        else:
            p3.addItem(region, ignoreBounds=True)
        #p3.addLegend()
        #p3.setMaximumHeight(200)
        #

        # plot to show intensity over time
        #win.nextRow()
        p4 = win.addPlot(2,1,rowspan=1,colspan=3)
        #p4.setMaximumHeight(200)

        #if self.app_ctr.time_scan:
        p5 = win.addPlot(3,1,rowspan=1,colspan=3)
        p5.setMaximumHeight(200)

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
        self.p5 = p5
        self.p2_mon = p2_mon       
        p1.autoRange()  
        self.legend_p2 = self.p2.addLegend()
        self.legend_p3 = self.p3.addLegend()
        self.legend_p4 = self.p4.addLegend()
        self.legend_p5 = self.p5.addLegend()
        #set axis labels
        self.p2.setLabel('left','Normalized intensity',units='c/s')
        self.p2.setLabel('bottom','2theta angle',units='°')
        self.p2_mon.setLabel('left','Normalized intensity',units='c/s')
        self.p2_mon.setLabel('bottom','2theta angle',units='°')
        self.p3.setLabel('left','Integrated peak intensity',units='c/s')
        #self.p3.setLabel('bottom','frame_number')
        self.p4.setLabel('left','Peak position',units='°')
        #self.p4.setLabel('bottom','frame_number')
        self.p5.setLabel('left','Potential wrt Ag/AgCl',units='V')
        self.p5.setLabel('bottom','frame_number')

        # Callbacks for handling user interaction
        def updatePlot():
            try:
                selected = roi.getArrayRegion(self.app_ctr.img, self.img_pyqtgraph)
            except:
                pass
            self.app_ctr.run_update()
            #remove legend and plot again
            for legend in [self.legend_p2,self.legend_p3,self.legend_p4,self.legend_p5]:
                try:
                    legend.scene().removeItem(legend)
                except:
                    pass
            self.legend_p2 = self.p2.addLegend()
            self.legend_p3 = self.p3.addLegend()
            self.legend_p4 = self.p4.addLegend()
            self.legend_p5 = self.p5.addLegend()

            if self.app_ctr.time_scan:
                plot_pxrd_fit_gui_pyqtgraph([self.p2_mon,self.p2], self.p3, self.p4,self.p5,self.app_ctr)
                self.p2.addItem(region, ignoreBounds=True)#re-add the region item
            else:
                plot_pxrd_fit_gui_pyqtgraph([self.p2_mon,self.p2], self.p3, None,self.p4, self.app_ctr)
                self.p3.addItem(region, ignoreBounds=True)
                #region.setRegion(self.region_bounds)#re-add the region item
            try:
                self.lcdNumber_potential.display(self.app_ctr.data[self.app_ctr.img_loader.scan_number]['potential'][-1])
                self.lcdNumber_current.display(self.app_ctr.data[self.app_ctr.img_loader.scan_number]['current'][-1])
                self.lcdNumber_intensity.display(self.app_ctr.data[self.app_ctr.img_loader.scan_number]['peak_intensity'][-1])
            except:
                pass
            self.lcdNumber_iso.display(isoLine.value())
            
        def updatePlot_after_remove_point():#not implemented for PXRD
            #global data
            try:
                selected = roi.getArrayRegion(self.app_ctr.bkg_sub.img, self.img_pyqtgraph)
            except:
                #selected = roi.getArrayRegion(data, self.img_pyqtgraph)
                pass
            p2.plot(selected.sum(axis=0), clear=True)
            #self.reset_peak_center_and_width()
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
            #self.timer_save_data.start(self.spinBox_save_frequency.value()*1000)
            with open(fileName,'r') as f:
                self.textEdit.setText(f.read())

    def locate_data_folder(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","All Files (*);;Python Files (*.py)", options=options)
        if fileName:
            self.lineEdit_data_file_path.setText(os.path.dirname(fileName))

    def relaunch_file(self):
        self.save_file()     
        try:
            self.app_ctr.run(self.lineEdit.text())
            #self.setup_image()
            self.timer_save_data.stop()
            self.timer_save_data.start(self.spinBox_save_frequency.value()*1000)
            self.plot_()
            self.statusbar.showMessage('Initialization succeed!')
        except:
            self.statusbar.showMessage('Initialization failed!')

    def launch_file(self):
        self.save_file()     
        try:
            self.app_ctr.run(self.lineEdit.text())
            self.setup_image()
            self.timer_save_data.stop()
            self.timer_save_data.start(self.spinBox_save_frequency.value()*1000)
            self.plot_()
            self.launch.setEnabled(False)
            self.statusbar.showMessage('Initialization succeed!')
        except:
            self.statusbar.showMessage('Initialization failed!')

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
        #auto plotting timmer
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.plot_)
        self.timer.start(5)

    def plot_(self):
        #self.app_ctr.set_fig(self.MplWidget.canvas.figure)
        if self.stop:
            self.timer.stop()
        else:
            return_value = self.app_ctr.run_script()
            if self.app_ctr.img is not None:
                #if self.current_scan_number == None:
                #    self.current_scan_number = self.app_ctr.img_loader.scan_number
                self.lcdNumber_scan_number.display(self.app_ctr.img_loader.scan_number)
                self.img_pyqtgraph.setImage(self.app_ctr.img)
                if self.app_ctr.img_loader.frame_number == 0:
                    self.p1.autoRange() 
                self.hist.setLevels(self.app_ctr.img.min(), self.app_ctr.img.mean()*5)
                self.updatePlot()

            if return_value:
                self.statusbar.clearMessage()
                self.statusbar.showMessage('Working on scan{}: we are now at frame{} of {} frames in total!'.format(self.app_ctr.img_loader.scan_number,self.app_ctr.img_loader.frame_number+1,self.app_ctr.img_loader.total_frame_number))
                self.progressBar.setValue((self.app_ctr.img_loader.frame_number+1)/float(self.app_ctr.img_loader.total_frame_number)*100)
                self.lcdNumber_frame_number.display(self.app_ctr.img_loader.frame_number+1)
            else:
                self.timer.stop()
                self.stop = False
                self.stopBtn.setText('Stop')
                self.statusbar.clearMessage()
                self.statusbar.showMessage('Run for scan{} is finished, {} frames in total have been processed!'.format(self.app_ctr.img_loader.scan_number,self.app_ctr.img_loader.total_frame_number))

    def update_plot(self):
        img = self.app_ctr.run_update()
        plot_pxrd_fit_gui_pyqtgraph(self.p2, self.p3, self.p4, self.app_ctr)
        #self.MplWidget.canvas.figure.tight_layout()
        #self.MplWidget.canvas.draw()

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