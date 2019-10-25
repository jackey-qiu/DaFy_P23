import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from open_file import *
from mplwidget import MplWidget
import random
import numpy as np
import matplotlib.pyplot as plt

class MyMainWindow(QMainWindow, Ui_openfile):
    def __init__(self, parent = None):
        super(MyMainWindow, self).__init__(parent)
        self.setupUi(self)
        self.open.clicked.connect(self.load_file)
        self.save.clicked.connect(self.save_file)
        self.plot.clicked.connect(self.plot_figure)

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
        timer.timeout.connect(self.plot_figure_)
        timer.start(1000)

    def plot_figure_(self):
        if self.lineEdit.text()=='True':
            fs = 500
            f = random.randint(1, 100)
            ts = 1/fs
            length_of_signal = 100
            t = np.linspace(0,1,length_of_signal)
            
            cosinus_signal = np.cos(2*np.pi*f*t)
            sinus_signal = np.sin(2*np.pi*f*t)
            #print(dir(self.MplWidget.canvas.axes.figure))
            #fig = plt.figure()
            #self.MplWidget.canvas.axes.set_figure(fig)
            self.MplWidget.canvas.axes.clear()
            self.MplWidget.canvas.axes.plot(t, cosinus_signal)
            self.MplWidget.canvas.axes.plot(t, sinus_signal)
            self.MplWidget.canvas.axes.legend(('cosinus', 'sinus'),loc='upper right')
            self.MplWidget.canvas.axes.set_title('Cosinus - Sinus Signal')
            self.MplWidget.canvas.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    myWin = MyMainWindow()
    myWin.show()
    sys.exit(app.exec_())