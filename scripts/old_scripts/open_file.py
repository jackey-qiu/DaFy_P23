# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\open_file.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_openfile(object):
    def setupUi(self, openfile):
        openfile.setObjectName("openfile")
        openfile.resize(1407, 1007)
        self.centralwidget = QtWidgets.QWidget(openfile)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.centralwidget.sizePolicy().hasHeightForWidth())
        self.centralwidget.setSizePolicy(sizePolicy)
        self.centralwidget.setObjectName("centralwidget")
        self.MplWidget = MplWidget(self.centralwidget)
        self.MplWidget.setGeometry(QtCore.QRect(600, 90, 751, 801))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.MplWidget.sizePolicy().hasHeightForWidth())
        self.MplWidget.setSizePolicy(sizePolicy)
        self.MplWidget.setObjectName("MplWidget")
        self.horizontalSlider = QtWidgets.QSlider(self.centralwidget)
        self.horizontalSlider.setGeometry(QtCore.QRect(430, 260, 141, 22))
        self.horizontalSlider.setMinimum(-50)
        self.horizontalSlider.setMaximum(50)
        self.horizontalSlider.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider.setObjectName("horizontalSlider")
        self.verticalSlider = QtWidgets.QSlider(self.centralwidget)
        self.verticalSlider.setGeometry(QtCore.QRect(560, 90, 22, 171))
        self.verticalSlider.setMinimum(-50)
        self.verticalSlider.setMaximum(50)
        self.verticalSlider.setOrientation(QtCore.Qt.Vertical)
        self.verticalSlider.setObjectName("verticalSlider")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(550, 70, 53, 16))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(370, 260, 53, 16))
        self.label_4.setObjectName("label_4")
        self.progressBar = QtWidgets.QProgressBar(self.centralwidget)
        self.progressBar.setGeometry(QtCore.QRect(600, 910, 751, 23))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.stopBtn = QtWidgets.QPushButton(self.centralwidget)
        self.stopBtn.setGeometry(QtCore.QRect(1150, 50, 141, 23))
        self.stopBtn.setObjectName("stopBtn")
        self.widget = QtWidgets.QWidget(self.centralwidget)
        self.widget.setGeometry(QtCore.QRect(10, 30, 351, 921))
        self.widget.setObjectName("widget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.open = QtWidgets.QPushButton(self.widget)
        self.open.setObjectName("open")
        self.horizontalLayout.addWidget(self.open)
        self.save = QtWidgets.QPushButton(self.widget)
        self.save.setObjectName("save")
        self.horizontalLayout.addWidget(self.save)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.lineEdit = QtWidgets.QLineEdit(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit.sizePolicy().hasHeightForWidth())
        self.lineEdit.setSizePolicy(sizePolicy)
        self.lineEdit.setObjectName("lineEdit")
        self.verticalLayout.addWidget(self.lineEdit)
        self.textEdit = QtWidgets.QTextEdit(self.widget)
        self.textEdit.setObjectName("textEdit")
        self.verticalLayout.addWidget(self.textEdit)
        self.widget1 = QtWidgets.QWidget(self.centralwidget)
        self.widget1.setGeometry(QtCore.QRect(720, 50, 411, 25))
        self.widget1.setObjectName("widget1")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.widget1)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.plot = QtWidgets.QPushButton(self.widget1)
        self.plot.setObjectName("plot")
        self.horizontalLayout_2.addWidget(self.plot)
        self.runstepwise = QtWidgets.QPushButton(self.widget1)
        self.runstepwise.setObjectName("runstepwise")
        self.horizontalLayout_2.addWidget(self.runstepwise)
        self.widget2 = QtWidgets.QWidget(self.centralwidget)
        self.widget2.setGeometry(QtCore.QRect(420, 140, 121, 81))
        self.widget2.setObjectName("widget2")
        self.formLayout = QtWidgets.QFormLayout(self.widget2)
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.formLayout.setObjectName("formLayout")
        self.label = QtWidgets.QLabel(self.widget2)
        self.label.setObjectName("label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.spinBox_ver = QtWidgets.QSpinBox(self.widget2)
        self.spinBox_ver.setMinimum(-50)
        self.spinBox_ver.setMaximum(50)
        self.spinBox_ver.setObjectName("spinBox_ver")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.spinBox_ver)
        self.label_2 = QtWidgets.QLabel(self.widget2)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.spinBox_hor = QtWidgets.QSpinBox(self.widget2)
        self.spinBox_hor.setMinimum(-50)
        self.spinBox_hor.setMaximum(50)
        self.spinBox_hor.setObjectName("spinBox_hor")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.spinBox_hor)
        openfile.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(openfile)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1407, 21))
        self.menubar.setObjectName("menubar")
        openfile.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(openfile)
        self.statusbar.setObjectName("statusbar")
        openfile.setStatusBar(self.statusbar)

        self.retranslateUi(openfile)
        QtCore.QMetaObject.connectSlotsByName(openfile)

    def retranslateUi(self, openfile):
        _translate = QtCore.QCoreApplication.translate
        openfile.setWindowTitle(_translate("openfile", "MainWindow"))
        self.label_3.setText(_translate("openfile", "Cwidth_off"))
        self.label_4.setText(_translate("openfile", "Rwidth_off"))
        self.stopBtn.setText(_translate("openfile", "Stop"))
        self.open.setText(_translate("openfile", "Load..."))
        self.save.setText(_translate("openfile", "Save as..."))
        self.plot.setText(_translate("openfile", "Run continually"))
        self.runstepwise.setText(_translate("openfile", "Run stepwise"))
        self.label.setText(_translate("openfile", "cen_off_ver"))
        self.label_2.setText(_translate("openfile", "cen_off_hor"))

from mplwidget import MplWidget

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    openfile = QtWidgets.QMainWindow()
    ui = Ui_openfile()
    ui.setupUi(openfile)
    openfile.show()
    sys.exit(app.exec_())

