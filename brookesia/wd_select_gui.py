#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:04:02 2020

@author: matynia
"""

import pandas as pd
import os


from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import time as timer



root_path   = sys.argv[1]
if os.name == 'nt': #different python path on windows
    python_path = root_path.split('lib')[0] + 'python.exe'
else:
    python_path = root_path.split('/lib/')[0] + '/bin/python'

os.chdir(root_path)

import write_conditions_input as wc
import write_kinetic_mech     as wk
import write_results_input    as wr


global workdir_table
global full_workdir_table
workdir_table = pd.read_table("working_dir.txt",sep =';',header = 0, index_col=[0])
full_workdir_table = workdir_table
total_nwd = len(workdir_table['directory'])



class Ui_MainWD(object):
    def setupUi(self, MainWD):
        global sz_w
        global sz_h
        global workdir_table
        global full_workdir_table

        _height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()
        _width = QtWidgets.QDesktopWidget().screenGeometry(-1).width()
        sz_w = _width/1366
        sz_h = _height/768

        _translate = QtCore.QCoreApplication.translate

        MainWD.setObjectName("Select the working directory")
        MainWD.resize(450*sz_w, 400*sz_h)
        MainWD.setWindowTitle(_translate("MainWD", "Please select your working directory"))

        self.centralwidget = QtWidgets.QWidget(MainWD)
        self.centralwidget.setObjectName("centralwidget")
        MainWD.setCentralWidget(self.centralwidget)


        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setGeometry(QtCore.QRect(10*sz_w, 10*sz_h, 430*sz_w, 365*sz_h))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")


        self.scrollArea = QtWidgets.QScrollArea(self.frame)
        self.scrollArea.setGeometry(QtCore.QRect(10*sz_w, 10*sz_h, 420*sz_w, 280*sz_h))
        self.scrollArea.setWidgetResizable(False)
        self.scrollArea.setObjectName("scrollArea")

        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")

        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0*sz_w, 0*sz_h, 400*sz_w, (10+35*total_nwd)*sz_h))

        self.scrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)

#        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 2500*sz_w, (15+181*(num_case+1))*sz_h))


        self.pb_wd_name, self.pb_wd_edit_n, self.pb_wd_edit_p, self.pb_wd_rm, self.rm_list = [], [], [], [], []

        for nwd in range(100):

            self.pb_wd_name.append(QtWidgets.QPushButton(self.scrollAreaWidgetContents))
            self.pb_wd_name[-1].setGeometry(QtCore.QRect(10*sz_w, (10+35*nwd)*sz_h, 190*sz_w, 25*sz_h))
            self.pb_wd_name[-1].setObjectName("pb_wd_name "+str(nwd))
            self.pb_wd_name[-1].clicked.connect(lambda checked, nwd=nwd: self.selected_wd(nwd,MainWD))
#lambda checked, key=key:
            self.pb_wd_edit_n.append(QtWidgets.QPushButton(self.scrollAreaWidgetContents))
            self.pb_wd_edit_n[-1].setGeometry(QtCore.QRect(205*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
            self.pb_wd_edit_n[-1].setObjectName("pb_wd_edit_n "+str(nwd))
            self.pb_wd_edit_n[-1].setText(_translate("MainWD", 'Edit name'))
            self.pb_wd_edit_n[-1].clicked.connect(lambda checked, nwd=nwd: self.edit_wd(nwd,'name'))

            self.pb_wd_edit_p.append(QtWidgets.QPushButton(self.scrollAreaWidgetContents))
            self.pb_wd_edit_p[-1].setGeometry(QtCore.QRect(290*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
            self.pb_wd_edit_p[-1].setObjectName("pb_wd_edit_p "+str(nwd))
            self.pb_wd_edit_p[-1].setText(_translate("MainWD", 'Edit path'))
            self.pb_wd_edit_p[-1].clicked.connect(lambda checked, nwd=nwd: self.edit_wd(nwd,'path'))

            self.pb_wd_rm.append(QtWidgets.QPushButton(self.scrollAreaWidgetContents))
            self.pb_wd_rm[-1].setGeometry(QtCore.QRect(375*sz_w, (10+35*nwd)*sz_h, 25*sz_w, 25*sz_h))
            self.pb_wd_rm[-1].setObjectName("pb_wd_rm "+str(nwd))
            self.pb_wd_rm[-1].setText(_translate("MainWD", '(-)'))#chr(128465)))
            self.pb_wd_rm[-1].clicked.connect(lambda checked, nwd=nwd: self.remove_wd(nwd))

            if nwd <total_nwd:
                self.pb_wd_name[-1].setGeometry(QtCore.QRect(10*sz_w, (10+35*nwd)*sz_h, 190*sz_w, 25*sz_h))
                self.pb_wd_name[-1].setText(_translate("MainWD", workdir_table['working_dir_name'].iloc[nwd]))
                self.pb_wd_edit_n[-1].setGeometry(QtCore.QRect(205*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
                self.pb_wd_edit_p[-1].setGeometry(QtCore.QRect(290*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
                self.pb_wd_rm[-1].setGeometry(QtCore.QRect(375*sz_w, (10+35*nwd)*sz_h, 25*sz_w, 25*sz_h))
                self.rm_list.append(False)
            else:
                self.pb_wd_name[-1].setGeometry(QtCore.QRect(1000*sz_w, (10+35*nwd)*sz_h, 190*sz_w, 25*sz_h))
                self.pb_wd_name[-1].setText(_translate("MainWD", 'not_attributed'))
                self.pb_wd_edit_n[-1].setGeometry(QtCore.QRect(2050*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
                self.pb_wd_edit_p[-1].setGeometry(QtCore.QRect(2900*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
                self.pb_wd_rm[-1].setGeometry(QtCore.QRect(3750*sz_w, (10+35*nwd)*sz_h, 25*sz_w, 25*sz_h))

        self.pB_setnew_wd = QtWidgets.QPushButton(self.frame)
        self.pB_setnew_wd.setGeometry(QtCore.QRect(115*sz_w, 310*sz_h, 200*sz_w, 35*sz_h))
        self.pB_setnew_wd.setObjectName("pB_setnew_wd")
        self.pB_setnew_wd.setText(_translate("MainWD", "New working directory"))
        self.pB_setnew_wd.clicked.connect(self.add_wd)


        self.menubar = QtWidgets.QMenuBar(MainWD)
        self.menubar.setGeometry(QtCore.QRect(0*sz_w, 0*sz_h, 800*sz_w, 20*sz_h))
        self.menubar.setObjectName("menubar")

        MainWD.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWD)

        self.statusbar.setObjectName("statusbar")
        MainWD.setStatusBar(self.statusbar)

        QtCore.QMetaObject.connectSlotsByName(MainWD)



    def selected_wd(self,nwd,MainWD):
        global workdir_table
        global full_workdir_table

        WD_name = full_workdir_table.iloc[nwd,0]
        WD_path = full_workdir_table.iloc[nwd,1]

        os.chdir(root_path)

        fd = open('wd_activ.txt', 'w')
        fd.write(WD_name + ';' + WD_path)

        timer.sleep(0.2)

        MainWD.close()



    def edit_wd(self,nwd,option):
        global workdir_table
        global full_workdir_table
        _translate = QtCore.QCoreApplication.translate

        if option == 'name':
            os.system(python_path + ' ' + root_path + "wd_manage_gui.py True False " + root_path)
        else:
            os.system(python_path + ' ' + root_path + "wd_manage_gui.py False True " + root_path)

        os.chdir(root_path)

        fs = open('wd_manage.txt', 'r')
        txt = fs.readline()
        WD_name = txt.split(';')[0]
        WD_path = txt.split(';')[1]

        if option == 'name':
            workdir_table['working_dir_name'].iloc[nwd] = WD_name
            full_workdir_table['working_dir_name'].iloc[nwd] = WD_name
            self.pb_wd_name[nwd].setText(_translate("MainWD", WD_name))
        else:
            workdir_table['directory'].iloc[nwd]        = WD_path
            full_workdir_table['directory'].iloc[nwd]   = WD_path

        workdir_table.to_csv('working_dir.txt', sep=';')


    def remove_wd(self,nwd):
        global sz_w
        global sz_h
        global workdir_table
        global full_workdir_table

        _translate = QtCore.QCoreApplication.translate

        os.chdir(root_path)

        if not self.rm_list[nwd]:
            #remove wd
            self.rm_list[nwd] = True
            self.pb_wd_name[nwd].setGeometry(QtCore.QRect(10000*sz_w, 10*sz_h, 190*sz_w, (25+35*nwd)*sz_h))
            self.pb_wd_edit_n[nwd].setGeometry(QtCore.QRect(20005*sz_w, 10*sz_h, 80*sz_w, (25+35*nwd)*sz_h))
            self.pb_wd_edit_p[nwd].setGeometry(QtCore.QRect(29000*sz_w, 10*sz_h, 80*sz_w, (25+35*nwd)*sz_h))
            self.pb_wd_rm[nwd].setText(_translate("MainWD", '(+)'))#chr(128465)))
        else:
            # cancel wd withdrawal
            self.rm_list[nwd] = False
            self.pb_wd_name[nwd].setGeometry(QtCore.QRect(10*sz_w, (10+35*nwd)*sz_h, 190*sz_w, 25*sz_h))
            self.pb_wd_edit_n[nwd].setGeometry(QtCore.QRect(205*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
            self.pb_wd_edit_p[nwd].setGeometry(QtCore.QRect(290*sz_w, (10+35*nwd)*sz_h, 80*sz_w, 25*sz_h))
            self.pb_wd_rm[nwd].setText(_translate("MainWD", '(-)'))#chr(128465)))

        new_wd = True

        for wd in range(len(self.rm_list)):
            if not self.rm_list[wd]:
                if new_wd:
                    WD_name, WD_path = full_workdir_table.iloc[wd,0], full_workdir_table.iloc[wd,1]
                    data = {'working_dir_name':[WD_name],
                            'directory': [WD_path]}
                    workdir_table = pd.DataFrame(data=data)
                    new_wd = False
                else:
                    workdir_table = workdir_table.append(full_workdir_table.iloc[wd])
        workdir_table.to_csv('working_dir.txt', sep=';')

    def add_wd(self):
        global sz_w
        global sz_h

        global workdir_table
        global full_workdir_table
        _translate = QtCore.QCoreApplication.translate

        os.system(python_path + ' ' + root_path + "wd_manage_gui.py True False " + root_path)
        fs = open('wd_manage.txt', 'r')
        txt = fs.readline()
        WD_name = txt.split(';')[0]
        os.system(python_path + ' ' + root_path + "wd_manage_gui.py False True " + root_path)
        fs = open('wd_manage.txt', 'r')
        txt = fs.readline()
        WD_path = txt.split(';')[1]

        os.chdir(WD_path)
        wc.write_conditions()
        os.chdir(WD_path)
        wk.write_kinetic()
        os.chdir(WD_path)
        wr.write_results()

        full_workdir_table = full_workdir_table.append(full_workdir_table.iloc[0])
        full_workdir_table['working_dir_name'].iloc[-1] = WD_name
        full_workdir_table['directory'].iloc[-1]        = WD_path

        workdir_table = workdir_table.append(full_workdir_table.iloc[-1])
#        workdir_table = workdir_table.append(full_workdir_table.iloc[0])
#        workdir_table['working_dir_name'].iloc[-1] = WD_name
#        workdir_table['directory'].iloc[-1]        = WD_path


        os.chdir(root_path)
        workdir_table.to_csv('working_dir.txt', sep=';')

        wdn = len(self.rm_list)

        self.pb_wd_name[wdn].setGeometry(QtCore.QRect(10*sz_w, (10+35*wdn)*sz_h, 190*sz_w, 25*sz_h))
        self.pb_wd_name[wdn].setText(_translate("MainWD", WD_name))
        self.pb_wd_edit_n[wdn].setGeometry(QtCore.QRect(205*sz_w, (10+35*wdn)*sz_h, 80*sz_w, 25*sz_h))
        self.pb_wd_edit_p[wdn].setGeometry(QtCore.QRect(290*sz_w, (10+35*wdn)*sz_h, 80*sz_w, 25*sz_h))
        self.pb_wd_rm[wdn].setGeometry(QtCore.QRect(375*sz_w, (10+35*wdn)*sz_h, 25*sz_w, 25*sz_h))
        self.pb_wd_rm[wdn].setText(_translate("MainWD", '(-)'))#chr(128465)))
        self.rm_list.append(False)
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0*sz_w, 0*sz_h, 400*sz_w, (10+35*(wdn+1))*sz_h))




class Files_windows(QtWidgets.QWidget):

    def __init__(self):
        global sz_w
        global sz_h
        super().__init__()
        self.title = 'PyQt5 file dialogs - pythonspot.com'
        self.left = 10
        self.top = 10
        self.width = 640*sz_w
        self.height = 480*sz_h

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.show()



#if __name__ == "__main__":
app = QtWidgets.QApplication(sys.argv)
MainWD = QtWidgets.QMainWindow()
ui = Ui_MainWD()
ui.setupUi(MainWD)
MainWD.show()
sys.exit(app.exec_())


#
#
#def get_wd_folder():
#    global sz_w
#    global sz_h
#    global WD_path
#
#    app = QtWidgets.QApplication(sys.argv)
#    Window_f = QtWidgets.QMainWindow()
#
#    _height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()
#    _width = QtWidgets.QDesktopWidget().screenGeometry(-1).width()
#    sz_w = _width/1366
#    sz_h = _height/768
#
#    Window_f.setObjectName("Brookesia")
#    Window_f.resize(790*sz_w, 671*sz_h)
#
#    A=Files_windows()
#    options = QtWidgets.QFileDialog.Options()
#    options |= QtWidgets.QFileDialog.DontUseNativeDialog
#    WD_path = QtWidgets.QFileDialog.getExistingDirectory(A, "Select the new working directory", options=options)
#    WD_name, ok = QtWidgets.QInputDialog.getText(A, 'Text Input Dialog', 'Enter the new name:')
#    return WD_path, WD_name



