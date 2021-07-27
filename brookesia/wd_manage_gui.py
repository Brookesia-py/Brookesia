#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 14:00:02 2020

@author: matynia
"""

from PyQt5 import QtCore, QtGui, QtWidgets

global sz_w
global sz_h
global WD_path
import sys
import time as timer
import os

WD_n      = sys.argv[1]
WD_p      = sys.argv[2]
root_path = sys.argv[3]
if os.name != 'nt': # for Linux or Mac
    pers_config_path = '~/.Brookesia'
    try:
        os.chdir(os.path.expanduser(pers_config_path))
    except:
        os.chdir(os.path.expanduser('~'))
        try:
            os.mkdir('.Brookesia')
        except:
            a=False
        os.chdir('.Brookesia')
        pers_config_path = os.getcwd()
else: # for windows
    pers_config_path = root_path


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


app = QtWidgets.QApplication(sys.argv)
Window_f = QtWidgets.QMainWindow()

_height = QtWidgets.QDesktopWidget().screenGeometry(-1).height()
_width = QtWidgets.QDesktopWidget().screenGeometry(-1).width()
sz_w = _width/1366
sz_h = _height/768

Window_f.setObjectName("Brookesia")
Window_f.resize(790*sz_w, 671*sz_h)

A=Files_windows()
options = QtWidgets.QFileDialog.Options()
options |= QtWidgets.QFileDialog.DontUseNativeDialog

WD_path = sys.argv[1]
WD_name = sys.argv[2]

WD_path = 'False' ; WD_name = 'False'
if WD_p=='True':
    try:    os.chdir(pers_config_path)
    except: os.chdir(os.path.expanduser(pers_config_path))
    fs = open('wd_activ.txt', 'r')
    txt = fs.readline()
    previous_WD_path = txt.split(';')[1]
    try:
        os.chdir(previous_WD_path)
        os.chdir('../')
    except:
        a = False
    WD_path = QtWidgets.QFileDialog.getExistingDirectory(A, "Select the new working directory", options=options)

if WD_n=='True':
    WD_name, ok = QtWidgets.QInputDialog.getText(A, 'Text Input Dialog', 'Enter the new name:')

try:    os.chdir(pers_config_path)
except: os.chdir(os.path.expanduser(pers_config_path))
fd = open('wd_manage.txt', 'w')
fd.write(WD_name + ';' + WD_path)
timer.sleep(0.2)

