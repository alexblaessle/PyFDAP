#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller, and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFDAP.

#PyFDAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program. If not, see <http://www.gnu.org/licenses/>.

#=====================================================================================================================================
#Module Description
#=====================================================================================================================================

#Module containing PyQT classes needed for PyFRAP GUI including:
#1) molecule_dialog: Dialog for editing molecule object
#2) dataset_dialog: Dialog for editing embryo object
#3) bkgd_dataset_dialog: Dialog for editing bkgd object
#4) pre_dataset_dialog: Dialog for editing pre object
#5) noise_dataset_dialog: Dialog for editing noise object
#6) fit_dialog: Dialog for editing fit object
#7) mult_fit_dialog: Dialog for editing multiple fit objects at once
#8) about_dialog: About dialog of the software
#9) analyze_all_prog: Dialog shown during analysis of molecule
#10) analyze_all_thread: QThread for molecule analysis
#11) analyze_prog: Dialog shown during analysis of embryo
#12) analyze_thread: QThread for embryo analysis
#13) analyze_bkgd_prog: Dialog shown during analysis of all bkgds
#14) analyze_bkgd_thread: QThread for bkgd analysis
#15) fitting_prog: Dialog shown during fitting process
#16) fitting_thread: QThread for fitting of single fit
#17) fitting_all_thread: QThread for Fitting all three data series
#18) fitting_mol_thread: QThread for fitting of complete molecule
#19) select_fits: Dialog to select fits out of list of fits
#20) select_ignored_frames: Dialog to select frames to be ignored for fitting

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

import sys
from numpy import *
from PyQt4 import QtGui, QtCore
import pyfdap_img_module as pyfdap_img
import pyfdap_misc_module as pyfdap_misc
import pyfdap_stats_module as pyfdap_stats
import pyfdap_fit_module as pyfdap_fit
from pyfdap_term import *
import code
from embryo import *
import time
import os, os.path
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as  NavigationToolbar
from matplotlib.figure import Figure
#import skimage.io as skiio
import matplotlib.image as mpimg

#===================================================================================================================================
#Dialog for select/edit molecule
#===================================================================================================================================

class molecule_dialog(QtGui.QDialog):
	def __init__(self,molecule,parent):
		super(molecule_dialog,self).__init__(parent)
		
		self.molecule=molecule
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		self.lbl_name = QtGui.QLabel("Name:", self)
		
		self.qle_name = QtGui.QLineEdit(self.molecule.name)
		self.qle_name.editingFinished.connect(self.set_name)
		
		grid = QtGui.QGridLayout()
		grid.addWidget(self.lbl_name,1,1)
		grid.addWidget(self.qle_name,2,1)
		grid.addWidget(self.btn_done,2,2)
		
		self.setLayout(grid)    
			
		self.setWindowTitle('Edit Molecule')    
		self.show()
		
	def set_name(self):
		text=self.qle_name.text()
		self.molecule.name=str(text)	
	
	def done_pressed(self):
	
		self.done(1)
		return self.molecule
	
#===================================================================================================================================
#Dialog for select/edit data set
#===================================================================================================================================

class dataset_dialog(QtGui.QDialog):
	def __init__(self,embryo,parent):
		super(dataset_dialog,self).__init__(parent)
			
		self.embryo=embryo
					
		self.dpi = 100
		self.setMinimumSize(1000,500) 
		self.resize(1300,500)
		#Some variables
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_set_datafolder=QtGui.QPushButton('Change')
		self.btn_set_maskfolder=QtGui.QPushButton('Change')
		
		
		self.btn_next=QtGui.QPushButton('Next')
		self.btn_prev=QtGui.QPushButton('Previous')
		
		self.btn_set_for_following=QtGui.QPushButton('Copy geometry for following images')
		
		#Button Actions
		self.btn_set_datafolder.connect(self.btn_set_datafolder, QtCore.SIGNAL('clicked()'), self.sel_datafolder)
		self.btn_set_maskfolder.connect(self.btn_set_maskfolder, QtCore.SIGNAL('clicked()'), self.sel_maskfolder)
		
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		self.btn_next.connect(self.btn_next, QtCore.SIGNAL('clicked()'), self.next_img)
		self.btn_prev.connect(self.btn_prev, QtCore.SIGNAL('clicked()'), self.prev_img)
		self.btn_set_for_following.connect(self.btn_set_for_following, QtCore.SIGNAL('clicked()'), self.set_follow)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_name = QtGui.QLabel("Name:", self)
		self.lbl_name_datafolder = QtGui.QLabel("Photoconverted Folder:", self)
		self.lbl_name_maskfolder = QtGui.QLabel("Counterlabeled Folder:", self)
		
		self.lbl_name_datafolder.setToolTip("red images")
		self.lbl_name_maskfolder.setToolTip("green images")
		
		self.lbl_name_ft = QtGui.QLabel("Filetype:", self)
		self.lbl_name_enc = QtGui.QLabel("Encoding:", self)
		self.lbl_name_res = QtGui.QLabel("Image dimension (pixels):", self)
		self.lbl_name_radius = QtGui.QLabel("Embryo radius (pixels):", self)
		self.lbl_name_center = QtGui.QLabel("Embryo center (pixels):", self)
		self.lbl_name_framerate = QtGui.QLabel("Interval between images (s):", self)
		self.lbl_name_nframes = QtGui.QLabel("Number of frames:", self)
		self.lbl_name_tstart = QtGui.QLabel("tstart (s):", self)
		self.lbl_name_tend = QtGui.QLabel("tend (s):", self)
		self.lbl_name_delay = QtGui.QLabel("post delay (s):", self)
		
		self.lbl_datafolder = QtGui.QLabel(self.embryo.fn_datafolder, self)
		self.lbl_maskfolder = QtGui.QLabel(self.embryo.fn_maskfolder, self)
		
		if len(self.embryo.fn_datafolder)>50:
			self.lbl_datafolder.setText("..."+self.embryo.fn_datafolder[-50:])
		else:
			self.lbl_datafolder.setText(self.embryo.fn_datafolder)
		
		if len(self.embryo.fn_maskfolder)>50:
			self.lbl_maskfolder.setText("..."+self.embryo.fn_maskfolder[-50:])
		else:
			self.lbl_maskfolder.setText(self.embryo.fn_maskfolder)
			
		self.lbl_tend = QtGui.QLabel(str(self.embryo.tend), self)
		self.lbl_nframes = QtGui.QLabel(str(self.embryo.nframes), self)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_ft = QtGui.QComboBox(self)
		self.combo_ft.addItem("tif")
		
		self.combo_enc = QtGui.QComboBox(self)
		self.combo_enc.addItem("8bit")
		self.combo_enc.addItem("16bit")
		self.combo_enc.setCurrentIndex(1) 
		
		self.combo_ft.activated[str].connect(self.sel_ft)   
		self.combo_enc.activated[str].connect(self.sel_enc)   
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_name = QtGui.QLineEdit(self.embryo.name)
		
		self.qle_framerate = QtGui.QLineEdit(str(self.embryo.framerate))
		self.qle_tstart = QtGui.QLineEdit(str(self.embryo.tstart))
		self.qle_delay = QtGui.QLineEdit(str(self.embryo.post_delay))
		self.qle_res = QtGui.QLineEdit(str(self.embryo.data_res_px))
		
		centergrid= QtGui.QGridLayout()
		if shape(self.embryo.centers_embr_px)[0]>0:
			self.qle_radius = QtGui.QLineEdit(str(self.embryo.radiuses_embr_px[0]))
			self.qle_center_x = QtGui.QLineEdit(str(self.embryo.centers_embr_px[0][0]))
			self.qle_center_y = QtGui.QLineEdit(str(self.embryo.centers_embr_px[0][1]))
		else:
			self.qle_radius = QtGui.QLineEdit("")
			self.qle_center_x = QtGui.QLineEdit("")
			self.qle_center_y = QtGui.QLineEdit("")
			
		centergrid.addWidget(self.qle_center_x,0,0)
		centergrid.addWidget(self.qle_center_y,0,1)
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_framerate.setValidator(self.double_valid)
		self.qle_tstart.setValidator(self.double_valid)
		self.qle_res.setValidator(self.double_valid)
		self.qle_radius.setValidator(self.double_valid)
		self.qle_center_x.setValidator(self.double_valid)
		self.qle_center_y.setValidator(self.double_valid)
		
		self.qle_name.editingFinished.connect(self.set_name)
		self.qle_framerate.editingFinished.connect(self.set_framerate)
		self.qle_tstart.editingFinished.connect(self.set_tstart)
		self.qle_delay.editingFinished.connect(self.set_delay)
		self.qle_res.editingFinished.connect(self.set_res)
		self.qle_radius.editingFinished.connect(self.set_radius)
		self.qle_center_x.editingFinished.connect(self.set_center_x)
		self.qle_center_y.editingFinished.connect(self.set_center_y)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_rindicator = QtGui.QCheckBox('Range indicator', self)
		
		self.cb_rindicator.setCheckState(0)	
		
		self.connect(self.cb_rindicator, QtCore.SIGNAL('stateChanged(int)'), self.check_rindicator)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Plot frame
		#-------------------------------------------------------------------------------------------------------------------
		
		self.plot_frame = QtGui.QWidget()
		self.plot_frame.setMaximumWidth(1)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_name,1,1)
		
		grid.addWidget(self.lbl_name_datafolder,4,1)
		grid.addWidget(self.lbl_name_maskfolder,5,1)
		
		grid.addWidget(self.lbl_name_ft,6,1)
		grid.addWidget(self.lbl_name_enc,7,1)
		grid.addWidget(self.lbl_name_framerate,8,1)
		grid.addWidget(self.lbl_name_nframes,9,1)
		grid.addWidget(self.lbl_name_tstart,10,1)
		grid.addWidget(self.lbl_name_tend,11,1)
		grid.addWidget(self.lbl_name_delay,12,1)
		grid.addWidget(self.lbl_name_res,13,1)
		grid.addWidget(self.lbl_name_radius,14,1)
		grid.addWidget(self.lbl_name_center,15,1)
		grid.addWidget(self.btn_done,20,1)
		
		grid.addWidget(self.cb_rindicator,19,2)
		
		grid.addWidget(self.qle_name,1,2)
		
		grid.addWidget(self.lbl_datafolder,4,2)
		grid.addWidget(self.lbl_maskfolder,5,2)
		grid.addWidget(self.combo_ft,6,2)
		grid.addWidget(self.combo_enc,7,2)
		grid.addWidget(self.qle_framerate,8,2)
		grid.addWidget(self.lbl_nframes,9,2)
		grid.addWidget(self.qle_tstart,10,2)
		grid.addWidget(self.lbl_tend,11,2)
		grid.addWidget(self.qle_delay,12,2)
		grid.addWidget(self.qle_res,13,2)
		grid.addWidget(self.qle_radius,14,2)
		grid.addLayout(centergrid,15,2)
		grid.addWidget(self.btn_set_for_following,18,2)
		
		grid.addWidget(self.btn_set_datafolder,4,3)
		grid.addWidget(self.btn_set_maskfolder,5,3)
		
		grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Create Canvas
		#-------------------------------------------------------------------------------------------------------------------
		
		self.create_canvas()
		
		#-------------------------------------------------------------------------------------------------------------------
		#Radius/Center List
		#-------------------------------------------------------------------------------------------------------------------
		
		self.prop_list=QtGui.QTreeWidget()
		self.prop_list.setHeaderLabels(["Img","Center","Radius"])
	
		self.prop_list.setColumnWidth(0,40)
		self.prop_list.setColumnWidth(1,100)
		self.prop_list.setColumnWidth(2,75)
		self.prop_list.itemClicked.connect(self.prop_list_click)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.hbox2 = QtGui.QHBoxLayout()
		self.hbox2.addWidget(self.btn_prev)
		self.hbox2.addWidget(self.btn_next)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.canvas)
		self.vbox.addLayout(self.hbox2)
		
		#Add everything to horizontal Box
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(grid)
		self.hbox.addWidget(self.prop_list)
		self.hbox.addLayout(self.vbox)
		self.setLayout(self.hbox)    
		
		#Check if there is already a first image
		self.update_prop_list()
		if shape(self.embryo.radiuses_embr_px)[0]>0:
			self.show_img(0)	
			self.draw_patches(0)
			
		self.setWindowTitle('Edit Dataset')    
		self.show()
	
	def create_frame(self,frame):
		
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		frame.setBackgroundRole(QtGui.QPalette.Light)        
		frame.setLineWidth (1)
		frame.setFrameShadow (frame.Sunken)
		frame.frameRect().setWidth(100)
		
		return frame
		
	def sel_datafolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Data Directory"))
		if folder=='':
			return
		self.embryo.fn_datafolder=folder
	
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
		self.embryo.radiuses_embr_px=[(self.embryo.data_res_px/2)]*shape(self.file_list)[0]
		self.embryo.centers_embr_px=[[self.embryo.data_res_px/2,self.embryo.data_res_px/2]]*shape(self.file_list)[0]
		
		self.init_prop_list()
		
		#If folder string ends without /, add /
		if self.embryo.fn_datafolder[-1]=="/":
			pass
		else:
			self.embryo.fn_datafolder=self.embryo.fn_datafolder+"/"
		
		if len(self.embryo.fn_datafolder)>50:
			self.lbl_datafolder.setText("..."+self.embryo.fn_datafolder[-50:])
		else:
			self.lbl_datafolder.setText(self.embryo.fn_datafolder)
		
		self.embryo.nframes=int(shape(self.file_list)[0])
		self.update_tvec()
		
		self.lbl_nframes.setText(str(self.embryo.nframes))
		self.lbl_tend.setText(str(self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)))
		self.show_img(0)
		self.update_qles(0)
		
		self.sel_other_folder(folder,"green")
			
			
	def sel_maskfolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Mask Data Directory"))
		if folder=='':
			return
		self.embryo.fn_maskfolder = folder
		
		#If folder string ends without /, add /
		if self.embryo.fn_maskfolder[-1]=="/":
			pass
		else:
			self.embryo.fn_maskfolder=self.embryo.fn_maskfolder+"/"
		
		if len(self.embryo.fn_maskfolder)>50:
			self.lbl_maskfolder.setText("..."+self.embryo.fn_maskfolder[-50:])
		else:
			self.lbl_maskfolder.setText(self.embryo.fn_maskfolder)
		
		self.sel_other_folder(folder,"red")
		
	def sel_other_folder(self,folder,color):
		
		#Delete last "/" for os.path.split to work
		if folder[-1]=="/":
			folder=folder[:-1]
		
		#Getting upper folder
		[upper_folder,last_folder]=os.path.split(folder)
		
		#Check if other folder is where it should be, select if exists and is not empty
		otherfolder=upper_folder+"/"+color
		if os.path.isdir(otherfolder):
			
			#Check of folder is empty
			if pyfdap_misc.check_folder_empty(otherfolder):
				pass
				#print "Did not automatically select ", otherfolder, "as folder. Folder is empty!"
			else: 	
				#Add / so analysis functions don't get confused
				otherfolder=otherfolder+"/"
				
				if color=="green":
					#Update maskfolder
					self.embryo.fn_maskfolder = otherfolder
					if len(self.embryo.fn_maskfolder)>50:
						self.lbl_maskfolder.setText("..."+self.embryo.fn_maskfolder[-50:])
					else:
						self.lbl_maskfolder.setText(self.embryo.fn_maskfolder)
				
				elif color=="red":
					#Check if otherfolder is different from datafolder so we don't override any settings
					if self.embryo.fn_datafolder!=otherfolder:
						
						self.embryo.fn_datafolder=otherfolder
						
						self.file_list = pyfdap_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
						self.embryo.radiuses_embr_px=[(self.embryo.data_res_px/2)]*shape(self.file_list)[0]
						self.embryo.centers_embr_px=[[self.embryo.data_res_px/2,self.embryo.data_res_px/2]]*shape(self.file_list)[0]
						
						self.init_prop_list()
						
						if len(self.embryo.fn_datafolder)>50:
							self.lbl_datafolder.setText("..."+self.embryo.fn_datafolder[-50:])
						else:
							self.lbl_datafolder.setText(self.embryo.fn_datafolder)
					
						self.embryo.nframes=int(shape(self.file_list)[0])
						self.update_tvec()
						
						self.lbl_nframes.setText(str(self.embryo.nframes))
						self.lbl_tend.setText(str(self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)))
						self.show_img(0)
						self.update_qles(0)
		else:
			pass
			#print "Did not automatically select ", otherfolder, "as folder. Folder does not exists!"
		
	def sel_ft(self,text):
		if text=="tif":
			self.embryo.data_ft=text
		
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
		self.embryo.nframes=int(shape(self.file_list[0]))
		self.update_tvec()
		
		self.lbl_nframes.setText(str(self.embryo.nframes))
		self.lbl_tend.setText(str(self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)))
		self.show_img(self.curr_img_ind)
		
	def sel_enc(self,text):
		if text=="8bit":
			self.embryo.data_enc="uint8"
		elif text=="16bit":	
			self.embryo.data_enc="uint16"
		self.show_img(self.curr_img_ind)
	
	def set_name(self):
		text=self.qle_name.text()
		self.embryo.name=str(text)
	
	def set_framerate(self):
		text=self.qle_framerate.text()
		self.embryo.framerate=float(str(text))
		self.lbl_tend.setText(str(self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)+self.embryo.post_delay))
		self.update_tvec()
		
	def set_tstart(self):
		text=self.qle_tstart.text()
		self.embryo.tstart=float(str(text))
		self.lbl_tend.setText(str(self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)))
		self.update_tvec()
	
	def set_delay(self):
		text=self.qle_delay.text()
		self.embryo.post_delay=float(str(text))
		self.lbl_tend.setText(str(self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)+self.embryo.post_delay))
		self.update_tvec()
	
	def set_res(self):
		text=self.qle_res.text()
		self.embryo.data_res_px=float(str(text))
		self.ax.set_xlim([0, self.embryo.data_res_px])
		self.ax.set_ylim([0, self.embryo.data_res_px])
		
	def set_radius(self):
		
		text=self.qle_radius.text()
		self.embryo.radiuses_embr_px[self.curr_img_ind]=float(str(text))	
		self.update_prop_list_entry(self.curr_img_ind)
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.embryo.centers_embr_px[self.curr_img_ind][0],self.embryo.centers_embr_px[self.curr_img_ind][1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.embryo.centers_embr_px[self.curr_img_ind][0])
		
		embr_circ=ptc.Circle(self.embryo.centers_embr_px[self.curr_img_ind],radius=self.embryo.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
		self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
		self.xcoords2.append(self.embryo.centers_embr_px[self.curr_img_ind][0])
	
		self.fig.canvas.draw()	
		
	def set_center_x(self):
		text=self.qle_center_x.text()
		self.embryo.centers_embr_px[self.curr_img_ind][0]=float(str(text))
		self.update_prop_list_entry(self.curr_img_ind)
		
		if hasattr(self,'embr_circ_in_ax'):
			had_circ=1
		else:
			had_circ=0
			
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.embryo.centers_embr_px[self.curr_img_ind][0],self.embryo.centers_embr_px[self.curr_img_ind][1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.embryo.centers_embr_px[self.curr_img_ind][0])
		
		if had_circ==1:
			embr_circ=ptc.Circle(self.embryo.centers_embr_px[self.curr_img_ind],radius=self.embryo.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(self.embryo.centers_embr_px[self.curr_img_ind][0])
		
		self.fig.canvas.draw()
		
	def set_center_y(self):
		text=self.qle_center_y.text()
		self.embryo.centers_embr_px[self.curr_img_ind][1]=float(str(text))
		self.update_prop_list_entry(self.curr_img_ind)
		
		if hasattr(self,'embr_circ_in_ax'):
			had_circ=1
		else:
			had_circ=0
			
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.embryo.centers_embr_px[self.curr_img_ind][0],self.embryo.centers_embr_px[self.curr_img_ind][1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.embryo.centers_embr_px[self.curr_img_ind][0])
		
		if had_circ==1:
			embr_circ=ptc.Circle(self.embryo.centers_embr_px[self.curr_img_ind],radius=self.embryo.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(self.embryo.centers_embr_px[self.curr_img_ind][0])
		
		self.fig.canvas.draw()
	
	def set_follow(self):
		reply = QtGui.QMessageBox.question(self, 'Message',"Are you sure that you want to copy the geometry properties to all following images?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
	
		if reply == QtGui.QMessageBox.Yes:
			for i in range(shape(self.file_list)[0]):
				if i>=self.curr_img_ind:
					
					self.embryo.radiuses_embr_px[i]=self.embryo.radiuses_embr_px[self.curr_img_ind]
					self.embryo.centers_embr_px[i]=self.embryo.centers_embr_px[self.curr_img_ind]
					self.update_prop_list_entry(i)
		else:
			return
	
	def prop_list_click(self):
		
		self.rem_rindicator()
		
		self.new_img_ind=self.prop_list.currentIndex().row()
		self.update_qles(self.new_img_ind)
		self.draw_patches(self.new_img_ind)
		self.show_img(self.new_img_ind)
	
	def next_img(self):
		if self.curr_img_ind==shape(self.file_list)[0]-1:
			self.prop_list.setCurrentItem(self.nodes[0])
			self.update_qles(0)
			self.draw_patches(0)
			self.show_img(0)
		else:	
			self.prop_list.setCurrentItem(self.nodes[self.curr_img_ind+1])
			self.update_qles(self.curr_img_ind+1)
			self.draw_patches(self.curr_img_ind+1)
			self.show_img(self.curr_img_ind+1)
					
	def prev_img(self):
		if self.curr_img_ind==0:
			self.prop_list.setCurrentItem(self.nodes[shape(self.file_list)[0]-1])
			self.update_qles(shape(self.file_list)[0]-1)
			self.draw_patches(shape(self.file_list)[0]-1)
			self.show_img(shape(self.file_list)[0]-1)
				
		else:	
			self.prop_list.setCurrentItem(self.nodes[self.curr_img_ind-1])
			self.update_qles(self.curr_img_ind-1)
			self.draw_patches(self.curr_img_ind-1)
			self.show_img(self.curr_img_ind-1)
	
	def update_tvec(self):
		
		self.embryo.tend=self.embryo.tstart+self.embryo.framerate*(self.embryo.nframes-1)
		self.embryo.steps_data=self.embryo.nframes
		self.embryo.tvec_data=linspace(self.embryo.tstart,self.embryo.tend,self.embryo.nframes)
		self.embryo.tvec_data[1:]=self.embryo.tvec_data[1:]+self.embryo.post_delay
		self.embryo.tend=self.embryo.tend+self.embryo.post_delay
		
		
	def update_qles(self,ind):
		
		self.qle_center_x.setText(str(self.embryo.centers_embr_px[ind][0]))
		self.qle_center_y.setText(str(self.embryo.centers_embr_px[ind][1]))
		self.qle_radius.setText(str(self.embryo.radiuses_embr_px[ind]))
	
	def update_prop_list_entry(self,ind):
		
		self.nodes[ind].setText(1,str(self.embryo.centers_embr_px[ind]))
		self.nodes[ind].setText(2,str(self.embryo.radiuses_embr_px[ind]))
	
	def init_prop_list(self):
		
		self.prop_list.clear()
		self.update_prop_list()
	
	def update_prop_list(self):
		
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
		self.nodes=[]
		
		if len(self.file_list)==len(self.embryo.radiuses_embr_px):
		
			for i in range(shape(self.file_list)[0]):
				
				self.nodes.append(QtGui.QTreeWidgetItem(self.prop_list,[str(i),str(self.embryo.centers_embr_px[i]),str(self.embryo.radiuses_embr_px[i])]))
			
			if shape(self.file_list)[0]>0:
				self.prop_list.setCurrentItem(self.nodes[0])
		
		else:
			pass
		
	def create_canvas(self):
			
		h=500/self.dpi
		v=500/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.plot_frame)
		
		self.ax = self.fig.add_subplot(111)
		self.ax.set_xlim([0, self.embryo.data_res_px])
		self.ax.set_ylim([0, self.embryo.data_res_px])
		
		#Connect to mouseclick event
		self.canvas.mpl_connect('button_press_event', self.get_mouse_embr_radius)
		
		self.canvas.draw()
		return 
	
	def draw_patches(self,ind):
		
		self.clear_patch_canvas()
		
		if shape(self.embryo.centers_embr_px[ind])[0]>0:
			pt=ptc.Circle([self.embryo.centers_embr_px[ind][0],self.embryo.centers_embr_px[ind][1]],radius=3,fill=True,color='r')
			self.pts_in_ax2.append(self.ax.add_patch(pt))
			self.xcoords2.append(self.embryo.centers_embr_px[ind])
			self.fig.canvas.draw()	
		if shape(self.embryo.radiuses_embr_px)[0]>0:
			embr_circ=ptc.Circle(self.embryo.centers_embr_px[ind],radius=self.embryo.radiuses_embr_px[ind],fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(100)
			self.fig.canvas.draw()		
				
	def show_img(self,ind):
		
		self.file_list= pyfdap_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
		
		#Check if there is a first image
		if shape(self.file_list)[0]>ind-1 and shape(self.file_list)[0]>0:
		
			#Grab first picture
			fn=self.embryo.fn_datafolder+self.file_list[ind]
			
			#Check if file exists
			if os.path.isfile(fn)==False:
				return 
			
			#Load img
			data_img = mpimg.imread(fn).astype(self.embryo.data_enc)
			data_vals=data_img.real
			data_vals=data_vals.astype('float')
			
			#Flip to be coherent with Fiji
			data_vals=flipud(data_vals)
			
			#Plot img
			self.ax.imshow(data_vals)
			self.canvas.draw()
			
			#Grab resolution
			self.embryo.data_res_px=shape(data_vals)[0]
			self.qle_res.setText(str(self.embryo.data_res_px))
			
			self.curr_img_ind=ind
			
			if self.cb_rindicator.isChecked():
				self.rindicator(data_vals)	
			
		else:
			#Clear if there is no first image
			self.ax.cla()
	
	def get_mouse_embr_radius(self,event):
		
		#Left click to define center and then radius
		if event.button==1:
			
			#Check if clicked withing axes
			if event.xdata==None:
				return
			
			#Check if there is already a center
			if shape(self.xcoords2)[0]<1:
					
				#No center yet
				self.xcoords2.append(event.xdata)
				self.ycoords2.append(event.ydata)
				
				#Plot center
				pt=ptc.Circle([event.xdata,event.ydata],radius=3,fill=True,color='r')
				self.pts_in_ax2.append(self.ax.add_patch(pt))
				
				#Set center_embr
				self.embryo.centers_embr_px[self.curr_img_ind]=[event.xdata,event.ydata]
				
				self.qle_center_x.setText(str(self.embryo.centers_embr_px[self.curr_img_ind][0]))
				self.qle_center_y.setText(str(self.embryo.centers_embr_px[self.curr_img_ind][1]))
				self.prop_list.currentItem().setText(1,str(self.embryo.centers_embr_px[self.curr_img_ind]))
				
			elif shape(self.xcoords2)[0]==1:
				
				#Center yet, so this click is for radius
				self.xcoords2.append(event.xdata)
				self.ycoords2.append(event.ydata)
					
				#Set radius_embr
				self.embryo.radiuses_embr_px[self.curr_img_ind]=sqrt((event.xdata-self.embryo.centers_embr_px[self.curr_img_ind][0])**2+(event.ydata-self.embryo.centers_embr_px[self.curr_img_ind][1])**2)
				
				#Plot circle
				embr_circ=ptc.Circle(self.embryo.centers_embr_px[self.curr_img_ind],radius=self.embryo.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
				self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
				
				self.qle_radius.setText(str(self.embryo.radiuses_embr_px[self.curr_img_ind]))
				self.prop_list.currentItem().setText(2,str(self.embryo.radiuses_embr_px[self.curr_img_ind]))
			
			else:
			
				self.clear_patch_canvas()
				
		self.fig.canvas.draw()			
	
	def clear_patch_canvas(self):
	
		#Clearing Canvas
		
		#Removing all points
		for pt in self.pts_in_ax2:
			pt.remove()
		
		#Removing circle
		if hasattr(self,'embr_circ_in_ax'):
			self.embr_circ_in_ax.remove()
			del self.embr_circ_in_ax

		#Setting back arrays
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		self.fig.canvas.draw()
	
	def check_rindicator(self,value,ind=None):
		
		#Get current index 
		if ind==None:
			ind=self.prop_list.currentIndex().row()
		
		#If checked, display rindicator
		if value==2:
		
			#Get file list
			self.file_list= pyfdap_misc.get_sorted_folder_list(self.embryo.fn_datafolder,self.embryo.data_ft)
		
			#Check if there is a first image
			if shape(self.file_list)[0]>ind-1 and shape(self.file_list)[0]>0:
			
				#Grab first picture
				fn=self.embryo.fn_datafolder+self.file_list[ind]
				
				#Check if file exists
				if os.path.isfile(fn)==False:
					return 
				
				#Load img
				data_img = mpimg.imread(fn).astype(self.embryo.data_enc)
				data_vals=data_img.real
				data_vals=data_vals.astype('float')
				
				#Flip to be coherent with Fiji
				data_vals=flipud(data_vals)
			
				self.rindicator(data_vals)
		
		#If unchecked remove rindicator
		else:	
			self.rem_rindicator()
					
	def rindicator(self,img):
		
		#Default value
		maxint=65535
		
		#Check what is maximum intensity
		if self.embryo.data_enc=='16bit':
			maxint=65535
		if self.embryo.data_enc=='8bit':
			maxint=255
			
		#Get indices of maximum intensity 
		max_pxs=pyfdap_img.get_max_int_pxs(img,maxint)
		
		#Make images with pxs filled
		labeled_image=label_pxs(img,max_pxs,fill=nan)
		
		#Plot into current canvas
		self.imgplt=self.ax.imshow(labeled_image)
		self.imgplt.set_cmap('Greys')
		self.canvas.draw()
		
	def rem_rindicator(self):
		if hasattr(self,'imgplt') and self.imgplt!=None:
			self.imgplt.remove()
			self.canvas.draw()
			self.imgplt=None
		
	def done_pressed(self):
	
		self.done(1)
		return self.embryo

#===================================================================================================================================
#Dialog for select/edit background data set
#===================================================================================================================================

class bkgd_dataset_dialog(QtGui.QDialog):
	def __init__(self,bkgd,parent):
		super(bkgd_dataset_dialog,self).__init__(parent)
		
		self.bkgd=bkgd
			
		self.dpi = 100
		self.setMinimumSize(1000,500) 
		self.resize(1300,500)
		#Some variables
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_set_bkgdfolder=QtGui.QPushButton('Change')
		self.btn_set_maskfolder=QtGui.QPushButton('Change')
		self.btn_next=QtGui.QPushButton('Next')
		self.btn_prev=QtGui.QPushButton('Previous')
		
		self.btn_set_for_following=QtGui.QPushButton('Copy geometry for following images')
		
		#Button Actions
		self.btn_set_bkgdfolder.connect(self.btn_set_bkgdfolder, QtCore.SIGNAL('clicked()'), self.sel_bkgdfolder)
		self.btn_set_maskfolder.connect(self.btn_set_maskfolder, QtCore.SIGNAL('clicked()'), self.sel_maskfolder)
		
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		self.btn_next.connect(self.btn_next, QtCore.SIGNAL('clicked()'), self.next_img)
		self.btn_prev.connect(self.btn_prev, QtCore.SIGNAL('clicked()'), self.prev_img)
		self.btn_set_for_following.connect(self.btn_set_for_following, QtCore.SIGNAL('clicked()'), self.set_follow)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_name = QtGui.QLabel("Name:", self)
		self.lbl_name_bkgdfolder = QtGui.QLabel("Photoconverted Folder:", self)
		self.lbl_name_maskfolder = QtGui.QLabel("Counterlabeled Folder:", self)
		
		self.lbl_name_bkgdfolder.setToolTip("red images")
		self.lbl_name_maskfolder.setToolTip("green images")
		self.lbl_name_ft = QtGui.QLabel("Filetype:", self)
		self.lbl_name_enc = QtGui.QLabel("Encoding:", self)
		self.lbl_name_res = QtGui.QLabel("RImage dimensions (pixels):", self)
		self.lbl_name_radius = QtGui.QLabel("Embryo radius (pixels):", self)
		self.lbl_name_center = QtGui.QLabel("Embryo center (pixels):", self)
		
		self.lbl_bkgdfolder = QtGui.QLabel(self.bkgd.fn_bkgdfolder, self)
		self.lbl_maskfolder = QtGui.QLabel(self.bkgd.fn_maskfolder, self)
		
		if len(self.bkgd.fn_bkgdfolder)>50:
			self.lbl_bkgdfolder.setText("..."+self.bkgd.fn_bkgdfolder[-50:])
		else:
			self.lbl_bkgdfolder.setText(self.bkgd.fn_bkgdfolder)
		
		if len(self.bkgd.fn_maskfolder)>50:
			self.lbl_maskfolder.setText("..."+self.bkgd.fn_maskfolder[-50:])
		else:
			self.lbl_maskfolder.setText(self.bkgd.fn_maskfolder)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_ft = QtGui.QComboBox(self)
		self.combo_ft.addItem("tif")
		
		self.combo_enc = QtGui.QComboBox(self)
		self.combo_enc.addItem("8bit")
		self.combo_enc.addItem("16bit")
		self.combo_enc.setCurrentIndex(1) 
		
		self.combo_ft.activated[str].connect(self.sel_ft)   
		self.combo_enc.activated[str].connect(self.sel_enc)   
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_name = QtGui.QLineEdit(self.bkgd.name)
		self.qle_res = QtGui.QLineEdit(str(self.bkgd.res_px))
		
		
		centergrid= QtGui.QGridLayout()
		if shape(self.bkgd.centers_embr_px)[0]>0:
			self.qle_radius = QtGui.QLineEdit(str(self.bkgd.radiuses_embr_px[0]))
			self.qle_center_x = QtGui.QLineEdit(str(self.bkgd.centers_embr_px[0][0]))
			self.qle_center_y = QtGui.QLineEdit(str(self.bkgd.centers_embr_px[0][1]))
		else:
			self.qle_radius = QtGui.QLineEdit("")
			self.qle_center_x = QtGui.QLineEdit("")
			self.qle_center_y = QtGui.QLineEdit("")
			
		centergrid.addWidget(self.qle_center_x,0,0)
		centergrid.addWidget(self.qle_center_y,0,1)
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_res.setValidator(self.double_valid)
		self.qle_radius.setValidator(self.double_valid)
		self.qle_center_x.setValidator(self.double_valid)
		self.qle_center_y.setValidator(self.double_valid)
		
		self.qle_name.editingFinished.connect(self.set_name)
		self.qle_res.editingFinished.connect(self.set_res)
		self.qle_radius.editingFinished.connect(self.set_radius)
		self.qle_center_x.editingFinished.connect(self.set_center_x)
		self.qle_center_y.editingFinished.connect(self.set_center_y)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Plot frame
		#-------------------------------------------------------------------------------------------------------------------
		
		self.plot_frame = QtGui.QWidget()
		self.plot_frame.setMaximumWidth(1)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_name,1,1)
		grid.addWidget(self.lbl_name_bkgdfolder,2,1)
		grid.addWidget(self.lbl_name_maskfolder,3,1)
		grid.addWidget(self.lbl_name_ft,4,1)
		grid.addWidget(self.lbl_name_enc,5,1)
		grid.addWidget(self.lbl_name_res,6,1)
		grid.addWidget(self.lbl_name_radius,7,1)
		grid.addWidget(self.lbl_name_center,8,1)
		grid.addWidget(self.btn_done,10,1)
		
		grid.addWidget(self.qle_name,1,2)
		grid.addWidget(self.lbl_bkgdfolder,2,2)
		grid.addWidget(self.lbl_maskfolder,3,2)
		grid.addWidget(self.combo_ft,4,2)
		grid.addWidget(self.combo_enc,5,2)
		grid.addWidget(self.qle_res,6,2)
		grid.addWidget(self.qle_radius,7,2)
		grid.addLayout(centergrid,8,2)
		grid.addWidget(self.btn_set_for_following,9,2)
		
		grid.addWidget(self.btn_set_bkgdfolder,2,3)
		grid.addWidget(self.btn_set_maskfolder,3,3)
		
		grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Create Canvas
		#-------------------------------------------------------------------------------------------------------------------
		
		self.create_canvas()
		
		#-------------------------------------------------------------------------------------------------------------------
		#Radius/Center List
		#-------------------------------------------------------------------------------------------------------------------
		
		self.prop_list=QtGui.QTreeWidget()
		self.prop_list.setHeaderLabels(["Img","Center","Radius"])
	
		self.prop_list.setColumnWidth(0,40)
		self.prop_list.setColumnWidth(1,100)
		self.prop_list.setColumnWidth(2,75)
		self.prop_list.itemClicked.connect(self.prop_list_click) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.hbox2 = QtGui.QHBoxLayout()
		self.hbox2.addWidget(self.btn_prev)
		self.hbox2.addWidget(self.btn_next)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.canvas)
		self.vbox.addLayout(self.hbox2)
		
		#Add everything to Horizontal Box
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(grid)
		self.hbox.addWidget(self.prop_list)
		self.hbox.addLayout(self.vbox)
		self.setLayout(self.hbox)    
		
		#Check if there is already a first image
		self.update_prop_list()
		self.show_img(0)	
		self.draw_patches(0)
		
		self.setWindowTitle('Edit Background Dataset')    
		self.show()
	
	def create_frame(self,frame):
		
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		frame.setBackgroundRole(QtGui.QPalette.Light)      
		frame.setLineWidth (1)
		frame.setFrameShadow (frame.Sunken)
		frame.frameRect().setWidth(100)
		
		return frame
		
	def sel_bkgdfolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Data Directory"))
		if folder=='':
			return
		self.bkgd.fn_bkgdfolder=folder
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.bkgd.fn_bkgdfolder,self.bkgd.data_ft)
		self.bkgd.centers_embr_px=[[self.bkgd.res_px/2,self.bkgd.res_px/2]]*shape(self.file_list)[0]
		self.bkgd.radiuses_embr_px=[253.]*shape(self.file_list)[0]
		
		self.init_prop_list()
		
		#If folder string ends without /, add /
		if self.bkgd.fn_bkgdfolder[-1]=="/":
			pass
		else:
			self.bkgd.fn_bkgdfolder=self.bkgd.fn_bkgdfolder+"/"
		
		if len(self.bkgd.fn_bkgdfolder)>50:
			self.lbl_bkgdfolder.setText("..."+self.bkgd.fn_bkgdfolder[-50:])
		else:
			self.lbl_bkgdfolder.setText(self.bkgd.fn_bkgdfolder)
		
		self.show_img(0)
		self.update_qles(0)
		
		self.sel_other_folder(folder,"green")
		
	def sel_maskfolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Mask Data Directory"))
		if folder=='':
			return
		self.bkgd.fn_maskfolder = folder
		
		#If folder string ends without /, add /
		if self.bkgd.fn_maskfolder[-1]=="/":
			pass
		else:
			self.bkgd.fn_maskfolder=self.bkgd.fn_maskfolder+"/"
		
		if len(self.bkgd.fn_maskfolder)>50:
			self.lbl_maskfolder.setText("..."+self.bkgd.fn_maskfolder[-50:])
		else:
			self.lbl_maskfolder.setText(self.bkgd.fn_maskfolder)
		
		
		self.sel_other_folder(folder,"red")
		
	def sel_other_folder(self,folder,color):
		
		#Delete last "/" for os.path.split to work
		if folder[-1]=="/":
			folder=folder[:-1]
		
		#Getting upper folder
		[upper_folder,last_folder]=os.path.split(folder)
		
		#Check if other folder is where it should be, select if exists and is not empty
		otherfolder=upper_folder+"/"+color
		if os.path.isdir(otherfolder):
			
			#Check of folder is empty
			if pyfdap_misc.check_folder_empty(otherfolder):
				pass
				#print "Did not automatically select ", otherfolder, "as folder. Folder is empty!"
			else: 	
				#Add / so analysis functions don't get confused
				otherfolder=otherfolder+"/"
				
				if color=="green":
					#Update maskfolder
					self.bkgd.fn_maskfolder = otherfolder
					if len(self.bkgd.fn_maskfolder)>50:
						self.lbl_maskfolder.setText("..."+self.bkgd.fn_maskfolder[-50:])
					else:
						self.lbl_maskfolder.setText(self.bkgd.fn_maskfolder)
				
				elif color=="red":
					#Check if otherfolder is different from datafolder so we don't override any settings
					if self.bkgd.fn_bkgdfolder!=otherfolder:
						
						self.bkgd.fn_bkgdfolder=otherfolder
						self.file_list = pyfdap_misc.get_sorted_folder_list(self.bkgd.fn_bkgdfolder,self.bkgd.data_ft)
						self.bkgd.centers_embr_px=[[self.bkgd.res_px/2,self.bkgd.res_px/2]]*shape(self.file_list)[0]
						self.bkgd.radiuses_embr_px=[253.]*shape(self.file_list)[0]
						
						self.init_prop_list()
						
						if len(self.bkgd.fn_bkgdfolder)>50:
							self.lbl_bkgdfolder.setText("..."+self.bkgd.fn_bkgdfolder[-50:])
						else:
							self.lbl_bkgdfolder.setText(self.bkgd.fn_bkgdfolder)
						
						self.show_img(0)
						self.update_qles(0)
		else:
			pass
			#print "Did not automatically select ", otherfolder, "as folder. Folder does not exists!"
	
	def sel_ft(self,text):
		if text=="tif":
			self.bkgd.data_ft=text
		
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.bkgd.fn_bkgdfolder,self.bkgd.data_ft)		
		self.show_img(self.curr_img_ind)
		
	def sel_enc(self,text):
		if text=="8bit":
			self.bkgd.data_enc="uint8"
		elif text=="16bit":	
			self.bkgd.data_enc="uint16"
		self.show_img(self.curr_img_ind)
	
	def set_name(self):
		text=self.qle_name.text()
		self.bkgd.name=str(text)
	
	def set_res(self):
		text=self.qle_res.text()
		self.bkgd.res_px=float(str(text))
		self.ax.set_xlim([0, self.bkgd.res_px])
		self.ax.set_ylim([0, self.bkgd.res_px])
		
	def set_radius(self):
		
		text=self.qle_radius.text()
		self.bkgd.radiuses_embr_px[self.curr_img_ind]=float(str(text))	
		self.update_prop_list_entry(self.curr_img_ind)
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.bkgd.centers_embr_px[self.curr_img_ind][0],self.bkgd.centers_embr_px[self.curr_img_ind][1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.bkgd.centers_embr_px[self.curr_img_ind][0])
		
		embr_circ=ptc.Circle(self.bkgd.centers_embr_px[self.curr_img_ind],radius=self.bkgd.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
		self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
		self.xcoords2.append(self.bkgd.centers_embr_px[self.curr_img_ind][0])
	
		self.fig.canvas.draw()	
		
	def set_center_x(self):
		text=self.qle_center_x.text()
		self.bkgd.centers_embr_px[self.curr_img_ind][0]=float(str(text))
		self.update_prop_list_entry(self.curr_img_ind)
		
		if hasattr(self,'embr_circ_in_ax'):
			had_circ=1
		else:
			had_circ=0
			
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.bkgd.centers_embr_px[self.curr_img_ind][0],self.bkgd.centers_embr_px[self.curr_img_ind][1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.bkgd.centers_embr_px[self.curr_img_ind][0])
		
		if had_circ==1:
			embr_circ=ptc.Circle(self.bkgd.centers_embr_px[self.curr_img_ind],radius=self.bkgd.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(self.bkgd.centers_embr_px[self.curr_img_ind][0])
		
		self.fig.canvas.draw()
		
	def set_center_y(self):
		text=self.qle_center_y.text()
		self.bkgd.centers_embr_px[self.curr_img_ind][1]=float(str(text))
		self.update_prop_list_entry(self.curr_img_ind)
		
		if hasattr(self,'embr_circ_in_ax'):
			had_circ=1
		else:
			had_circ=0
			
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.bkgd.centers_embr_px[self.curr_img_ind][0],self.bkgd.centers_embr_px[self.curr_img_ind][1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.bkgd.centers_embr_px[self.curr_img_ind][0])
		
		if had_circ==1:
			embr_circ=ptc.Circle(self.bkgd.centers_embr_px[self.curr_img_ind],radius=self.bkgd.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(self.bkgd.centers_embr_px[self.curr_img_ind][0])
		
		self.fig.canvas.draw()
	
	def set_follow(self):
		reply = QtGui.QMessageBox.question(self, 'Message',"Are you sure that you want to copy the geometry properties to all following images?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
	
		if reply == QtGui.QMessageBox.Yes:
			for i in range(shape(self.file_list)[0]):
				if i>=self.curr_img_ind:
					
					self.bkgd.radiuses_embr_px[i]=self.bkgd.radiuses_embr_px[self.curr_img_ind]
					self.bkgd.centers_embr_px[i]=self.bkgd.centers_embr_px[self.curr_img_ind]
					self.update_prop_list_entry(i)
		else:
			return
	
	def prop_list_click(self):
		
		self.new_img_ind=self.prop_list.currentIndex().row()
		self.update_qles(self.new_img_ind)
		self.draw_patches(self.new_img_ind)
		self.show_img(self.new_img_ind)
		
	
	def next_img(self):
		if self.curr_img_ind==shape(self.file_list)[0]-1:
			self.prop_list.setCurrentItem(self.nodes[0])
			self.update_qles(0)
			self.draw_patches(0)
			self.show_img(0)
		else:	
			self.prop_list.setCurrentItem(self.nodes[self.curr_img_ind+1])
			self.update_qles(self.curr_img_ind+1)
			self.draw_patches(self.curr_img_ind+1)
			self.show_img(self.curr_img_ind+1)
					
	def prev_img(self):
		if self.curr_img_ind==0:
			self.prop_list.setCurrentItem(self.nodes[shape(self.file_list)[0]-1])
			self.update_qles(shape(self.file_list)[0]-1)
			self.draw_patches(shape(self.file_list)[0]-1)
			self.show_img(shape(self.file_list)[0]-1)
				
		else:	
			self.prop_list.setCurrentItem(self.nodes[self.curr_img_ind-1])
			self.update_qles(self.curr_img_ind-1)
			self.draw_patches(self.curr_img_ind-1)
			self.show_img(self.curr_img_ind-1)
	
	def update_qles(self,ind):
		
		self.qle_center_x.setText(str(self.bkgd.centers_embr_px[ind][0]))
		self.qle_center_y.setText(str(self.bkgd.centers_embr_px[ind][1]))
		self.qle_radius.setText(str(self.bkgd.radiuses_embr_px[ind]))
	
	def update_prop_list_entry(self,ind):
		
		self.nodes[ind].setText(1,str(self.bkgd.centers_embr_px[ind]))
		self.nodes[ind].setText(2,str(self.bkgd.radiuses_embr_px[ind]))
	
	def init_prop_list(self):
		
		self.prop_list.clear()
		self.update_prop_list()
	
	def update_prop_list(self):
		
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.bkgd.fn_bkgdfolder,self.bkgd.data_ft)
		self.nodes=[]
		
		if len(self.file_list)==len(self.bkgd.radiuses_embr_px):
			pass
		else:
			self.bkgd.centers_embr_px=[[self.bkgd.res_px/2,self.bkgd.res_px/2]]*shape(self.file_list)[0]
			self.bkgd.radiuses_embr_px=[253.]*shape(self.file_list)[0]
		
		for i in range(shape(self.file_list)[0]):
			
			self.nodes.append(QtGui.QTreeWidgetItem(self.prop_list,[str(i),str(self.bkgd.centers_embr_px[i]),str(self.bkgd.radiuses_embr_px[i])]))
		
		if shape(self.file_list)[0]>0:
			self.prop_list.setCurrentItem(self.nodes[0])
			
	def create_canvas(self):
			
		h=500/self.dpi
		v=500/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.plot_frame)
		
		self.ax = self.fig.add_subplot(111)
		self.ax.set_xlim([0, self.bkgd.res_px])
		self.ax.set_ylim([0, self.bkgd.res_px])
		
		#Connect to mouse click event
		self.canvas.mpl_connect('button_press_event', self.get_mouse_embr_radius)
		
		self.canvas.draw()
		return 
	
	def draw_patches(self,ind):
		if shape(self.file_list)[0]>ind:
			self.clear_patch_canvas()
			
			if shape(self.bkgd.centers_embr_px[ind])[0]>0:
				pt=ptc.Circle([self.bkgd.centers_embr_px[ind][0],self.bkgd.centers_embr_px[ind][1]],radius=3,fill=True,color='r')
				self.pts_in_ax2.append(self.ax.add_patch(pt))
				self.xcoords2.append(self.bkgd.centers_embr_px[ind])
				self.fig.canvas.draw()	
			if shape(self.bkgd.radiuses_embr_px)[0]>0:
				embr_circ=ptc.Circle(self.bkgd.centers_embr_px[ind],radius=self.bkgd.radiuses_embr_px[ind],fill=False,color='r',linewidth=3)
				self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
				self.xcoords2.append(100)
				self.fig.canvas.draw()		
				
	def show_img(self,ind):
		
		self.file_list= pyfdap_misc.get_sorted_folder_list(self.bkgd.fn_bkgdfolder,self.bkgd.data_ft)
	
		#Check if there is a first image
		if shape(self.file_list)[0]>ind-1 and shape(self.file_list)[0]>0:
		
			#Grab first picture
			fn=self.bkgd.fn_bkgdfolder+self.file_list[ind]
			
			#Load img
			data_img = mpimg.imread(fn).astype(self.bkgd.data_enc)
			#data_img = skiio.imread(fn).astype(self.bkgd.data_enc)
			data_vals=data_img.real
			data_vals=data_vals.astype('float')
			
			#Flip to be coherent with Fiji
			data_vals=flipud(data_vals)
			
			#Plot img
			self.ax.imshow(data_vals)
			self.canvas.draw()
			
			#Grab resolution
			self.bkgd.res_px=shape(data_vals)[0]
			self.qle_res.setText(str(self.bkgd.res_px))
			
			self.curr_img_ind=ind
			
		else:
			#Clear if there is no first image
			self.ax.cla()
	
	def get_mouse_embr_radius(self,event):
		
		#Left click to define center and then radius
		if event.button==1:
			
			#Check if clicked within axes
			if event.xdata==None:
				return
			
			#Check if there is already a center
			if shape(self.xcoords2)[0]<1:
					
				#No center yet
				self.xcoords2.append(event.xdata)
				self.ycoords2.append(event.ydata)
				
				#Plot center
				pt=ptc.Circle([event.xdata,event.ydata],radius=3,fill=True,color='r')
				self.pts_in_ax2.append(self.ax.add_patch(pt))
				
				#Set center_embr
				self.bkgd.centers_embr_px[self.curr_img_ind]=[event.xdata,event.ydata]
				
				self.qle_center_x.setText(str(self.bkgd.centers_embr_px[self.curr_img_ind][0]))
				self.qle_center_y.setText(str(self.bkgd.centers_embr_px[self.curr_img_ind][1]))
				self.prop_list.currentItem().setText(1,str(self.bkgd.centers_embr_px[self.curr_img_ind]))
				
			elif shape(self.xcoords2)[0]==1:
				
				#Center yet, so this click is for radius
				self.xcoords2.append(event.xdata)
				self.ycoords2.append(event.ydata)
					
				#Set radius_embr
				self.bkgd.radiuses_embr_px[self.curr_img_ind]=sqrt((event.xdata-self.bkgd.centers_embr_px[self.curr_img_ind][0])**2+(event.ydata-self.bkgd.centers_embr_px[self.curr_img_ind][1])**2)
				
				#Plot circle
				embr_circ=ptc.Circle(self.bkgd.centers_embr_px[self.curr_img_ind],radius=self.bkgd.radiuses_embr_px[self.curr_img_ind],fill=False,color='r',linewidth=3)
				self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
				
				self.qle_radius.setText(str(self.bkgd.radiuses_embr_px[self.curr_img_ind]))
				self.prop_list.currentItem().setText(2,str(self.bkgd.radiuses_embr_px[self.curr_img_ind]))
			
			else:
			
				self.clear_patch_canvas()			
		
		self.fig.canvas.draw()			
	
	def clear_patch_canvas(self):
	
		#Clearing Canvas
		
		#Removing all points
		for pt in self.pts_in_ax2:
			pt.remove()
		
		#Removing circle
		if hasattr(self,'embr_circ_in_ax'):
			self.embr_circ_in_ax.remove()
			del self.embr_circ_in_ax

		#Setting back arrays
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		self.fig.canvas.draw()
		
	def done_pressed(self):
	
		self.done(1)
		return self.bkgd
	
#===================================================================================================================================
#Dialog for select/edit pre data set
#===================================================================================================================================

class pre_dataset_dialog(QtGui.QDialog):
	def __init__(self,pre,parent):
		super(pre_dataset_dialog,self).__init__(parent)
		
		self.pre=pre
			
		self.dpi = 100
		self.setMinimumSize(1000,500) 
		self.resize(1300,500)
		
		#Some variables
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_set_datafolder=QtGui.QPushButton('Change')
		self.btn_set_maskfolder=QtGui.QPushButton('Change')
		
		#Button Actions
		self.btn_set_datafolder.connect(self.btn_set_datafolder, QtCore.SIGNAL('clicked()'), self.sel_datafolder)
		self.btn_set_maskfolder.connect(self.btn_set_maskfolder, QtCore.SIGNAL('clicked()'), self.sel_maskfolder)
		
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_name_datafolder = QtGui.QLabel("Photoconverted Folder:", self)
		self.lbl_name_maskfolder = QtGui.QLabel("Counterlabeled Folder:", self)
		
		self.lbl_name_datafolder.setToolTip("red images")
		self.lbl_name_maskfolder.setToolTip("green images")
		self.lbl_name_ft = QtGui.QLabel("Filetype:", self)
		self.lbl_name_enc = QtGui.QLabel("Encoding:", self)
		self.lbl_name_res = QtGui.QLabel("Image dimensions (pixels):", self)
		self.lbl_name_radius = QtGui.QLabel("Embryo radius (pixels):", self)
		self.lbl_name_center = QtGui.QLabel("Embryo center (pixels):", self)
		
		self.lbl_datafolder = QtGui.QLabel(self.pre.fn_datafolder, self)
		self.lbl_maskfolder = QtGui.QLabel(self.pre.fn_maskfolder, self)
		
		if len(self.pre.fn_datafolder)>50:
			self.lbl_datafolder.setText("..."+self.pre.fn_datafolder[-50:])
		else:
			self.lbl_datafolder.setText(self.pre.fn_datafolder)
		
		if len(self.pre.fn_maskfolder)>50:
			self.lbl_maskfolder.setText("..."+self.pre.fn_maskfolder[-50:])
		else:
			self.lbl_maskfolder.setText(self.pre.fn_maskfolder)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_ft = QtGui.QComboBox(self)
		self.combo_ft.addItem("tif")
		
		self.combo_enc = QtGui.QComboBox(self)
		self.combo_enc.addItem("8bit")
		self.combo_enc.addItem("16bit")
		self.combo_enc.setCurrentIndex(1) 
		
		self.combo_ft.activated[str].connect(self.sel_ft)   
		self.combo_enc.activated[str].connect(self.sel_enc)   
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_res = QtGui.QLineEdit(str(self.pre.res_px))
		
		centergrid= QtGui.QGridLayout()
		if self.pre.center_embr_px[0]>0:
			self.qle_radius = QtGui.QLineEdit(str(self.pre.radius_embr_px))
			self.qle_center_x = QtGui.QLineEdit(str(self.pre.center_embr_px[0]))
			self.qle_center_y = QtGui.QLineEdit(str(self.pre.center_embr_px[1]))
		else:
			self.qle_radius = QtGui.QLineEdit("")
			self.qle_center_x = QtGui.QLineEdit("")
			self.qle_center_y = QtGui.QLineEdit("")
			
		centergrid.addWidget(self.qle_center_x,0,0)
		centergrid.addWidget(self.qle_center_y,0,1)
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_res.setValidator(self.double_valid)
		self.qle_radius.setValidator(self.double_valid)
		self.qle_center_x.setValidator(self.double_valid)
		self.qle_center_y.setValidator(self.double_valid)
		
		self.qle_res.editingFinished.connect(self.set_res)
		self.qle_radius.editingFinished.connect(self.set_radius)
		self.qle_center_x.editingFinished.connect(self.set_center_x)
		self.qle_center_y.editingFinished.connect(self.set_center_y)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Plot frame
		#-------------------------------------------------------------------------------------------------------------------
		
		self.plot_frame = QtGui.QWidget()
		self.plot_frame.setMaximumWidth(1)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_name_datafolder,2,1)
		grid.addWidget(self.lbl_name_maskfolder,3,1)
		grid.addWidget(self.lbl_name_ft,4,1)
		grid.addWidget(self.lbl_name_enc,5,1)
		grid.addWidget(self.lbl_name_res,6,1)
		grid.addWidget(self.lbl_name_radius,7,1)
		grid.addWidget(self.lbl_name_center,8,1)
		grid.addWidget(self.btn_done,10,1)
		
		grid.addWidget(self.lbl_datafolder,2,2)
		grid.addWidget(self.lbl_maskfolder,3,2)
		grid.addWidget(self.combo_ft,4,2)
		grid.addWidget(self.combo_enc,5,2)
		grid.addWidget(self.qle_res,6,2)
		grid.addWidget(self.qle_radius,7,2)
		grid.addLayout(centergrid,8,2)
		
		grid.addWidget(self.btn_set_datafolder,2,3)
		grid.addWidget(self.btn_set_maskfolder,3,3)
		
		grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Create Canvas
		#-------------------------------------------------------------------------------------------------------------------
		
		self.create_canvas()
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.canvas)
		
		#Add everything to Horizontal Box
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(grid)
		self.hbox.addWidget(self.canvas)
		
		self.setLayout(self.hbox)    
		
		#Check if there is already a first image
		if self.pre.fn_datafolder!="":
			self.file_list = pyfdap_misc.get_sorted_folder_list(self.pre.fn_datafolder,self.pre.data_ft)
			self.show_img()	
			self.draw_patches()
		else:
			self.sel_auto_pre_folders()
		
		self.setWindowTitle('Edit Preconversion Dataset')    
		self.show()
	
	def create_frame(self,frame):
		
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		frame.setBackgroundRole(QtGui.QPalette.Light)      
		frame.setLineWidth (1)
		frame.setFrameShadow (frame.Sunken)
		frame.frameRect().setWidth(100)
		
		return frame
		
	def sel_datafolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Data Directory"))
		if folder=='':
			return
		
		self.pre.fn_datafolder=folder
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.pre.fn_datafolder,self.pre.data_ft)
		self.pre.center_embr_px=[self.pre.res_px/2,self.pre.res_px/2]
		self.pre.radiuses_embr_px=253.
		
		#If folder string ends without /, add /
		if self.pre.fn_datafolder[-1]=="/":
			pass
		else:
			self.pre.fn_datafolder=self.pre.fn_datafolder+"/"
		
		if len(self.pre.fn_datafolder)>50:
			self.lbl_datafolder.setText("..."+self.pre.fn_datafolder[-50:])
		else:
			self.lbl_datafolder.setText(self.pre.fn_datafolder)
		
		self.show_img()
		self.update_qles()
		
		self.sel_other_folder(folder,"green")
		
	def sel_maskfolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Mask Data Directory"))
		if folder=='':
			return
		self.pre.fn_maskfolder = folder
		
		#If folder string ends without /, add /
		if self.pre.fn_maskfolder[-1]=="/":
			pass
		else:
			self.pre.fn_maskfolder=self.pre.fn_maskfolder+"/"
		
		if len(self.pre.fn_maskfolder)>50:
			self.lbl_maskfolder.setText("..."+self.pre.fn_maskfolder[-50:])
		else:
			self.lbl_maskfolder.setText(self.pre.fn_maskfolder)
		
		self.sel_other_folder(folder,"red")
		
	def sel_other_folder(self,folder,color):
		
		#Delete last "/" for os.path.split to work
		if folder[-1]=="/":
			folder=folder[:-1]
		
		#Getting upper folder
		[upper_folder,last_folder]=os.path.split(folder)
		
		#Check if other folder is where it should be, select if exists and is not empty
		otherfolder=upper_folder+"/"+color
		if os.path.isdir(otherfolder):
			
			#Check of folder is empty
			if pyfdap_misc.check_folder_empty(otherfolder):
				pass
				#print "Did not automatically select ", otherfolder, "as folder. Folder is empty!"
			else: 	
				#Add / so analysis functions don't get confused
				otherfolder=otherfolder+"/"
				
				if color=="green":
					#Update maskfolder
					self.pre.fn_maskfolder = otherfolder
					if len(self.pre.fn_maskfolder)>50:
						self.lbl_maskfolder.setText("..."+self.pre.fn_maskfolder[-50:])
					else:
						self.lbl_maskfolder.setText(self.pre.fn_maskfolder)
				
				elif color=="red":
					#Check if otherfolder is different from datafolder so we don't override any settings
					if self.pre.fn_datafolder!=otherfolder:
						
						self.pre.fn_datafolder=otherfolder
						self.file_list = pyfdap_misc.get_sorted_folder_list(self.pre.fn_datafolder,self.pre.data_ft)
						self.pre.center_embr_px=[self.pre.res_px/2,self.pre.res_px/2]
						self.pre.radiuses_embr_px=253.
						
						if len(self.pre.fn_datafolder)>50:
							self.lbl_datafolder.setText("..."+self.pre.fn_datafolder[-50:])
						else:
							self.lbl_datafolder.setText(self.pre.fn_datafolder)
						
						self.show_img()
						self.update_qles()
		else:
			pass
			#print "Did not automatically select ", otherfolder, "as folder. Folder does not exists!"
	
	def sel_auto_pre_folders(self):
		
		if hasattr(self.pre,'embryo'):
			datafolder="/"+self.pre.embryo.fn_datafolder.strip("red/")+"/"
		elif hasattr(self.pre,'bkgd'):
			datafolder="/"+self.pre.bkgd.fn_bkgdfolder.strip("red/")+"/"
		
		if os.path.isdir(datafolder):
			prefolder=datafolder+'pre/red/'
			
			if os.path.isdir(prefolder):
				
				if not pyfdap_misc.check_folder_empty(prefolder):
					
					self.pre.fn_datafolder=prefolder
					self.file_list = pyfdap_misc.get_sorted_folder_list(self.pre.fn_datafolder,self.pre.data_ft)
					self.pre.center_embr_px=[self.pre.res_px/2,self.pre.res_px/2]
					self.pre.radiuses_embr_px=253.
					
					if len(self.pre.fn_datafolder)>50:
						self.lbl_datafolder.setText("..."+self.pre.fn_datafolder[-50:])
					else:
						self.lbl_datafolder.setText(self.pre.fn_datafolder)
					
					self.sel_other_folder(prefolder,"green")
					
					self.show_img()
					self.update_qles()
				
	def sel_ft(self,text):
		if text=="tif":
			self.pre.data_ft=text
		
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.pre.fn_datafolder,self.pre.data_ft)		
		self.show_img()
		
	def sel_enc(self,text):
		if text=="8bit":
			self.pre.data_enc="uint8"
		elif text=="16bit":	
			self.pre.data_enc="uint16"
		self.show_img()
	
	def set_name(self):
		text=self.qle_name.text()
		self.pre.name=str(text)
	
	def set_res(self):
		text=self.qle_res.text()
		self.pre.res_px=float(str(text))
		self.ax.set_xlim([0, self.pre.res_px])
		self.ax.set_ylim([0, self.pre.res_px])
		
	def set_radius(self):
		
		text=self.qle_radius.text()
		self.pre.radius_embr_px=float(str(text))	
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.pre.center_embr_px[0],self.pre.center_embr_px[1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.pre.center_embr_px[0])
		
		embr_circ=ptc.Circle(self.pre.center_embr_px,radius=self.pre.radius_embr_px,fill=False,color='r',linewidth=3)
		self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
		self.xcoords2.append(self.pre.center_embr_px[0])
	
		self.fig.canvas.draw()	
		
	def set_center_x(self):
		text=self.qle_center_x.text()
		self.pre.center_embr_px[0]=float(str(text))
			
		if hasattr(self,'embr_circ_in_ax'):
			had_circ=1
		else:
			had_circ=0
			
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.pre.center_embr_px[0],self.pre.center_embr_px[1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.pre.center_embr_px[0])
		
		if had_circ==1:
			embr_circ=ptc.Circle(self.pre.center_embr_px,radius=self.pre.radius_embr_px,fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(self.pre.center_embr_px[0])
		
		self.fig.canvas.draw()
		
	def set_center_y(self):
		text=self.qle_center_y.text()
		self.pre.center_embr_px[1]=float(str(text))
		
		if hasattr(self,'embr_circ_in_ax'):
			had_circ=1
		else:
			had_circ=0
			
		self.clear_patch_canvas()
		
		pt=ptc.Circle([self.pre.center_embr_px[0],self.pre.center_embr_px[1]],radius=3,fill=True,color='r')
		self.pts_in_ax2.append(self.ax.add_patch(pt))
		self.xcoords2.append(self.pre.center_embr_px[0])
		
		if had_circ==1:
			embr_circ=ptc.Circle(self.pre.center_embr_px,radius=self.pre.radius_embr_px,fill=False,color='r',linewidth=3)
			self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
			self.xcoords2.append(self.pre.center_embr_px[0])
		
		self.fig.canvas.draw()
		
	def update_qles(self):
		
		self.qle_center_x.setText(str(self.pre.center_embr_px[0]))
		self.qle_center_y.setText(str(self.pre.center_embr_px[1]))
		self.qle_radius.setText(str(self.pre.radius_embr_px))
							
	def create_canvas(self):
			
		h=500/self.dpi
		v=500/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.plot_frame)
		
		self.ax = self.fig.add_subplot(111)
		self.ax.set_xlim([0, self.pre.res_px])
		self.ax.set_ylim([0, self.pre.res_px])
		
		#Connect to mouse click event
		self.canvas.mpl_connect('button_press_event', self.get_mouse_embr_radius)
		
		self.canvas.draw()
		return 
	
	def draw_patches(self):
		if shape(self.file_list)[0]>0:
			self.clear_patch_canvas()
			
			if self.pre.center_embr_px>0:
				pt=ptc.Circle([self.pre.center_embr_px[0],self.pre.center_embr_px[1]],radius=3,fill=True,color='r')
				self.pts_in_ax2.append(self.ax.add_patch(pt))
				self.xcoords2.append(self.pre.center_embr_px)
				self.fig.canvas.draw()	
			if self.pre.radius_embr_px>0:
				embr_circ=ptc.Circle(self.pre.center_embr_px,radius=self.pre.radius_embr_px,fill=False,color='r',linewidth=3)
				self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
				self.xcoords2.append(100)
				self.fig.canvas.draw()		
				
	def show_img(self):
		
		#Check if there is a first image
		if shape(self.file_list)[0]>0:
			
			#Grab first picture
			fn=self.pre.fn_datafolder+self.file_list[0]
			
			#Load img
			#data_img = skiio.imread(fn).astype(self.pre.data_enc)
			data_img = mpimg.imread(fn).astype(self.pre.data_enc)
			data_vals=data_img.real
			data_vals=data_vals.astype('float')
			
			#Flip to be coherent with Fiji
			data_vals=flipud(data_vals)
			
			#Plot img
			self.ax.imshow(data_vals)
			self.canvas.draw()
			
			#Grab resolution
			self.pre.res_px=shape(data_vals)[0]
			self.qle_res.setText(str(self.pre.res_px))
			
		else:
			#Clear if there is no first image
			self.ax.cla()
	
	
	def get_mouse_embr_radius(self,event):
		
		#Left click to define center and then radius
		if event.button==1:
			
			#Check if clicked withing axes
			if event.xdata==None:
				return
			
			#Check if there is already a center
			if shape(self.xcoords2)[0]<1:
					
				#No center yet
				self.xcoords2.append(event.xdata)
				self.ycoords2.append(event.ydata)
				
				#Plot center
				pt=ptc.Circle([event.xdata,event.ydata],radius=3,fill=True,color='r')
				self.pts_in_ax2.append(self.ax.add_patch(pt))
				
				#Set center_embr
				self.pre.center_embr_px=[event.xdata,event.ydata]
				
				self.qle_center_x.setText(str(self.pre.center_embr_px[0]))
				self.qle_center_y.setText(str(self.pre.center_embr_px[1]))
				
			elif shape(self.xcoords2)[0]==1:
				
				#Center yet, so this click is for radius
				self.xcoords2.append(event.xdata)
				self.ycoords2.append(event.ydata)
					
				#Set radius_embr
				self.pre.radius_embr_px=sqrt((event.xdata-self.pre.center_embr_px[0])**2+(event.ydata-self.pre.center_embr_px[1])**2)
				
				#Plot circle
				embr_circ=ptc.Circle(self.pre.center_embr_px,radius=self.pre.radius_embr_px,fill=False,color='r',linewidth=3)
				self.embr_circ_in_ax=self.ax.add_patch(embr_circ)
				
				self.qle_radius.setText(str(self.pre.radius_embr_px))
			
			else:
			
				self.clear_patch_canvas()
				
		self.fig.canvas.draw()			
	
	def clear_patch_canvas(self):
	
		#Clearing Canvas
		
		#Removing all points
		for pt in self.pts_in_ax2:
			pt.remove()
		
		#Removing circle
		if hasattr(self,'embr_circ_in_ax'):
			self.embr_circ_in_ax.remove()
			del self.embr_circ_in_ax

		#Setting back arrays
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		self.fig.canvas.draw()
		
	def done_pressed(self):
	
		self.done(1)
		return self.pre

#===================================================================================================================================
#Dialog for select/edit noise data set
#===================================================================================================================================

class noise_dataset_dialog(QtGui.QDialog):
	def __init__(self,noise,parent):
		super(noise_dataset_dialog,self).__init__(parent)
		
		self.noise=noise
			
		self.dpi = 100
		self.setMinimumSize(400,300) 
		self.resize(400,300)
		
		#Some variables
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_set_datafolder=QtGui.QPushButton('Change')
		
		#Button Actions
		self.btn_set_datafolder.connect(self.btn_set_datafolder, QtCore.SIGNAL('clicked()'), self.sel_datafolder)
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_name_mode = QtGui.QLabel("Mode:", self)
		self.lbl_name_datafolder = QtGui.QLabel("Datafolder:", self)
		self.lbl_name_noise = QtGui.QLabel("Noise level:", self)
		self.lbl_name_ft = QtGui.QLabel("Filetype:", self)
		self.lbl_name_enc = QtGui.QLabel("Encoding:", self)
		self.lbl_name_res = QtGui.QLabel("Resolution:", self)
		
		self.lbl_datafolder = QtGui.QLabel(self.noise.fn_datafolder, self)
		
		if len(self.noise.fn_datafolder)>50:
			self.lbl_datafolder.setText("..."+self.noise.fn_datafolder[-50:])
		else:
			self.lbl_datafolder.setText(self.noise.fn_datafolder)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_mode = QtGui.QComboBox(self)
		self.combo_mode.addItem("Outside")
		self.combo_mode.addItem("Predefined")
		self.combo_mode.addItem("Seperate Dataset")
		self.combo_mode.setCurrentIndex(0) 
		
		self.combo_ft = QtGui.QComboBox(self)
		self.combo_ft.addItem("tif")
		
		self.combo_enc = QtGui.QComboBox(self)
		self.combo_enc.addItem("8bit")
		self.combo_enc.addItem("16bit")
		self.combo_enc.setCurrentIndex(1) 
		
		self.combo_mode.activated[str].connect(self.sel_mode)   
		self.combo_ft.activated[str].connect(self.sel_ft)   
		self.combo_enc.activated[str].connect(self.sel_enc)   
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_res = QtGui.QLineEdit(str(self.noise.res_px))
		self.qle_noise	= QtGui.QLineEdit(str(self.noise.noise))
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_res.setValidator(self.double_valid)
		self.qle_noise.setValidator(self.double_valid)
		
		self.qle_res.editingFinished.connect(self.set_res)
		self.qle_noise.editingFinished.connect(self.set_noise)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_name_mode,1,1)
		grid.addWidget(self.lbl_name_datafolder,2,1)
		grid.addWidget(self.lbl_name_noise,3,1)
		grid.addWidget(self.lbl_name_ft,4,1)
		grid.addWidget(self.lbl_name_enc,5,1)
		grid.addWidget(self.lbl_name_res,6,1)
		grid.addWidget(self.btn_done,10,1)
		
		grid.addWidget(self.combo_mode,1,2)
		grid.addWidget(self.lbl_datafolder,2,2)
		grid.addWidget(self.qle_noise,3,2)
		grid.addWidget(self.combo_ft,4,2)
		grid.addWidget(self.combo_enc,5,2)
		grid.addWidget(self.qle_res,6,2)
	
		grid.addWidget(self.btn_set_datafolder,2,3)
		
		grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.setLayout(grid)    
		self.update_mode("Outside")	
		self.setWindowTitle('Edit Noise Dataset')    
		self.show()
		
	def sel_datafolder(self):
		folder = str(QFileDialog.getExistingDirectory(self, "Select Noise Directory"))
		if folder=='':
			return
		
		self.noise.fn_datafolder=folder
		
		#If folder string ends without /, add /
		if self.noise.fn_datafolder[-1]=="/":
			pass
		else:
			self.noise.fn_datafolder=self.noise.fn_datafolder+"/"
		
		if len(self.noise.fn_datafolder)>50:
			self.lbl_datafolder.setText("..."+self.noise.fn_datafolder[-50:])
		else:
			self.lbl_datafolder.setText(self.noise.fn_datafolder)
		
		self.update_mode("Seperate Dataset")
		
	def sel_ft(self,text):
		if text=="tif":
			self.noise.data_ft=text
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.noise.fn_datafolder,self.noise.data_ft)		
			
	def sel_enc(self,text):
		if text=="8bit":
			self.noise.data_enc="uint8"
		elif text=="16bit":	
			self.noise.data_enc="uint16"
	
	def sel_mode(self,text):
		self.update_mode(text)
	
	def update_mode(self,text):
		if text=="Outside":
			self.noise.mode="outside"
			self.qle_res.setReadOnly(True)
			self.qle_noise.setReadOnly(True)
			self.qle_res.setText(str(self.noise.embryo.data_res_px))
			self.qle_noise.setText("")
			self.lbl_datafolder.setHidden(True)
			
			if self.noise.embryo.data_enc=="uint8":
				self.combo_enc.setCurrentIndex(0)
			elif self.noise.embryo.data_enc=="uint16":	
				self.combo_enc.setCurrentIndex(1)
			
			self.combo_enc.setEnabled(False)
			
			if self.noise.embryo.data_ft=="tif":
				self.combo_ft.setCurrentIndex(0)
			
			self.combo_ft.setEnabled(False)
			self.btn_set_datafolder.setEnabled(False)
			
		elif text=="Predefined":	
			self.noise.mode="predefined"
			self.qle_res.setReadOnly(True)
			self.qle_noise.setReadOnly(False)
			self.qle_noise.setText(str(self.noise.noise))
			self.qle_res.setText("")
			self.lbl_datafolder.setText("")
			self.lbl_datafolder.setHidden(True)
			self.btn_set_datafolder.setEnabled(False)
			self.combo_enc.setEnabled(False)
			self.combo_ft.setEnabled(False)
			
		elif text=="Seperate Dataset":	
			self.noise.mode="seperate"	
			self.qle_res.setReadOnly(False)
			self.qle_noise.setReadOnly(True)
			self.qle_noise.setText("")
			self.qle_res.setText("")
			self.btn_set_datafolder.setEnabled(True)
			
			if len(self.noise.fn_datafolder)>50:
				self.lbl_datafolder.setText("..."+self.noise.fn_datafolder[-50:])
			else:
				self.lbl_datafolder.setText(self.noise.fn_datafolder)
			
			self.file_list = pyfdap_misc.get_sorted_folder_list(self.noise.fn_datafolder,self.noise.data_ft)
			if shape(self.file_list)[0]>0:
				fn=self.noise.fn_datafolder+self.file_list[0]
				data_img = mpimg.imread(fn).astype(self.noise.data_enc)
				#data_img = skiio.imread(fn).astype(self.noise.data_enc)
				data_vals=data_img.real
				data_vals=data_vals.astype('float')
				
				self.noise.res_px=shape(data_vals)[0]
				
				self.qle_res.setText(str(self.noise.res_px))
			
			self.lbl_datafolder.setHidden(False)
			self.combo_enc.setEnabled(True)
			self.combo_ft.setEnabled(True)
	
	def set_res(self):
		text=self.qle_res.text()
		if text!="":
			self.noise.res_px=float(str(text))
		
	def set_noise(self):	
		text=self.qle_noise.text()
		if text!="":
			self.noise.noise=float(str(text))	
										
	def done_pressed(self):
	
		self.done(1)
		return self.noise

#===================================================================================================================================
#Dialog for select/edit fit
#===================================================================================================================================

class fit_dialog(QtGui.QDialog):
	
	def __init__(self,fit,molecule,embryo,parent):
		super(fit_dialog,self).__init__(parent)
		
		self.fit=fit
		self.molecule=molecule
		self.embryo=embryo
		
		self.temp_LB_k=[]
		self.temp_UB_k=[]
		self.temp_LB_cnaught=[]
		self.temp_UB_cnaught=[]
		self.temp_LB_ynaught=[]
		self.temp_UB_ynaught=[]
		
		self.setMinimumSize(600,300) 
		
		#Make sure fit object has model
		if hasattr(self.fit,"model"):
			pass
		else:
			self.fit.model="exp"
			self.fit.npower=1
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lbl_opt_parms = QtGui.QLabel("Optimization Parameters", self)
		self.lbl_opt_parms.setFont(boldfont)
		self.lbl_guess = QtGui.QLabel("Initial Guess", self)
		self.lbl_guess.setFont(boldfont)
		self.lbl_bounds = QtGui.QLabel("Variable bounds", self)
		self.lbl_bounds.setFont(boldfont)
		self.lbl_bounded = QtGui.QLabel("Bounded?", self)
		self.lbl_bounded.setFont(boldfont)
		self.lbl_fitting = QtGui.QLabel("Fitting options", self)
		self.lbl_fitting.setFont(boldfont)
		self.lbl_model_head = QtGui.QLabel("Model", self)
		self.lbl_model_head.setFont(boldfont)
		
		self.lbl_name = QtGui.QLabel("Name:", self)
		self.lbl_opt_meth = QtGui.QLabel("opt_meth:", self)
		self.lbl_opt_tol = QtGui.QLabel("opt_tol", self)
		self.lbl_maxfun = QtGui.QLabel("maxfun:", self)
		self.lbl_save_track = QtGui.QLabel("save_track:", self)
		self.lbl_debug_fit = QtGui.QLabel("debug_fit:", self)
		
		self.lbl_x0_k = QtGui.QLabel("x0_k:", self)
		self.lbl_x0_c0 = QtGui.QLabel("x0_c0:", self)
		self.lbl_x0_y0 = QtGui.QLabel("x0_y0:", self)
		
		self.lbl_npower = QtGui.QLabel("n:", self) 
		self.lbl_model =  QtGui.QLabel("Model:", self)
		
		self.lbl_LB_k = QtGui.QLabel("LB_k:", self)
		self.lbl_UB_k = QtGui.QLabel("UB_k:", self)
		self.lbl_LB_y0 = QtGui.QLabel("LB_y0:", self)
		self.lbl_UB_y0 = QtGui.QLabel("UB_y0:", self)
		self.lbl_LB_c0 = QtGui.QLabel("LB_c0:", self)
		self.lbl_UB_c0 = QtGui.QLabel("UB_c0:", self)
		
		self.lbl_fit_ext = QtGui.QLabel("fit_ext:", self)
		self.lbl_fit_int = QtGui.QLabel("fit_int:", self)
		self.lbl_fit_slice = QtGui.QLabel("fit_slice:", self)
		self.lbl_fit_y0 = QtGui.QLabel("fit_y0:", self)
		self.lbl_fit_c0 = QtGui.QLabel("fit_c0:", self)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_meth = QtGui.QComboBox(self)
		self.combo_meth.addItem("Constrained Nelder-Mead")
		self.combo_meth.addItem("TNC")
		self.combo_meth.addItem("Nelder-Mead")
		self.combo_meth.addItem("L-BFGS-B")
		self.combo_meth.addItem("SLSQP")
		self.combo_meth.addItem("brute")
		self.combo_meth.addItem("BFGS")
		self.combo_meth.addItem("CG")
		self.combo_meth.activated[str].connect(self.sel_meth)   
		
		self.combo_x0_c0 = QtGui.QComboBox(self)
		self.combo_x0_c0.addItem("Custom")
		self.combo_x0_c0.addItem("I_post(tstart)-I_pre")
		self.combo_x0_c0.activated[str].connect(self.sel_x0_c0)   
		
		self.combo_x0_y0 = QtGui.QComboBox(self)
		self.combo_x0_y0.addItem("Custom")
		self.combo_x0_y0.addItem("I_post(tstart)")
		self.combo_x0_y0.addItem("I_post(tend)")
		self.combo_x0_y0.addItem("I_pre")
		self.combo_x0_y0.activated[str].connect(self.sel_x0_y0)   
		
		self.combo_UB_y0 = QtGui.QComboBox(self)
		self.combo_UB_y0.addItem("Custom")
		self.combo_UB_y0.addItem("Max(I_int(t))")
		self.combo_UB_y0.activated[str].connect(self.sel_UB_y0)   
		
		self.combo_LB_y0 = QtGui.QComboBox(self)
		self.combo_LB_y0.addItem("Custom")
		self.combo_LB_y0.addItem("Noise")
		self.combo_LB_y0.addItem("Bkgd_pre")
		self.combo_LB_y0.addItem("Bkgd")
		self.combo_LB_y0.addItem("F")
		self.curr_LBmode="Custom"
		self.combo_LB_y0.activated[str].connect(self.sel_LB_y0) 
		
		self.combo_model = QtGui.QComboBox(self)
		self.combo_model.addItem("exp")
		self.combo_model.addItem("power")
		self.combo_model.activated[str].connect(self.sel_model) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_debug_fit = QtGui.QCheckBox('', self)
		self.cb_save_track = QtGui.QCheckBox('', self)
		
		self.cb_bound_LB_k = QtGui.QCheckBox('', self)
		self.cb_bound_LB_c0 = QtGui.QCheckBox('', self)
		self.cb_bound_LB_y0 = QtGui.QCheckBox('', self)
		
		self.cb_bound_UB_k = QtGui.QCheckBox('', self)
		self.cb_bound_UB_c0 = QtGui.QCheckBox('', self)
		self.cb_bound_UB_y0 = QtGui.QCheckBox('', self)
		
		self.cb_fit_ext = QtGui.QCheckBox('', self)
		self.cb_fit_int = QtGui.QCheckBox('', self)
		self.cb_fit_slice = QtGui.QCheckBox('', self)
		
		self.cb_fit_c0 = QtGui.QCheckBox('', self)
		self.cb_fit_y0 = QtGui.QCheckBox('', self)
		
		self.cb_debug_fit.setCheckState(2*self.fit.debug_fit)	
		self.cb_save_track.setCheckState(2*self.fit.save_track)	
	
		if self.fit.LB_k==None:
			self.cb_bound_LB_k.setCheckState(0)
		else: 
			self.cb_bound_LB_k.setCheckState(2)
		
		if self.fit.UB_k==None:
			self.cb_bound_UB_k.setCheckState(0)
		else: 
			self.cb_bound_UB_k.setCheckState(2)
		
		if self.fit.LB_cnaught==None:
			self.cb_bound_LB_c0.setCheckState(0)
		else: 
			self.cb_bound_LB_c0.setCheckState(2)
		
		if self.fit.UB_cnaught==None:
			self.cb_bound_UB_c0.setCheckState(0)
		else: 
			self.cb_bound_UB_c0.setCheckState(2)
		
		if self.fit.LB_ynaught==None:
			self.cb_bound_LB_y0.setCheckState(0)
		else: 
			self.cb_bound_LB_y0.setCheckState(2)
		
		if self.fit.UB_ynaught==None:
			self.cb_bound_UB_y0.setCheckState(0)
		else: 
			self.cb_bound_UB_y0.setCheckState(2)
		
		self.connect(self.cb_bound_LB_k, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_k)
		self.connect(self.cb_bound_LB_c0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_c0)
		self.connect(self.cb_bound_LB_y0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_y0)
		self.connect(self.cb_bound_UB_k, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_k)
		self.connect(self.cb_bound_UB_c0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_c0)
		self.connect(self.cb_bound_UB_y0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_y0)
		
		self.cb_fit_ext.setCheckState(2*self.fit.fit_ext)
		self.cb_fit_int.setCheckState(2*self.fit.fit_int)
		self.cb_fit_slice.setCheckState(2*self.fit.fit_slice)
		self.cb_fit_c0.setCheckState(2*self.fit.fit_cnaught)
		self.cb_fit_y0.setCheckState(2*self.fit.fit_ynaught)
		
		self.connect(self.cb_debug_fit, QtCore.SIGNAL('stateChanged(int)'), self.check_debug_fit)
		self.connect(self.cb_save_track, QtCore.SIGNAL('stateChanged(int)'), self.check_save_track)
		
		self.connect(self.cb_fit_ext, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_ext)
		self.connect(self.cb_fit_int, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_int)
		self.connect(self.cb_fit_slice, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_slice)
		self.connect(self.cb_fit_c0, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_c0)
		self.connect(self.cb_fit_y0, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_y0)
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_name = QtGui.QLineEdit(self.fit.name)
		self.qle_opt_tol = QtGui.QLineEdit(str(self.fit.opt_tol))
		self.qle_maxfun = QtGui.QLineEdit(str(self.fit.maxfun))
		
		self.qle_x0_k = QtGui.QLineEdit(str(self.fit.x0[0]))
		self.qle_x0_c0 = QtGui.QLineEdit(str(self.fit.x0[1]))
		self.qle_x0_y0 = QtGui.QLineEdit(str(self.fit.x0[2]))
		self.qle_npower = QtGui.QLineEdit(str(self.fit.npower))
	
		if self.fit.LB_k==None:
			self.qle_LB_k = QtGui.QLineEdit("0.0")
			self.qle_LB_k.setReadOnly(True)
		else:	
			
			self.qle_LB_k = QtGui.QLineEdit(str(self.fit.LB_k))
			
		if self.fit.UB_k==None:
			self.qle_UB_k = QtGui.QLineEdit("0.0")
			self.qle_UB_k.setReadOnly(True)
		else:
			self.qle_UB_k = QtGui.QLineEdit(str(self.fit.UB_k))
			
		if self.fit.LB_cnaught==None:
			self.qle_LB_c0 = QtGui.QLineEdit("0.0")
			self.qle_LB_c0.setReadOnly(True)
		else:
			self.qle_LB_c0 = QtGui.QLineEdit(str(self.fit.LB_cnaught))
			
		if self.fit.UB_cnaught==None:
			self.qle_UB_c0 = QtGui.QLineEdit("0.0")
			self.qle_UB_c0.setReadOnly(True)
		else:
			self.qle_UB_c0 = QtGui.QLineEdit(str(self.fit.UB_cnaught))
			
		if self.fit.LB_ynaught==None:
			self.qle_LB_y0 = QtGui.QLineEdit("0.0")
			self.qle_LB_y0.setReadOnly(True)
		else:
			self.qle_LB_y0 = QtGui.QLineEdit(str(self.fit.LB_ynaught))
			
		if self.fit.UB_ynaught==None:			
			self.qle_UB_y0 = QtGui.QLineEdit("0.0")
			self.qle_UB_y0.setReadOnly(True)
		else:
			self.qle_UB_y0 = QtGui.QLineEdit(str(self.fit.UB_ynaught))
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_LB_k.setValidator(self.double_valid)
		self.qle_LB_c0.setValidator(self.double_valid)
		self.qle_LB_y0.setValidator(self.double_valid)
		self.qle_UB_k.setValidator(self.double_valid)
		self.qle_UB_c0.setValidator(self.double_valid)
		self.qle_UB_y0.setValidator(self.double_valid)
		
		self.int_valid=QtGui.QIntValidator()
		self.qle_npower.setValidator(self.int_valid)
		
		self.qle_name.textChanged[str].connect(self.set_name)
		self.qle_opt_tol.textChanged[str].connect(self.set_opt_tol)
		self.qle_maxfun.textChanged[str].connect(self.set_maxfun)
		
		self.qle_x0_k.textChanged[str].connect(self.set_x0_k)
		self.qle_x0_c0.textChanged[str].connect(self.set_x0_c0)
		self.qle_x0_y0.textChanged[str].connect(self.set_x0_y0)
		self.qle_npower.textChanged[str].connect(self.set_npower)
		
		self.qle_UB_k.textChanged[str].connect(self.set_UB_k)
		self.qle_UB_c0.textChanged[str].connect(self.set_UB_c0)
		self.qle_UB_y0.textChanged[str].connect(self.set_UB_y0)
		
		self.qle_LB_k.textChanged[str].connect(self.set_LB_k)
		self.qle_LB_c0.textChanged[str].connect(self.set_LB_c0)
		self.qle_LB_y0.textChanged[str].connect(self.set_LB_y0)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_opt_parms,0,1,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_model_head,0,3,1,3,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_guess,0,3,9,3,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_bounds,0,6,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounded,0,9,1,1,Qt.AlignHCenter)
		grid.addWidget(self.lbl_fitting,0,10,1,2,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_name,1,1)
		grid.addWidget(self.lbl_opt_meth,2,1)
		grid.addWidget(self.lbl_opt_tol,3,1)
		grid.addWidget(self.lbl_maxfun,4,1)
		grid.addWidget(self.lbl_debug_fit,5,1)
		grid.addWidget(self.lbl_save_track,6,1)
		
		grid.addWidget(self.qle_name,1,2)
		grid.addWidget(self.combo_meth,2,2)
		grid.addWidget(self.qle_opt_tol,3,2)
		grid.addWidget(self.qle_maxfun,4,2)
		grid.addWidget(self.cb_debug_fit,5,2)
		grid.addWidget(self.cb_save_track,6,2)
		
		grid.addWidget(self.lbl_model,1,3)
		grid.addWidget(self.lbl_npower,2,3)
		
		grid.addWidget(self.lbl_x0_k,5,3)
		grid.addWidget(self.lbl_x0_c0,6,3)
		grid.addWidget(self.lbl_x0_y0,7,3)
		
		grid.addWidget(self.combo_model,1,4)
		grid.addWidget(self.qle_npower,2,4)
		
		grid.addWidget(self.qle_x0_k,5,4)
		grid.addWidget(self.qle_x0_c0,6,4)
		grid.addWidget(self.qle_x0_y0,7,4)
		
		grid.addWidget(self.combo_x0_c0,6,5)
		grid.addWidget(self.combo_x0_y0,7,5)
		
		grid.addWidget(self.lbl_LB_k,1,6)
		grid.addWidget(self.lbl_LB_c0,3,6)
		grid.addWidget(self.lbl_LB_y0,5,6)
		
		grid.addWidget(self.lbl_UB_k,2,6)
		grid.addWidget(self.lbl_UB_c0,4,6)
		grid.addWidget(self.lbl_UB_y0,6,6)
		
		grid.addWidget(self.qle_LB_k,1,7)
		grid.addWidget(self.qle_LB_c0,3,7)
		grid.addWidget(self.qle_LB_y0,5,7)
		
		grid.addWidget(self.qle_UB_k,2,7)
		grid.addWidget(self.qle_UB_c0,4,7)
		grid.addWidget(self.qle_UB_y0,6,7)
		
		grid.addWidget(self.combo_LB_y0,5,8)
		grid.addWidget(self.combo_UB_y0,6,8)
		
		grid.addWidget(self.cb_bound_LB_k,1,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_c0,3,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_y0,5,9,Qt.AlignHCenter)
		
		grid.addWidget(self.cb_bound_UB_k,2,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_c0,4,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_y0,6,9,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_fit_ext,1,10)
		grid.addWidget(self.lbl_fit_int,2,10)
		grid.addWidget(self.lbl_fit_slice,3,10)
		grid.addWidget(self.lbl_fit_c0,4,10)
		grid.addWidget(self.lbl_fit_y0,5,10)
		
		grid.addWidget(self.cb_fit_ext,1,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_int,2,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_slice,3,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_c0,4,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_y0,5,11,Qt.AlignHCenter)
		
		grid.addWidget(self.btn_done,8,11)
		
		grid.setColumnStretch(0,1)
		grid.setColumnStretch(12,1)
	
		grid.setRowStretch(9,1)
		
		self.setLayout(grid)    
		self.setWindowTitle('Edit Fit')    
		self.power_vis()
		self.show()
	
	def set_name(self,text):
		self.fit.name=str(text)
	
	def set_opt_tol(self,text):
		self.fit.opt_tol=float(str(text))
		
	def set_maxfun(self,text):
		self.fit.maxfun=int(str(text))
		
	def set_x0_k(self,text):
		self.fit.x0[0]=float(str(text))
		
	def set_x0_c0(self,text):
		self.fit.x0[1]=float(str(text))
		self.combo_x0_c0.setCurrentIndex(0)
		
	def set_x0_y0(self,text):
		self.fit.x0[2]=float(str(text))
		self.combo_x0_y0.setCurrentIndex(0)
		
	def set_UB_k(self,text):
		self.fit.UB_k=float(str(text))
		
	def set_UB_c0(self,text):
		self.fit.UB_cnaught=float(str(text))
	
	def set_UB_y0(self,text):
		self.fit.UB_ynaught=float(str(text))
		self.combo_UB_y0.setCurrentIndex(0)
		
	def set_LB_k(self,text):
		self.fit.LB_k=float(str(text))
		
	def set_LB_c0(self,text):
		self.fit.LB_cnaught=float(str(text))
	
	def set_LB_y0(self,text):
		self.fit.LB_ynaught=float(str(text))
		self.combo_LB_y0.setCurrentIndex(0)
		
	def sel_model(self,text):
		self.fit.model=str(text)
		self.power_vis()
			
	def set_npower(self,text):
		if str(text)!='':
			if float(str(text))<=1:
				print "Warning, n should be >1 to obtain a non-linear decay. If you want n=1, choose the exponential model."
			else:
				self.fit.npower=float(str(text))
			
	def sel_LB_y0(self,text):
		self.curr_LBmode=str(text)
		self.update_LB_y0()
		self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Checked)
	
	def power_vis(self):
		if self.fit.model=="power":
			self.qle_npower.setVisible(True)
			self.lbl_npower.setVisible(True)
		else: 	
			self.qle_npower.setVisible(False)
			self.lbl_npower.setVisible(False)
	
	def sel_meth(self,text):
		self.fit.opt_meth=str(text)
		self.update_bounds_after_meth(str(text))
	
	def update_bounds_after_meth(self,text):
	
		if text not in ["Constrained Nelder-Mead","TNC","L-BFGS-B","SLSQP","brute"]:
			
			self.temp_LB_k=self.fit.LB_k
			self.temp_UB_k=self.fit.UB_k
			self.temp_LB_cnaught=self.fit.LB_cnaught
			self.temp_UB_cnaught=self.fit.UB_cnaught
			self.temp_LB_ynaught=self.fit.LB_ynaught
			self.temp_UB_ynaught=self.fit.UB_ynaught
			
			self.check_bound_LB_k(0)
			self.check_bound_UB_k(0)
			self.check_bound_LB_c0(0)
			self.check_bound_UB_c0(0)
			self.check_bound_LB_y0(0)
			self.check_bound_UB_y0(0)
			
			self.cb_bound_LB_k.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_LB_c0.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_k.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_c0.setCheckState(QtCore.Qt.Unchecked)
			self.cb_bound_UB_y0.setCheckState(QtCore.Qt.Unchecked)		
			
		else:
			if self.temp_LB_k==[] or text=="brute":
				
				self.check_bound_LB_k(2)
				self.check_bound_UB_k(2)
				self.check_bound_LB_c0(2)
				self.check_bound_UB_c0(2)
				self.check_bound_LB_y0(2)
				self.check_bound_UB_y0(2)
				
				self.cb_bound_LB_k.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_LB_c0.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_k.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_c0.setCheckState(QtCore.Qt.Checked)
				self.cb_bound_UB_y0.setCheckState(QtCore.Qt.Checked)
				
				if text=="brute":
					self.bounds_for_brute()
				
			else: 
				if self.temp_LB_k==None:
					self.check_bound_LB_k(0)
					self.cb_bound_LB_k.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_k(2)
					self.cb_bound_LB_k.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_k==None:
					self.check_bound_UB_k(0)
					self.cb_bound_UB_k.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_k(2)
					self.cb_bound_UB_k.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_LB_cnaught==None:
					self.check_bound_LB_c0(0)
					self.cb_bound_LB_c0.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_c0(2)
					self.cb_bound_LB_c0.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_cnaught==None:
					self.check_bound_UB_c0(0)
					self.cb_bound_UB_c0.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_c0(2)
					self.cb_bound_UB_c0.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_LB_ynaught==None:
					self.check_bound_LB_y0(0)
					self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_LB_y0(2)
					self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Checked)
				
				if self.temp_UB_ynaught==None:
					self.check_bound_UB_y0(0)
					self.cb_bound_UB_y0.setCheckState(QtCore.Qt.Unchecked)
				else: 	
					self.check_bound_UB_c0(2)
					self.cb_bound_UB_y0.setCheckState(QtCore.Qt.Checked)
		
		return
		
	def bounds_for_brute(self):

		if self.fit.LB_k>=self.fit.UB_k:
			print "LBs need to be larger than UBs for k, going to fix this."
			if self.fit.UB_k>0:
				self.fit.UB_k==100*self.fit.LB_k
			else:
				self.fit.UB_k=0.01
			
			self.qle_UB_k.setText(str(self.fit.UB_k))
			
		if self.fit.LB_cnaught>=self.fit.UB_cnaught:
			print "LBs need to be larger than UBs for c0, going to fix this."
			self.fit.UB_cnaught=1.5*max(self.fit.embryo.ext_av_data_d)
			self.qle_UB_c0.setText(str(self.fit.UB_cnaught))
			
		if self.fit.LB_ynaught>=self.fit.UB_ynaught:
			print "LBs need to be larger than UBs for y0, going to fix this."
			self.sel_UB_y0("Max(I_int(t))")	
			self.qle_UB_y0.setText(str(self.fit.UB_ynaught))
			self.combo_UB_y0.setCurrentIndex(1)
		return
			
	def sel_LB_y0(self,text):
		self.curr_LBmode=str(text)
		self.update_LB_y0()
		self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Checked)
	
	def sel_UB_y0(self,text):
		
		self.cb_bound_UB_y0.setCheckState(QtCore.Qt.Checked)	
		if text=="Custom":
			pass
		if text=="Max(I_int(t))":
			self.fit.UB_ynaught=max(self.fit.embryo.int_av_data_d)
		self.qle_UB_y0.setText(str(self.fit.UB_ynaught))	
			
	def sel_x0_c0(self,text):

		self.update_x0_c0(text)
		
	def sel_x0_y0(self,text):

		self.update_x0_y0(text)
		
		
	def update_x0_c0(self,text):
			
		if text=="Custom":
			pass
		if text=="I_post(tstart)-I_pre":
			if self.fit.fit_ext==1:
				self.fit.x0[1]=self.fit.embryo.ext_av_data_d[0]-self.fit.embryo.pre.pre_ext
			elif self.fit.fit_int==1:
				self.fit.x0[1]=self.fit.embryo.int_av_data_d[0]-self.fit.embryo.pre.pre_int
			elif self.fit.fit_slice==1:
				self.fit.x0[1]=self.fit.embryo.slice_av_data_d[0]-self.fit.embryo.pre.pre_slice
		
		self.qle_x0_c0.setText(str(self.fit.x0[1]))
		
	def update_x0_y0(self,text):
		
		if text=="Custom":
			pass
		if text=="I_post(tstart)":
			if self.fit.fit_ext==1:
				self.fit.x0[2]=self.fit.embryo.ext_av_data_d[0]
			elif self.fit.fit_int==1:
				self.fit.x0[2]=self.fit.embryo.int_av_data_d[0]
			elif self.fit.fit_slice==1:
				self.fit.x0[2]=self.fit.embryo.slice_av_data_d[0]
		if text=="I_post(tend)":
			if self.fit.fit_ext==1:
				self.fit.x0[2]=self.fit.embryo.ext_av_data_d[-1]
			elif self.fit.fit_int==1:
				self.fit.x0[2]=self.fit.embryo.int_av_data_d[-1]
			elif self.fit.fit_slice==1:
				self.fit.x0[2]=self.fit.embryo.slice_av_data_d[-1]
		if text=="I_pre":
			if self.fit.fit_ext==1:
				self.fit.x0[2]=self.fit.embryo.pre.pre_ext
			elif self.fit.fit_int==1:
				self.fit.x0[2]=self.fit.embryo.pre.pre_int
			elif self.fit.fit_slice==1:
				self.fit.x0[2]=self.fit.embryo.pre.pre_slice
	
		self.qle_x0_y0.setText(str(self.fit.x0[2]))
			
	def update_LB_y0(self):
		
		text=self.curr_LBmode
		
		if text=="Custom":
			pass
		
		elif text=="Noise":
			self.fit.LB_ynaught=self.embryo.noise.noise
			
		elif text=="Bkgd":
		
			if self.fit.fit_ext==1:
				self.fit.LB_ynaught=self.molecule.bkgd_ext
			elif self.fit.fit_int==1:
				self.fit.LB_ynaught=self.molecule.bkgd_int
			elif self.fit.fit_slice==1:
				self.fit.LB_ynaught=self.molecule.bkgd_slice

		elif text=="Bkgd_pre":	
			if self.fit.fit_ext==1:
				self.fit.LB_ynaught=self.molecule.bkgd_pre_ext
			elif self.fit.fit_int==1:
				self.fit.LB_ynaught=self.molecule.bkgd_pre_int
			elif self.fit.fit_slice==1:
				self.fit.LB_ynaught=self.molecule.bkgd_pre_slice
		elif text=="F":
			self.embryo=pyfdap_fit.LB_ynaught_func(self.molecule,self.embryo,self.fit.fit_number)
			
		self.qle_LB_y0.setText(str(self.fit.LB_ynaught))
			
	def check_fit_ext(self, value):
		
		if value==2:
			self.fit.fit_ext=int(value/2)
			self.fit.fit_slice=0
			self.fit.fit_int=0
			self.cb_fit_slice.setCheckState(QtCore.Qt.Unchecked)	
			self.cb_fit_int.setCheckState(QtCore.Qt.Unchecked)	
			self.update_x0_c0(str(self.combo_x0_c0.currentText()))
			self.update_x0_y0(str(self.combo_x0_y0.currentText()))
			self.update_LB_y0()
		else:
			pass
		
	def check_fit_slice(self, value):
		
		if value==2:
			self.fit.fit_slice=int(value/2)
			self.fit.fit_ext=0
			self.fit.fit_int=0
			self.cb_fit_ext.setCheckState(QtCore.Qt.Unchecked)	
			self.cb_fit_int.setCheckState(QtCore.Qt.Unchecked)	
			self.update_x0_c0(str(self.combo_x0_c0.currentText()))
			self.update_x0_y0(str(self.combo_x0_y0.currentText()))
			self.update_LB_y0()
		else: 
			pass
		
	def check_fit_int(self, value):
		
		if value==2:
			self.fit.fit_int=int(value/2)
			self.fit.fit_slice=0
			self.fit.fit_ext=0
			self.cb_fit_slice.setCheckState(QtCore.Qt.Unchecked)	
			self.cb_fit_ext.setCheckState(QtCore.Qt.Unchecked)	
			self.update_x0_c0(str(self.combo_x0_c0.currentText()))
			self.update_x0_y0(str(self.combo_x0_y0.currentText()))
			self.update_LB_y0()
		else:
			pass
			
	def check_fit_c0(self, value):
		self.fit.fit_cnaught=int(value/2)
	
	def check_fit_y0(self, value):
		self.fit.fit_ynaught=int(value/2)
	
	def check_debug_fit(self, value):
		self.fit.debug_fit=int(value/2)
		
	def check_save_track(self, value):
		self.fit.save_track=int(value/2)	
	
	def check_bound_LB_k(self,value):
		
		if value==2:
			self.qle_LB_k.setReadOnly(False)
			self.fit.LB_k=float(self.qle_LB_k.text())
			
		elif value==0:	
			self.qle_LB_k.setReadOnly(True)
			self.fit.LB_k=None
			
	def check_bound_UB_k(self,value):
		
		if value==2:
			self.qle_UB_k.setReadOnly(False)
			self.fit.UB_k=float(self.qle_UB_k.text())
			
		elif value==0:	
			self.qle_UB_k.setReadOnly(True)
			self.fit.UB_k=None		
	
	def check_bound_LB_c0(self,value):
		
		if value==2:
			self.qle_LB_c0.setReadOnly(False)
			self.fit.LB_cnaught=float(self.qle_LB_c0.text())
			
		elif value==0:	
			self.qle_LB_c0.setReadOnly(True)
			self.fit.LB_cnaught=None
			
	def check_bound_UB_c0(self,value):
		
		if value==2:
			self.qle_UB_c0.setReadOnly(False)
			self.fit.UB_cnaught=float(self.qle_UB_c0.text())
			
		elif value==0:	
			self.qle_UB_c0.setReadOnly(True)
			self.fit.UB_cnaught=None	
			
	def check_bound_LB_y0(self,value):
		
		if value==2:
			self.qle_LB_y0.setReadOnly(False)
			self.fit.LB_ynaught=float(self.qle_LB_y0.text())
			
		elif value==0:	
			self.qle_LB_y0.setReadOnly(True)
			self.fit.LB_ynaught=None
			
	def check_bound_UB_y0(self,value):
		
		if value==2:
			self.qle_UB_y0.setReadOnly(False)
			self.fit.UB_ynaught=float(self.qle_UB_y0.text())
			
		elif value==0:	
			self.qle_UB_y0.setReadOnly(True)
			self.fit.UB_ynaught=None			
			
	def done_pressed(self):
		if self.fit.fit_ext==0 and self.fit.fit_int==0 and self.fit.fit_slice==0:
			self.fit.fit_ext=1
			
		self.done(1)
		return self.fit

#===================================================================================================================================
#Dialog to edit multiple fits
#===================================================================================================================================

class mult_fit_dialog(QtGui.QDialog):
	
	def __init__(self,molecule,sel_fits,parent):
		super(mult_fit_dialog,self).__init__(parent)
		
		self.sel_fits=sel_fits
		self.molecule=molecule
		
		self.setMinimumSize(600,300) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		
		#Button Actions
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lbl_opt_parms = QtGui.QLabel("Optimization Parameters", self)
		self.lbl_opt_parms.setFont(boldfont)
		self.lbl_guess = QtGui.QLabel("Initial Guess", self)
		self.lbl_guess.setFont(boldfont)
		self.lbl_bounds = QtGui.QLabel("Variable bounds", self)
		self.lbl_bounds.setFont(boldfont)
		self.lbl_bounded = QtGui.QLabel("Bounded?", self)
		self.lbl_bounded.setFont(boldfont)
		self.lbl_fitting = QtGui.QLabel("Fitting options", self)
		self.lbl_fitting.setFont(boldfont)
		
		self.lbl_opt_meth = QtGui.QLabel("opt_meth:", self)
		self.lbl_opt_tol = QtGui.QLabel("opt_tol", self)
		self.lbl_maxfun = QtGui.QLabel("maxfun:", self)
		self.lbl_save_track = QtGui.QLabel("save_track:", self)
		self.lbl_debug_fit = QtGui.QLabel("debug_fit:", self)
		
		self.lbl_x0_k = QtGui.QLabel("x0_k:", self)
		self.lbl_x0_c0 = QtGui.QLabel("x0_c0:", self)
		self.lbl_x0_y0 = QtGui.QLabel("x0_y0:", self)
		
		self.lbl_LB_k = QtGui.QLabel("LB_k:", self)
		self.lbl_UB_k = QtGui.QLabel("UB_k:", self)
		self.lbl_LB_y0 = QtGui.QLabel("LB_y0:", self)
		self.lbl_UB_y0 = QtGui.QLabel("UB_y0:", self)
		self.lbl_LB_c0 = QtGui.QLabel("LB_c0:", self)
		self.lbl_UB_c0 = QtGui.QLabel("UB_c0:", self)
		
		self.lbl_fit_ext = QtGui.QLabel("fit_ext:", self)
		self.lbl_fit_int = QtGui.QLabel("fit_int:", self)
		self.lbl_fit_slice = QtGui.QLabel("fit_slice:", self)
		self.lbl_fit_y0 = QtGui.QLabel("fit_y0:", self)
		self.lbl_fit_c0 = QtGui.QLabel("fit_c0:", self)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_meth = QtGui.QComboBox(self)
		self.combo_meth.addItem("Constrained Nelder-Mead")
		self.combo_meth.addItem("TNC")
		self.combo_meth.addItem("Nelder-Mead")
		self.combo_meth.addItem("Anneal")
		self.combo_meth.addItem("brute")
		self.combo_meth.addItem("COBYLA")
		self.combo_meth.addItem("BFGS")
		self.combo_meth.addItem("CG")
		self.combo_meth.activated[str].connect(self.sel_meth)   
		
		self.combo_x0_c0 = QtGui.QComboBox(self)
		self.combo_x0_c0.addItem("Custom")
		self.combo_x0_c0.addItem("I_post(tstart)-I_pre")
		self.combo_x0_c0.activated[str].connect(self.sel_x0_c0)   
		
		self.combo_x0_y0 = QtGui.QComboBox(self)
		self.combo_x0_y0.addItem("Custom")
		self.combo_x0_y0.addItem("I_post(tstart)")
		self.combo_x0_y0.addItem("I_post(tend)")
		self.combo_x0_y0.addItem("I_pre")
		self.combo_x0_y0.activated[str].connect(self.sel_x0_y0)   
		
		self.combo_UB_y0 = QtGui.QComboBox(self)
		self.combo_UB_y0.addItem("Custom")
		self.combo_UB_y0.addItem("Max(I_int(t))")
		self.combo_UB_y0.activated[str].connect(self.sel_UB_y0)   
		
		self.combo_LB_y0 = QtGui.QComboBox(self)
		self.combo_LB_y0.addItem("Custom")
		self.combo_LB_y0.addItem("Noise")
		self.combo_LB_y0.addItem("Bkgd_pre")
		self.combo_LB_y0.addItem("Bkgd")
		self.combo_LB_y0.addItem("F")
		self.curr_LBmode="Custom"
		self.combo_LB_y0.activated[str].connect(self.sel_LB_y0) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_debug_fit = QtGui.QCheckBox('', self)
		self.cb_save_track = QtGui.QCheckBox('', self)
		
		self.cb_bound_LB_k = QtGui.QCheckBox('', self)
		self.cb_bound_LB_c0 = QtGui.QCheckBox('', self)
		self.cb_bound_LB_y0 = QtGui.QCheckBox('', self)
		
		self.cb_bound_UB_k = QtGui.QCheckBox('', self)
		self.cb_bound_UB_c0 = QtGui.QCheckBox('', self)
		self.cb_bound_UB_y0 = QtGui.QCheckBox('', self)
		
		self.cb_fit_ext = QtGui.QCheckBox('', self)
		self.cb_fit_int = QtGui.QCheckBox('', self)
		self.cb_fit_slice = QtGui.QCheckBox('', self)
		
		self.cb_fit_c0 = QtGui.QCheckBox('', self)
		self.cb_fit_y0 = QtGui.QCheckBox('', self)
		
		self.cb_debug_fit.setCheckState(2*self.sel_fits[0].debug_fit)	
		self.cb_save_track.setCheckState(2*self.sel_fits[0].save_track)	
	
		if self.sel_fits[0].LB_k==None:
			self.cb_bound_LB_k.setCheckState(0)
		else: 
			self.cb_bound_LB_k.setCheckState(2)
		
		if self.sel_fits[0].UB_k==None:
			self.cb_bound_UB_k.setCheckState(0)
		else: 
			self.cb_bound_UB_k.setCheckState(2)
		
		if self.sel_fits[0].LB_cnaught==None:
			self.cb_bound_LB_c0.setCheckState(0)
		else: 
			self.cb_bound_LB_c0.setCheckState(2)
		
		if self.sel_fits[0].UB_cnaught==None:
			self.cb_bound_UB_c0.setCheckState(0)
		else: 
			self.cb_bound_UB_c0.setCheckState(2)
		
		if self.sel_fits[0].LB_ynaught==None:
			self.cb_bound_LB_y0.setCheckState(0)
		else: 
			self.cb_bound_LB_y0.setCheckState(2)
		
		if self.sel_fits[0].UB_ynaught==None:
			self.cb_bound_UB_y0.setCheckState(0)
		else: 
			self.cb_bound_UB_y0.setCheckState(2)
		
		self.connect(self.cb_bound_LB_k, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_k)
		self.connect(self.cb_bound_LB_c0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_c0)
		self.connect(self.cb_bound_LB_y0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_LB_y0)
		self.connect(self.cb_bound_UB_k, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_k)
		self.connect(self.cb_bound_UB_c0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_c0)
		self.connect(self.cb_bound_UB_y0, QtCore.SIGNAL('stateChanged(int)'), self.check_bound_UB_y0)
		
		self.cb_fit_ext.setCheckState(2*self.sel_fits[0].fit_ext)
		self.cb_fit_int.setCheckState(2*self.sel_fits[0].fit_int)
		self.cb_fit_slice.setCheckState(2*self.sel_fits[0].fit_slice)
		self.cb_fit_c0.setCheckState(2*self.sel_fits[0].fit_cnaught)
		self.cb_fit_y0.setCheckState(2*self.sel_fits[0].fit_ynaught)
		
		self.connect(self.cb_debug_fit, QtCore.SIGNAL('stateChanged(int)'), self.check_debug_fit)
		self.connect(self.cb_save_track, QtCore.SIGNAL('stateChanged(int)'), self.check_save_track)
		
		self.connect(self.cb_fit_ext, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_ext)
		self.connect(self.cb_fit_int, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_int)
		self.connect(self.cb_fit_slice, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_slice)
		self.connect(self.cb_fit_c0, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_c0)
		self.connect(self.cb_fit_y0, QtCore.SIGNAL('stateChanged(int)'), self.check_fit_y0)
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_opt_tol = QtGui.QLineEdit(str(self.sel_fits[0].opt_tol))
		self.qle_maxfun = QtGui.QLineEdit(str(self.sel_fits[0].maxfun))
		
		self.qle_x0_k = QtGui.QLineEdit("")
		self.qle_x0_c0 = QtGui.QLineEdit("")
		self.qle_x0_y0 = QtGui.QLineEdit("")
		
		if self.sel_fits[0].LB_k==None:
			self.qle_LB_k = QtGui.QLineEdit("0.0")
			self.qle_LB_k.setReadOnly(True)
		else:	
			self.qle_LB_k = QtGui.QLineEdit(str(self.sel_fits[0].LB_k))
			
		if self.sel_fits[0].UB_k==None:
			self.qle_UB_k = QtGui.QLineEdit("0.0")
			self.qle_UB_k.setReadOnly(True)
		else:
			self.qle_UB_k = QtGui.QLineEdit(str(self.sel_fits[0].UB_k))
			
		if self.sel_fits[0].LB_cnaught==None:
			self.qle_LB_c0 = QtGui.QLineEdit("0.0")
			self.qle_LB_c0.setReadOnly(True)
		else:
			self.qle_LB_c0 = QtGui.QLineEdit(str(self.sel_fits[0].LB_cnaught))
			
		if self.sel_fits[0].UB_cnaught==None:
			self.qle_UB_c0 = QtGui.QLineEdit("0.0")
			self.qle_UB_c0.setReadOnly(True)
		else:
			self.qle_UB_c0 = QtGui.QLineEdit(str(self.sel_fits[0].UB_cnaught))
			
		if self.sel_fits[0].LB_ynaught==None:
			self.qle_LB_y0 = QtGui.QLineEdit("0.0")
			self.qle_LB_y0.setReadOnly(True)
		else:
			self.qle_LB_y0 = QtGui.QLineEdit(str(self.sel_fits[0].LB_ynaught))
			
		if self.sel_fits[0].UB_ynaught==None:			
			self.qle_UB_y0 = QtGui.QLineEdit("0.0")
			self.qle_UB_y0.setReadOnly(True)
		else:
			self.qle_UB_y0 = QtGui.QLineEdit(str(self.sel_fits[0].UB_ynaught))
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_LB_k.setValidator(self.double_valid)
		self.qle_LB_c0.setValidator(self.double_valid)
		self.qle_LB_y0.setValidator(self.double_valid)
		self.qle_UB_k.setValidator(self.double_valid)
		self.qle_UB_c0.setValidator(self.double_valid)
		self.qle_UB_y0.setValidator(self.double_valid)
		
		self.qle_opt_tol.textChanged[str].connect(self.set_opt_tol)
		self.qle_maxfun.textChanged[str].connect(self.set_maxfun)
		
		self.qle_x0_k.textChanged[str].connect(self.set_x0_k)
		self.qle_x0_c0.textChanged[str].connect(self.set_x0_c0)
		self.qle_x0_y0.textChanged[str].connect(self.set_x0_y0)
		
		self.qle_UB_k.textChanged[str].connect(self.set_UB_k)
		self.qle_UB_c0.textChanged[str].connect(self.set_UB_c0)
		self.qle_UB_y0.textChanged[str].connect(self.set_UB_y0)
		
		self.qle_LB_k.textChanged[str].connect(self.set_LB_k)
		self.qle_LB_c0.textChanged[str].connect(self.set_LB_c0)
		self.qle_LB_y0.textChanged[str].connect(self.set_LB_y0)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_opt_parms,0,1,1,2,Qt.AlignHCenter)
		grid.addWidget(self.lbl_guess,0,3,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounds,0,6,1,3,Qt.AlignHCenter)
		grid.addWidget(self.lbl_bounded,0,9,1,1,Qt.AlignHCenter)
		grid.addWidget(self.lbl_fitting,0,10,1,2,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_opt_meth,2,1)
		grid.addWidget(self.lbl_opt_tol,3,1)
		grid.addWidget(self.lbl_maxfun,4,1)
		grid.addWidget(self.lbl_debug_fit,5,1)
		grid.addWidget(self.lbl_save_track,6,1)
		
		grid.addWidget(self.combo_meth,2,2)
		grid.addWidget(self.qle_opt_tol,3,2)
		grid.addWidget(self.qle_maxfun,4,2)
		grid.addWidget(self.cb_debug_fit,5,2)
		grid.addWidget(self.cb_save_track,6,2)
		
		grid.addWidget(self.lbl_x0_k,1,3)
		grid.addWidget(self.lbl_x0_c0,2,3)
		grid.addWidget(self.lbl_x0_y0,3,3)
		
		grid.addWidget(self.qle_x0_k,1,4)
		grid.addWidget(self.qle_x0_c0,2,4)
		grid.addWidget(self.qle_x0_y0,3,4)
		
		grid.addWidget(self.combo_x0_c0,2,5)
		grid.addWidget(self.combo_x0_y0,3,5)
		
		grid.addWidget(self.lbl_LB_k,1,6)
		grid.addWidget(self.lbl_LB_c0,3,6)
		grid.addWidget(self.lbl_LB_y0,5,6)
		
		grid.addWidget(self.lbl_UB_k,2,6)
		grid.addWidget(self.lbl_UB_c0,4,6)
		grid.addWidget(self.lbl_UB_y0,6,6)
		
		grid.addWidget(self.qle_LB_k,1,7)
		grid.addWidget(self.qle_LB_c0,3,7)
		grid.addWidget(self.qle_LB_y0,5,7)
		
		grid.addWidget(self.qle_UB_k,2,7)
		grid.addWidget(self.qle_UB_c0,4,7)
		grid.addWidget(self.qle_UB_y0,6,7)
		
		grid.addWidget(self.combo_LB_y0,5,8)
		grid.addWidget(self.combo_UB_y0,6,8)
		
		grid.addWidget(self.cb_bound_LB_k,1,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_c0,3,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_LB_y0,5,9,Qt.AlignHCenter)
		
		grid.addWidget(self.cb_bound_UB_k,2,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_c0,4,9,Qt.AlignHCenter)
		grid.addWidget(self.cb_bound_UB_y0,6,9,Qt.AlignHCenter)
		
		grid.addWidget(self.lbl_fit_ext,1,10)
		grid.addWidget(self.lbl_fit_int,2,10)
		grid.addWidget(self.lbl_fit_slice,3,10)
		grid.addWidget(self.lbl_fit_c0,4,10)
		grid.addWidget(self.lbl_fit_y0,5,10)
		
		grid.addWidget(self.cb_fit_ext,1,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_int,2,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_slice,3,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_c0,4,11,Qt.AlignHCenter)
		grid.addWidget(self.cb_fit_y0,5,11,Qt.AlignHCenter)
		
		grid.addWidget(self.btn_done,8,11)
		
		grid.setColumnStretch(0,1)
		grid.setColumnStretch(12,1)
	
		grid.setRowStretch(9,1)
		
		self.setLayout(grid)    
		self.setWindowTitle('Edit multiple Fits')    
		self.show()
		
	def set_opt_tol(self,text):
		for fit in self.sel_fits:
			fit.opt_tol=float(str(text))
		
	def set_maxfun(self,text):
		for fit in self.sel_fits:
			fit.maxfun=int(str(text))
		
	def set_x0_k(self,text):
		for fit in self.sel_fits:
			fit.x0[0]=float(str(text))
		
	def set_x0_c0(self,text):
		for fit in self.sel_fits:
			fit.x0[1]=float(str(text))
		self.combo_x0_c0.setCurrentIndex(0)
		
	def set_x0_y0(self,text):
		for fit in self.sel_fits:
			fit.x0[2]=float(str(text))
		self.combo_x0_y0.setCurrentIndex(0)
		
	def set_UB_k(self,text):
		for fit in self.sel_fits:	
			fit.UB_k=float(str(text))
		
	def set_UB_c0(self,text):
		for fit in self.sel_fits:
			fit.UB_cnaught=float(str(text))
	
	def set_UB_y0(self,text):
		for fit in self.sel_fits:
			fit.UB_ynaught=float(str(text))
		self.combo_UB_y0.setCurrentIndex(0)
		
	def set_LB_k(self,text):
		for fit in self.sel_fits:
			fit.LB_k=float(str(text))
		
	def set_LB_c0(self,text):
		for fit in self.sel_fits:
			fit.LB_cnaught=float(str(text))
	
	def set_LB_y0(self,text):
		for fit in self.sel_fits:
			fit.LB_ynaught=float(str(text))
		self.combo_LB_y0.setCurrentIndex(0)
		
	def sel_meth(self,text):
		for fit in self.sel_fits:
			fit.opt_meth=str(text)
	
	def sel_LB_y0(self,text):
		for fit in self.sel_fits:
			self.curr_LBmode=str(text)
			self.update_LB_y0(fit)
		self.cb_bound_LB_y0.setCheckState(QtCore.Qt.Checked)
	
	def sel_UB_y0(self,text):
		
		self.cb_bound_UB_y0.setCheckState(QtCore.Qt.Checked)	
		for fit in self.sel_fits:
			self.update_UB_y0(text,fit)
			
	def update_UB_y0(self,text,fit):
	
		if text=="Custom":
			fit.UB_ynaught=(float(self.qle_LB_y0.text()))
		if text=="Max(I_int(t))":
			fit.UB_ynaught=max(fit.embryo.int_av_data_d)
				
	def sel_x0_c0(self,text):
		for fit in self.sel_fits:
			self.update_x0_c0(text,fit)
		
	def sel_x0_y0(self,text):
		for fit in self.sel_fits:
			self.update_x0_y0(text,fit)
		
	def update_x0_c0(self,text,fit):
		
		if text=="Custom":
			pass
		if text=="I_post(tstart)-I_pre":
			if fit.fit_ext==1:
				fit.x0[1]=fit.embryo.ext_av_data_d[0]-fit.embryo.pre.pre_ext
			elif fit.fit_int==1:
				fit.x0[1]=fit.embryo.int_av_data_d[0]-fit.embryo.pre.pre_int
			elif fit.fit_slice==1:
				fit.x0[1]=fit.embryo.slice_av_data_d[0]-fit.embryo.pre.pre_slice
			
	def update_x0_y0(self,text,fit):
		
		if text=="Custom":
			pass
		if text=="I_post(tstart)":
			if fit.fit_ext==1:
				fit.x0[2]=fit.embryo.ext_av_data_d[0]
			elif fit.fit_int==1:
				fit.x0[2]=fit.embryo.int_av_data_d[0]
			elif fit.fit_slice==1:
				fit.x0[2]=fit.embryo.slice_av_data_d[0]
		if text=="I_post(tend)":
			if fit.fit_ext==1:
				fit.x0[2]=fit.embryo.ext_av_data_d[-1]
			elif fit.fit_int==1:
				fit.x0[2]=fit.embryo.int_av_data_d[-1]
			elif fit.fit_slice==1:
				fit.x0[2]=fit.embryo.slice_av_data_d[-1]	
		if text=="I_pre":
			if fit.fit_ext==1:
				fit.x0[2]=fit.embryo.pre.pre_ext
			elif fit.fit_int==1:
				fit.x0[2]=fit.embryo.pre.pre_int
			elif fit.fit_slice==1:
				fit.x0[2]=fit.embryo.pre.pre_slice
				
			
	def update_LB_y0(self,fit):
		
		text=self.curr_LBmode
		
		if text=="Custom":
			fit.LB_ynaught=(float(self.qle_LB_y0.text()))
		
		elif text=="Noise":
			fit.LB_ynaught=self.embryo.noise.noise
			
		elif text=="Bkgd":
		
			if fit.fit_ext==1:
				fit.LB_ynaught=self.molecule.bkgd_ext
			elif fit.fit_int==1:
				fit.LB_ynaught=self.molecule.bkgd_int
			elif fit.fit_slice==1:
				fit.LB_ynaught=self.molecule.bkgd_slice
		
		elif text=="Bkgd_pre":	
			if fit.fit_ext==1:
				fit.LB_ynaught=self.molecule.bkgd_pre_ext
			elif fit.fit_int==1:
				fit.LB_ynaught=self.molecule.bkgd_pre_int
			elif fit.fit_slice==1:
				fit.LB_ynaught=self.molecule.bkgd_pre_slice
		elif text=="F":
			self.embryo=pyfdap_fit.LB_ynaught_func(self.molecule,fit.embryo,fit.fit_number)
	
	def update_all_bounds(self,fit):
		fit.LB_k=(float(self.qle_LB_k.text()))
		fit.UB_k=(float(self.qle_UB_k.text()))
		fit.LB_cnaught=(float(self.qle_LB_c0.text()))
		fit.UB_cnaught=(float(self.qle_UB_c0.text()))
		self.update_LB_y0(fit)
		self.update_UB_y0(str(self.combo_UB_y0.currentText()),fit)
		
	def check_fit_ext(self, value):
		
		if value==2:
			for fit in self.sel_fits:
				fit.fit_ext=int(value/2)
				fit.fit_slice=0
				fit.fit_int=0
				self.update_LB_y0(fit)
				
			self.cb_fit_slice.setCheckState(QtCore.Qt.Unchecked)	
			self.cb_fit_int.setCheckState(QtCore.Qt.Unchecked)	
			
		else:
			pass
		
	def check_fit_slice(self, value):
		
		if value==2:
			for fit in self.sel_fits:
				fit.fit_slice=int(value/2)
				fit.fit_ext=0
				fit.fit_int=0
				self.update_LB_y0(fit)
				
			self.cb_fit_ext.setCheckState(QtCore.Qt.Unchecked)	
			self.cb_fit_int.setCheckState(QtCore.Qt.Unchecked)	
			
		else: 
			pass
		
	def check_fit_int(self, value):
		
		if value==2:
			for fit in self.sel_fits:
				fit.fit_int=int(value/2)
				fit.fit_slice=0
				fit.fit_ext=0
				self.update_LB_y0(fit)
				
			self.cb_fit_slice.setCheckState(QtCore.Qt.Unchecked)	
			self.cb_fit_ext.setCheckState(QtCore.Qt.Unchecked)	
			
		else:
			pass
			
	def check_fit_c0(self, value):
		for fit in self.sel_fits:
			fit.fit_cnaught=int(value/2)
	
	def check_fit_y0(self, value):
		for fit in self.sel_fits:
			fit.fit_ynaught=int(value/2)
	
	def check_debug_fit(self, value):
		for fit in self.sel_fits:
			fit.debug_fit=int(value/2)
		
	def check_save_track(self, value):
		for fit in self.sel_fits:
			fit.save_track=int(value/2)	
	
	def check_bound_LB_k(self,value):
		
		if value==2:
			self.qle_LB_k.setReadOnly(False)
			for fit in self.sel_fits:
				self.update_all_bounds(fit)		
		elif value==0:	
			self.qle_LB_k.setReadOnly(True)
			for fit in self.sel_fits:
				fit.LB_k=None
			
	def check_bound_UB_k(self,value):
		
		if value==2:
			self.qle_UB_k.setReadOnly(False)
			for fit in self.sel_fits:
				self.update_all_bounds(fit)		
		elif value==0:	
			self.qle_UB_k.setReadOnly(True)
			for fit in self.sel_fits:
				fit.UB_k=None		
	
	def check_bound_LB_c0(self,value):
		
		if value==2:
			self.qle_LB_c0.setReadOnly(False)
			for fit in self.sel_fits:
				self.update_all_bounds(fit)	
		elif value==0:	
			self.qle_LB_c0.setReadOnly(True)
			for fit in self.sel_fits:
				fit.LB_cnaught=None
			
	def check_bound_UB_c0(self,value):
		
		if value==2:
			self.qle_UB_c0.setReadOnly(False)	
			for fit in self.sel_fits:
				self.update_all_bounds(fit)	
		elif value==0:	
			self.qle_UB_c0.setReadOnly(True)
			for fit in self.sel_fits:
				fit.UB_cnaught=None	
			
	def check_bound_LB_y0(self,value):
		
		if value==2:
			self.qle_LB_y0.setReadOnly(False)
			for fit in self.sel_fits:
				self.update_all_bounds(fit)	
		elif value==0:		
			self.qle_LB_y0.setReadOnly(True)
			for fit in self.sel_fits:
				fit.LB_ynaught=None
			
	def check_bound_UB_y0(self,value):
		
		if value==2:
			self.qle_UB_y0.setReadOnly(False)
			for fit in self.sel_fits:
				self.update_all_bounds(fit)	
		elif value==0:	
			self.qle_UB_y0.setReadOnly(True)
			for fit in self.sel_fits:
				fit.UB_ynaught=None			
			
	def done_pressed(self):
	
		self.done(1)
	
#===================================================================================================================================
#Dialog for About
#===================================================================================================================================

class about_dialog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(about_dialog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("This is PyFDAP version "+ parent.version, self)
		self.lbl_author = QtGui.QLabel("Author: Alexander Blaessle", self)
		self.lbl_website = QtGui.QLabel("Website: <a href="+parent.website+">"+parent.website+"</a>", self)
		self.connect(self.lbl_website, SIGNAL("linkActivated(QString)"), self.OpenURL) 

		self.btn_cancel=QtGui.QPushButton('Close')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.lbl_author)
		self.vbox.addWidget(self.lbl_website)		
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def OpenURL(self, URL): 
		QtGui.QDesktopServices().openUrl(QUrl(URL)) 
	
	def cancel(self):
	
		self.done(1)
					
#===================================================================================================================================
#Dialog for anaylze all progress
#===================================================================================================================================

class analyze_all_prog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(analyze_all_prog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("Analyzing complete molecule ...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel analysis')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_analysis)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def cancel_analysis(self):
	
		self.accepted.emit()
	
class analyze_all_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, molecule=None, parent=None):
		QtCore.QThread.__init__(self)
		self.molecule=molecule
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.molecule==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			for embryo in self.molecule.embryos:
				embryo=pyfdap_img.analyze_fdap_data(embryo)
			
			if len(self.molecule.bkgds)>0:
				self.molecule=pyfdap_img.analyze_fdap_bkgds(self.molecule)
				
			self.taskFinished.emit()

#===================================================================================================================================
#Dialog for anaylze progress
#===================================================================================================================================

class analyze_prog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(analyze_prog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("Data analysis in progress...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel analysis')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_analysis)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def cancel_analysis(self):
	
		self.accepted.emit()
		
class analyze_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, embryo=None, parent=None):
		QtCore.QThread.__init__(self)
		self.embryo=embryo
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.embryo==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			self.embryo=pyfdap_img.analyze_fdap_data(self.embryo)
			self.taskFinished.emit()
			
#===================================================================================================================================
#Dialog for bkgd analyze
#===================================================================================================================================

class analyze_bkgd_prog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(analyze_bkgd_prog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("Background analysis in progress...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel analysis')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_analysis)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def cancel_analysis(self):
	
		self.accepted.emit()
	
class analyze_bkgd_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, molecule=None, parent=None):
		QtCore.QThread.__init__(self)
		self.molecule=molecule
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.molecule==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			self.molecule=pyfdap_img.analyze_fdap_bkgds(self.molecule)
			self.taskFinished.emit()
			
#===================================================================================================================================
#Dialog for fitting progress
#===================================================================================================================================

class fitting_prog(QtGui.QDialog):
	
	def __init__(self,parent):
		super(fitting_prog,self).__init__(parent)
		self.lbl_name = QtGui.QLabel("Fitting in progress...", self)
		self.btn_cancel=QtGui.QPushButton('Cancel fitting')
		self.btn_cancel.connect(self.btn_cancel, QtCore.SIGNAL('clicked()'), self.cancel_fitting)	
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lbl_name)
		self.vbox.addWidget(self.btn_cancel)
		
		self.setLayout(self.vbox)
		self.show()	
	
	def cancel_fitting(self):
	
		self.accepted.emit()

class fitting_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, embryo=None, fit=None, gui=None, parent=None):
		QtCore.QThread.__init__(self)
		self.embryo=embryo
		self.fit=fit
		self.gui=gui
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.embryo==None or self.fit==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			self.embryo=pyfdap_fit.fdap_fitting(self.embryo,self.fit.fit_number,gui=self.gui)
			self.taskFinished.emit()

class fitting_all_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, embryo=None, fits=None, gui=None, parent=None):
		QtCore.QThread.__init__(self)
		self.embryo=embryo
		self.fits=fits
		self.gui=gui
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.embryo==None or self.fits==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			for fit in self.fits:
				self.embryo=pyfdap_fit.fdap_fitting(self.embryo,fit.fit_number,gui=self.gui)
				print "Fitted", fit.name
			self.taskFinished.emit()
			
class fitting_mol_thread(QtCore.QThread):
	taskFinished = QtCore.pyqtSignal()
    
	def __init__(self, molecule=None, gui=None, parent=None):
		QtCore.QThread.__init__(self)
		self.molecule=molecule
		self.gui=gui
		
	def __del__(self):
		self.wait()
    
	def run(self):
		
		if self.molecule==None:
			self.terminate()
			self.taskFinished.emit() 	
		else:
			for embryo in self.molecule.embryos:
				print "Fitting all fits of embryo", embryo.name
				for fit in embryo.fits:
					self.embryo=pyfdap_fit.fdap_fitting(embryo,fit.fit_number,gui=self.gui)	
					print "Fitted", fit.name
			
			self.taskFinished.emit()			

#===================================================================================================================================
#Dialog for selecting fits for averaging molecule
#===================================================================================================================================

class select_fits(QtGui.QDialog):
	
	def __init__(self,molecule,single_fit,parent):
		super(select_fits,self).__init__(parent)
		
		self.molecule=molecule
		self.embr_in_right_list=[]
		self.fit_in_right_list=[]
		self.single_fit=single_fit
		
		self.btn_add=QtGui.QToolButton()
		self.btn_add.connect(self.btn_add, QtCore.SIGNAL('clicked()'), self.add_fit)
		self.btn_add.setArrowType(QtCore.Qt.RightArrow)
		
		self.btn_remove=QtGui.QToolButton()
		self.btn_remove.connect(self.btn_remove, QtCore.SIGNAL('clicked()'), self.remove_fit)
		self.btn_remove.setArrowType(QtCore.Qt.LeftArrow)
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		self.left_list=QtGui.QTreeWidget()
		self.left_list.setHeaderLabels(["Available Fits"])
		self.left_list.setColumnWidth(0,200)
		self.left_list.setColumnWidth(1,75)
		self.left_list.itemDoubleClicked.connect(self.add_fit)
				
		self.right_list=QtGui.QTreeWidget()
		self.right_list.setHeaderLabels(["Selected Fits"])
		self.right_list.setColumnWidth(0,200)
		self.right_list.setColumnWidth(1,75)
		self.right_list.itemDoubleClicked.connect(self.remove_fit)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.btn_add)
		self.vbox.addWidget(self.btn_remove)
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.left_list)
		self.hbox.addLayout(self.vbox)
		self.hbox.addWidget(self.right_list)
		
		self.vbox2 = QtGui.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btn_done)
		
		self.init_left_list()
		self.resize(400,500)
		self.setLayout(self.vbox2)
		self.show()	
	
	def init_left_list(self):
		
		for emb in self.molecule.embryos:
			
			if self.single_fit==1:
			
				#Add to left bar
				fitted=0
				if shape(emb.fits)>0:
					for fit in emb.fits:
						if shape(fit.fit_av_d)[0]>1:	
							fitted=1
				
				#Adding embryo to sidebar
				if fitted==1:
					self.curr_embr_node=QtGui.QTreeWidgetItem(self.left_list,[emb.name])
			
			else:
				self.curr_embr_node=QtGui.QTreeWidgetItem(self.left_list,[emb.name])
			
			
			#Add fits if they exist
			for fit in emb.fits:
				if shape(fit.fit_av_d)[0]>0:		
					QtGui.QTreeWidgetItem(self.curr_embr_node,[fit.name])
			
			if self.single_fit==1:
				if fitted==1:
					self.left_list.expandItem(self.curr_embr_node)
			else:
				self.left_list.expandItem(self.curr_embr_node)

	def add_fit(self):
		self.get_left_selections()
		if self.left_list.currentItem()==None or self.left_list.currentItem().parent()==None:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			#Check if current embryo is already in right list
			if self.curr_embr in self.embr_in_right_list:
				#Check if fit is already in right list
				if self.curr_fit in self.fit_in_right_list:
					pass
				else:
					#Check if target embr already has a fit; if so, remove it
					self.get_right_ind(embryo_name=self.curr_embr.name)
					if self.single_fit==1:
						if self.curr_target_embr_node.childCount()>0:
							self.get_right_ind(fit_name=self.curr_target_embr_node.child(0).data(0,0).toString())
							self.curr_target_embr_node.takeChild(0)
							self.fit_in_right_list.remove(self.curr_target_fit)
							
					#Add to list of fits on right side	
					self.fit_in_right_list.append(self.curr_fit)
					new_node=QtGui.QTreeWidgetItem(self.curr_target_embr_node,[self.curr_fit.name])
					self.right_list.expandItem(self.curr_target_embr_node)
				
			else:
				#Add curr_embr to list
				self.embr_in_right_list.append(self.curr_embr)
				new_node=QtGui.QTreeWidgetItem(self.right_list,[self.curr_embr.name])
				
				#Add to list of fits on right side	
				self.fit_in_right_list.append(self.curr_fit)
				QtGui.QTreeWidgetItem(new_node,[self.curr_fit.name])
				self.right_list.expandItem(new_node)
		

	def remove_fit(self):
		self.get_right_selections()
		if self.right_list.currentItem()==None or self.right_list.currentItem().parent()==None:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			ind=self.curr_embr_node.indexOfChild(self.right_list.currentItem())
			self.curr_embr_node.takeChild(ind)
			self.fit_in_right_list.remove(self.curr_fit)
		
	
	def get_left_selections(self):

		#Get current embryo selection
		for emb in self.molecule.embryos:
			if self.left_list.currentItem().parent()==None: 
				#This is an embryo
				if self.left_list.currentItem().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.left_list.currentItem()
					self.curr_fit=None
					self.curr_fit_node=None
					break
			else:
				#This is a fit
				if self.left_list.currentItem().parent().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.left_list.currentItem().parent()
					for fit in self.curr_embr.fits:
						if fit.name==self.left_list.currentItem().data(0,0).toString():
							self.curr_fit=fit
							self.curr_fit_node=self.left_list.currentItem()
							break
		
	def get_right_selections(self):

		#Get current embryo selection
		for emb in self.molecule.embryos:
			if self.right_list.currentItem().parent()==None: 
				#This is an embryo
				if self.right_list.currentItem().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.right_list.currentItem()
					self.curr_fit=None
					self.curr_fit_node=None
					break
			else:
				#This is a fit
				if self.right_list.currentItem().parent().data(0,0).toString()==emb.name:
					self.curr_embr=emb
					self.curr_embr_node=self.right_list.currentItem().parent()
					for fit in self.curr_embr.fits:
						if fit.name==self.right_list.currentItem().data(0,0).toString():
							self.curr_fit=fit
							self.curr_fit_node=self.right_list.currentItem()
							break	
						
	def get_left_ind(self,embryo_name=None,fit_name=None):	
		if embryo_name!=None:
			self.curr_target_embr_node=self.left_list.findItems(embryo_name,Qt.MatchExactly,0)
			self.curr_target_embr_node=self.curr_target_embr_node[0]
		if fit_name!=None:
			self.curr_target_fit_node=self.left_list.findItems(fit_name,Qt.MatchExactly,0)
			for fit in self.curr_embr.fits:
				if fit.name==fit_name:
					self.curr_target_fit=fit
				
	def get_right_ind(self,embryo_name=None,fit_name=None):	
		if embryo_name!=None:
			self.curr_target_embr_node=self.right_list.findItems(embryo_name,Qt.MatchExactly,0)
			self.curr_target_embr_node=self.curr_target_embr_node[0]
		if fit_name!=None:
			self.curr_target_fit_node=self.right_list.findItems(fit_name,Qt.MatchExactly,0)
			for fit in self.curr_embr.fits:
				if fit.name==fit_name:
					self.curr_target_fit=fit
					
	def done_pressed(self):
		self.molecule.sel_fits=self.fit_in_right_list
		
			
		self.done(1)	
			
#===================================================================================================================================
#Dialog for selecting fits for averaging molecule
#===================================================================================================================================

class select_ignored_frames(QtGui.QDialog):
	
	def __init__(self,obj,parent):
		super(select_ignored_frames,self).__init__(parent)
		
		self.obj=obj
		
		self.btn_add=QtGui.QToolButton()
		self.btn_add.connect(self.btn_add, QtCore.SIGNAL('clicked()'), self.add_frame)
		self.btn_add.setArrowType(QtCore.Qt.RightArrow)
		
		self.btn_remove=QtGui.QToolButton()
		self.btn_remove.connect(self.btn_remove, QtCore.SIGNAL('clicked()'), self.remove_frame)
		self.btn_remove.setArrowType(QtCore.Qt.LeftArrow)
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		self.left_list=QtGui.QTreeWidget()
		self.left_list.setHeaderLabels(["Available Frames"])
		self.left_list.setColumnWidth(0,200)
		self.left_list.setColumnWidth(1,75)
		self.left_list.itemDoubleClicked.connect(self.add_frame)
			
		self.right_list=QtGui.QTreeWidget()
		self.right_list.setHeaderLabels(["Ignored Frames"])
		self.right_list.setColumnWidth(0,200)
		self.right_list.setColumnWidth(1,75)
		self.right_list.itemDoubleClicked.connect(self.remove_frame)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.btn_add)
		self.vbox.addWidget(self.btn_remove)
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.left_list)
		self.hbox.addLayout(self.vbox)
		self.hbox.addWidget(self.right_list)
		
		self.vbox2 = QtGui.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btn_done)
		
		self.init_left_list()
		self.resize(400,500)
		self.setLayout(self.vbox2)
		self.show()	
	
	def append_nulls(self,i):
		dec=len(str(len(self.obj.centers_embr_px)))
		nullstr=''	
		for j in range(dec-len(str(i))):
			nullstr=nullstr+'0'
			
		return 	nullstr+str(i)
	
	def init_left_list(self):
		
		#Collecting ignored frames
		if hasattr(self.obj,'ignored'):
			self.obj.ignored=list(self.obj.ignored)
		else:
			self.obj.ignored=[]
		
		for i in range(len(self.obj.centers_embr_px)):
			
			if i in self.obj.ignored:
				QtGui.QTreeWidgetItem(self.right_list,[self.append_nulls(i)])
			else:
				QtGui.QTreeWidgetItem(self.left_list,[self.append_nulls(i)])
			
	def add_frame(self):
		self.get_left_selections()
		self.left_list.takeTopLevelItem(self.curr_left_ind)
		
		#Insert new node at the right place
		new_node=QtGui.QTreeWidgetItem(self.right_list,[self.append_nulls(self.curr_frame)])
		self.right_list.sortItems(0,Qt.AscendingOrder)
		
		self.obj.ignored.append(self.curr_frame)			

	def remove_frame(self):
		self.get_right_selections()
		self.right_list.takeTopLevelItem(self.curr_right_ind)
		
		#Insert new node at the right place
		new_node=QtGui.QTreeWidgetItem(self.left_list,[self.append_nulls(self.curr_frame)])
		self.left_list.sortItems(0,Qt.AscendingOrder)
		
		self.obj.ignored.remove(self.curr_frame)
		
	def get_left_selections(self):
			
		#Get current frame selection
		for i in range(len(self.obj.centers_embr_px)):
			if int(self.left_list.currentItem().data(0,0).toString())==i:
				self.curr_frame=i
				self.curr_left_ind=self.left_list.indexFromItem(self.left_list.currentItem())
				self.curr_left_ind=self.curr_left_ind.row()
				break
		
	def get_right_selections(self):

		#Get current frame selection
		for i in range(len(self.obj.centers_embr_px)):
			if int(self.right_list.currentItem().data(0,0).toString())==i:
				self.curr_frame=i
				self.curr_right_ind=self.right_list.indexFromItem(self.right_list.currentItem())
				self.curr_right_ind=self.curr_right_ind.row()
			
				break
						
	def done_pressed(self):
		
		self.obj.ignored.sort()
		
		if 'embryo' == self.obj.__class__.__name__:
			self.obj.tvec_ignored=delete(self.obj.tvec_data,self.obj.ignored)	
		
		if 'embryo' == self.obj.__class__.__name__:
			if shape(self.obj.ext_av_data_d)[0]>0:
				self.obj.ext_av_data_ign=delete(self.obj.ext_av_data_d,self.obj.ignored)
				self.obj.int_av_data_ign=delete(self.obj.int_av_data_d,self.obj.ignored)
				self.obj.slice_av_data_ign=delete(self.obj.slice_av_data_d,self.obj.ignored)
		
		if 'bkgd' in self.obj.__class__.__name__:
			if shape(self.obj.bkgd_ext_vec)[0]>0:
				self.obj.bkgd_ext_vec_ign=delete(self.obj.bkgd_ext_vec,self.obj.ignored)
				self.obj.bkgd_int_vec_ign=delete(self.obj.bkgd_int_vec,self.obj.ignored)
				self.obj.bkgd_slice_vec_ign=delete(self.obj.bkgd_slice_vec,self.obj.ignored)
				
				self.obj.bkgd_ext
				
		self.done(1)

#===================================================================================================================================
#Dialog for selecting properties of an object that need to be handled
#===================================================================================================================================

class select_obj_props(QtGui.QDialog):
	
	def __init__(self,obj,parent):
		super(select_obj_props,self).__init__(parent)
		
		self.obj=obj
		
		self.btn_add=QtGui.QToolButton()
		self.btn_add.connect(self.btn_add, QtCore.SIGNAL('clicked()'), self.add_item)
		self.btn_add.setArrowType(QtCore.Qt.RightArrow)
		
		self.btn_remove=QtGui.QToolButton()
		self.btn_remove.connect(self.btn_remove, QtCore.SIGNAL('clicked()'), self.remove_item)
		self.btn_remove.setArrowType(QtCore.Qt.LeftArrow)
		
		self.btn_up=QtGui.QToolButton()
		self.btn_up.connect(self.btn_up, QtCore.SIGNAL('clicked()'), self.up_item)
		self.btn_up.setArrowType(QtCore.Qt.UpArrow)
		
		self.btn_down=QtGui.QToolButton()
		self.btn_down.connect(self.btn_down, QtCore.SIGNAL('clicked()'), self.down_item)
		self.btn_down.setArrowType(QtCore.Qt.DownArrow)
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		self.left_list=QtGui.QTreeWidget()
		self.left_list.setHeaderLabels(["Available Properties"])
		self.left_list.setColumnWidth(0,200)
		self.left_list.setColumnWidth(1,75)
		self.left_list.itemDoubleClicked.connect(self.add_item)
			
		self.right_list=QtGui.QTreeWidget()
		self.right_list.setHeaderLabels(["Selected Properties"])
		self.right_list.setColumnWidth(0,200)
		self.right_list.setColumnWidth(1,75)
		self.right_list.itemDoubleClicked.connect(self.remove_item)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.btn_up)
		self.vbox.addWidget(self.btn_add)
		self.vbox.addWidget(self.btn_remove)
		self.vbox.addWidget(self.btn_down)
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.left_list)
		self.hbox.addLayout(self.vbox)
		self.hbox.addWidget(self.right_list)
		
		self.vbox2 = QtGui.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btn_done)
		
		self.init_left_list()
		self.init_right_list()
		
		self.resize(600,500)
		self.setLayout(self.vbox2)
		self.show()	
	
	def init_left_list(self):
		
		#Check if object has selected props
		if hasattr(self.obj,'selected_props'):
			pass
		else:
			self.obj.selected_props=[]
		
		#Put properties in left list
		for item in vars(self.obj):
			
			#Don't put long arrays and stuff
			if isinstance(vars(self.obj)[str(item)],(int,float,str,type(None))):
				
				#If not already in selected props, put it in left list
				if item not in self.obj.selected_props:
					QtGui.QTreeWidgetItem(self.left_list,[item])
				
	def init_right_list(self):
		
		#Check if object has selected props
		if hasattr(self.obj,'selected_props'):
			pass
		else:
			self.obj.selected_props=[]
			
		#Put properties in selected props in right list
		for item in self.obj.selected_props:
			
			#If already in selected props, put it in right list
			if item in self.obj.selected_props:
				QtGui.QTreeWidgetItem(self.right_list,[item])
				
	def add_item(self):

		#Determine selected item
		self.curr_left_item=self.left_list.currentItem()
		self.curr_item=str(self.left_list.currentItem().data(0,0).toString())
		
		#Insert new node in right list
		new_node=QtGui.QTreeWidgetItem(self.right_list,[self.curr_item])
		
		#Remove node in left list
		self.curr_left_ind=self.left_list.indexFromItem(self.left_list.currentItem()).row()
		self.left_list.takeTopLevelItem(self.curr_left_ind)
		
		self.obj.selected_props.append(self.curr_item)
	
		
		
	def remove_item(self):
	
		#Determine selected item
		self.curr_right_item=self.right_list.currentItem()
		self.curr_item=str(self.right_list.currentItem().data(0,0).toString())
		
		#Insert new node in left list
		new_node=QtGui.QTreeWidgetItem(self.left_list,[self.curr_item])
		
		#Remove node in right list
		self.curr_right_ind=self.right_list.indexFromItem(self.right_list.currentItem()).row()
		self.right_list.takeTopLevelItem(self.curr_right_ind)
		
		self.obj.selected_props.remove(self.curr_item)	
		
		
	def up_item(self):
		
		#Determine selected item
		self.curr_right_item=self.right_list.currentItem()
		self.curr_item=str(self.right_list.currentItem().data(0,0).toString())
		
		#Determine index in sel prop list
		ind=self.obj.selected_props.index(self.curr_item)
		
		#Swap in list
		if ind>0:
			self.obj.selected_props[ind-1], self.obj.selected_props[ind] = self.obj.selected_props[ind], self.obj.selected_props[ind-1]
		
		#Clear list and recreate it
		self.right_list.clear()
		self.init_right_list()
		self.curr_right_item=self.right_list.topLevelItem(ind-1)
		self.right_list.setCurrentItem(self.curr_right_item)
		
	def down_item(self):
		
		#Determine selected item
		self.curr_right_item=self.right_list.currentItem()
		self.curr_item=str(self.right_list.currentItem().data(0,0).toString())
		
		#Determine index in sel prop list
		ind=self.obj.selected_props.index(self.curr_item)
		
		#Swap in list
		if ind<len(self.obj.selected_props)-1:
			self.obj.selected_props[ind+1], self.obj.selected_props[ind] = self.obj.selected_props[ind], self.obj.selected_props[ind+1]
		
		#Clear list and recreate it
		self.right_list.clear()
		self.init_right_list()
		self.curr_right_item=self.right_list.topLevelItem(ind+1)
		self.right_list.setCurrentItem(self.curr_right_item)
	
	def done_pressed(self):
		
		self.done(1)		
		
#===================================================================================================================================
#Dialog for selecting and modifying masking threshholds
#===================================================================================================================================

class select_threshhold(QtGui.QDialog):
	
	def __init__(self,emb,parent):
		super(select_threshhold,self).__init__(parent)
			
		self.emb=emb
			
		self.dpi = 100
		self.setMinimumSize(1000,500) 
		self.resize(1300,500)
		
		#Some variables
		self.xcoords2=[]
		self.ycoords2=[]
		self.pts_in_ax2=[]
		
		self.curr_img_ind=0
		
		#Making sure that embryo has all the new properties it needs
		if not hasattr(self.emb,"thresh_meth"):
			self.emb.thresh_meth="Otsu"
			self.emb.thresh=0
			
		self.curr_meth=self.emb.thresh_meth
		self.curr_thresh=self.emb.thresh
		
		
		self.methods=["Otsu","Adaptive","Manual"]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_next=QtGui.QPushButton('Next')
		self.btn_prev=QtGui.QPushButton('Prev')
		
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		self.btn_next.connect(self.btn_next, QtCore.SIGNAL('clicked()'), self.next_img)
		self.btn_prev.connect(self.btn_prev, QtCore.SIGNAL('clicked()'), self.prev_img)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_thresh = QtGui.QLabel("Threshhold = ", self)
		self.lbl_meth = QtGui.QLabel("Method:", self)
		self.lbl_max_name = QtGui.QLabel("Maximal Intensity = ", self)
		self.lbl_max = QtGui.QLabel(str(0), self)
			
		#-------------------------------------------------------------------------------------------------------------------
		#Dropdowns
		#-------------------------------------------------------------------------------------------------------------------
		
		self.combo_meth = QtGui.QComboBox(self)
		self.init_meth()
		self.update_meth(self.emb.thresh_meth)
		
		self.combo_meth.activated[str].connect(self.sel_meth)   
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_thresh = QtGui.QLineEdit(str(self.curr_thresh))
		
		self.double_valid=QtGui.QDoubleValidator()
		self.qle_thresh.setValidator(self.double_valid)
		
		self.qle_thresh.editingFinished.connect(self.set_thresh)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Plot frame
		#-------------------------------------------------------------------------------------------------------------------
		
		self.plot_frame = QtGui.QWidget()
		self.plot_frame.setMaximumWidth(1)
		
		self.hist_frame = QtGui.QWidget()
		self.hist_frame.setMaximumWidth(1)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Slider
		#-------------------------------------------------------------------------------------------------------------------
		
		self.slider=QtGui.QSlider(parent=self)
		self.slider.setOrientation(Qt.Horizontal)
		
		self.slider.setRange(0,1)
		self.slider.setSingleStep(1)
		self.connect(self.slider, QtCore.SIGNAL('sliderReleased()'), self.slider_moved)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_fill = QtGui.QCheckBox('Fill NaN', self)
		
		self.cb_fill.setCheckState(0)	
		
		self.connect(self.cb_fill, QtCore.SIGNAL('stateChanged(int)'), self.check_fill)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout Left hand side
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		vbox2 = QtGui.QVBoxLayout()
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_thresh,3,1)
		grid.addWidget(self.qle_thresh,3,2)
		grid.addWidget(self.lbl_meth,2,1)
		grid.addWidget(self.combo_meth,2,2)
		grid.setColumnMinimumWidth(2,200) 
		vbox2.addLayout(grid)
		vbox2.addWidget(self.slider)
		
		grid = QtGui.QGridLayout()
		grid.addWidget(self.lbl_max_name,1,1)
		grid.addWidget(self.lbl_max,1,2)
		grid.addWidget(self.cb_fill,2,1)
		grid.addWidget(self.btn_done,4,1)
		vbox2.addLayout(grid)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Image navigation list
		#-------------------------------------------------------------------------------------------------------------------
		
		self.prop_list=QtGui.QTreeWidget()
		self.prop_list.setHeaderLabels(["Img"])
	
		self.prop_list.setColumnWidth(0,40)
		self.prop_list.itemClicked.connect(self.prop_list_click) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Create Canvas
		#-------------------------------------------------------------------------------------------------------------------
		
		self.fig_img,self.canvas_img,self.ax_img=self.create_canvas(500,500,self.plot_frame)
		self.fig_hist,self.canvas_hist,self.ax_hist=self.create_canvas(500,500,self.hist_frame)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.hbox2 = QtGui.QHBoxLayout()
		self.hbox2.addWidget(self.btn_prev)
		self.hbox2.addWidget(self.btn_next)
		
		self.hbox3 = QtGui.QHBoxLayout()
		self.hbox3.addWidget(self.canvas_img)
		self.hbox3.addWidget(self.canvas_hist)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addLayout(self.hbox3)
		self.vbox.addLayout(self.hbox2)
		
		#Add everything to Horizontal Box
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(vbox2)
		self.hbox.addWidget(self.prop_list)
		self.hbox.addLayout(self.vbox)
		self.setLayout(self.hbox)    
		
		#Check if there is already a first image
		self.update_prop_list()
		self.show_img(self.curr_img_ind)	
		#self.draw_patches(0)
		
		self.setWindowTitle('Threshhold Editor')    
		self.show()
		
	def prop_list_click(self):
		
		self.curr_img_ind=self.prop_list.currentIndex().row()
		#self.draw_patches(self.curr_img_ind)
		self.show_img(self.curr_img_ind)
	
	def update_prop_list(self):
		
		self.file_list = pyfdap_misc.get_sorted_folder_list(self.emb.fn_maskfolder,self.emb.data_ft)
		self.nodes=[]
		
		for i in range(shape(self.file_list)[0]):
			
			self.nodes.append(QtGui.QTreeWidgetItem(self.prop_list,[str(i)]))
		
		if shape(self.file_list)[0]>0:
			self.prop_list.setCurrentItem(self.nodes[0])
		
		else:
			pass
	
	def slider_moved(self):
		value=self.slider.value()
		self.qle_thresh.setText(str(value))
		self.set_thresh()
		
	
	def next_img(self):
		if self.curr_img_ind==shape(self.file_list)[0]-1:
			self.prop_list.setCurrentItem(self.nodes[0])
			#self.draw_patches(0)
			self.show_img(0)
		else:	
			self.prop_list.setCurrentItem(self.nodes[self.curr_img_ind+1])
			#self.draw_patches(self.curr_img_ind+1)
			self.show_img(self.curr_img_ind+1)
					
	def prev_img(self):
		if self.curr_img_ind==0:
			self.prop_list.setCurrentItem(self.nodes[shape(self.file_list)[0]-1])
			#self.draw_patches(shape(self.file_list)[0]-1)
			self.show_img(shape(self.file_list)[0]-1)	
		else:	
			self.prop_list.setCurrentItem(self.nodes[self.curr_img_ind-1])
			#self.draw_patches(self.curr_img_ind-1)
			self.show_img(self.curr_img_ind-1)
			
	def create_canvas(self,h,v,parent):
			
		h=h/self.dpi
		v=v/self.dpi
		fig = Figure( dpi=self.dpi)
		fig.set_size_inches(h,v,forward=True)
		canvas = FigureCanvas(fig)
		canvas.setParent(parent)
		
		ax = fig.add_subplot(111)
		ax.set_xlim([0, self.emb.data_res_px])
		ax.set_ylim([0, self.emb.data_res_px])
		
		canvas.draw()
		return fig,canvas,ax 
	
	def show_img(self,ind):
		
		self.file_list=pyfdap_misc.get_sorted_folder_list(self.emb.fn_maskfolder,self.emb.data_ft)
		
		#Check if there is a first image
		if shape(self.file_list)[0]>ind-1 and shape(self.file_list)[0]>0:
		
			#Grab first picture
			fn=self.emb.fn_maskfolder+self.file_list[ind]
			
			#Load img
			data_img = mpimg.imread(fn).astype(self.emb.data_enc)
			data_vals=data_img.real
			data_vals=data_vals.astype('float')
			
			#Flip to be coherent with Fiji
			data_vals=flipud(data_vals)
			
			#Updating max intensity of slider
			self.slider.setRange(0,data_vals.max())
			self.lbl_max.setText(str(data_vals.max()))
			
			#Perform threshholding
			self.curr_img,self.curr_thresh=self.perform_thresh(self.curr_meth,data_vals,thresh=self.curr_thresh)
			self.update_qles()
			
			#Plot img
			self.ax_img.imshow(self.curr_img)
			self.canvas_img.draw()
			
			#Plot histogram
			self.draw_hist(data_vals)
			
		else:
			#Clear if there is no first image
			self.ax_img.cla()	
	
	def draw_hist(self,img):
		
		bins=linspace(nanmin(img),nanmax(img),256+1)
		data,bin_edges=histogram(img,bins)
		
		bins=linspace(nanmin(img),nanmax(img),256)
		self.ax_hist.cla()
		self.ax_hist.bar(bins,data,color='b')
		y=arange(data.max())
		x=ones(shape(y))*self.curr_thresh
		
		self.ax_hist.plot(x,y,'r')
		self.ax_hist.autoscale()
		self.canvas_hist.draw()
		
	def sel_meth(self,text):
		self.curr_meth=str(text)
		self.show_img(self.curr_img_ind)	
	
	def set_thresh(self):
		text=self.qle_thresh.text()
		self.curr_thresh=float(str(text))
		self.update_meth("Manual")
		self.curr_meth="Manual"
		self.show_img(self.curr_img_ind)
		
	def perform_thresh(self,meth,img,thresh=0):
		
		if meth=="Otsu":
			thresh,mask=pyfdap_img.otsu_imagej(img,1,0,0)
		elif meth=="Adaptive":
			mask = pyfdap_img.adaptive_thresh(img,5)
		elif meth=="Manual":
			mask = pyfdap_img.fixed_thresh(img,thresh)
			
		return mask,thresh	
	
	def update_qles(self):
		self.qle_thresh.setText(str(self.curr_thresh))
		self.slider.setValue(self.curr_thresh)
		
	def update_meth(self,meth):
		
		meth_ind=self.methods.index(meth)
		self.combo_meth.setCurrentIndex(meth_ind)	
		
	def init_meth(self):
		for m in self.methods:
			self.combo_meth.addItem(m)
	
	def check_fill(self,value):
		
		if value==2:
			self.emb.fill_mask=nan	
		else:	
			self.emb.fill_mask=0
			
		self.show_img(self.curr_img_ind)
		
	def done_pressed(self):
		
		self.emb.thresh_meth=self.curr_meth
		self.emb.tresh=self.curr_thresh
		
		
#===================================================================================================================================
#Dialog for importing embryos via csv sheet
#===================================================================================================================================
		

class csv_read_dialog(QtGui.QDialog):
	def __init__(self,mol,parent):
		super(csv_read_dialog,self).__init__(parent)
		
		self.mol=mol
		self.default_ident=["time","ext","int","slice"]
		self.parent=parent
		self.default_delim=','
		self.default_idx_pre=[0,1]
		self.default_insamecol=True
		self.default_noise=1.
			
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_set_csvfile=QtGui.QPushButton('Change')
		
		#Button Actions
		self.btn_set_csvfile.connect(self.btn_set_csvfile, QtCore.SIGNAL('clicked()'), self.sel_csvfile)
		
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_ident = QtGui.QLabel("Identifiers:", self)
		self.lbl_csvfile = QtGui.QLabel("CSV File:", self)
		self.lbl_csvpath = QtGui.QLabel("", self)
		self.lbl_delim = QtGui.QLabel("Delimiter", self)
		self.lbl_idx_pre = QtGui.QLabel("Indices with Pre-values", self)
		self.lbl_insamecol = QtGui.QLabel("Pre-values in same column?", self)
		self.lbl_noise = QtGui.QLabel("Noise:", self)
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_ident = QtGui.QLineEdit(str(self.default_ident))
		self.qle_delim = QtGui.QLineEdit(str(self.default_delim))
		self.qle_idx_pre = QtGui.QLineEdit(str(self.default_idx_pre))
		self.qle_noise = QtGui.QLineEdit(str(self.default_noise))
		
		
		self.qle_ident.editingFinished.connect(self.set_ident)
		self.qle_delim.editingFinished.connect(self.set_delim)
		self.qle_idx_pre.editingFinished.connect(self.set_idx_pre)
		self.qle_noise.editingFinished.connect(self.set_noise)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_insamecol = QtGui.QCheckBox('', self)
		self.cb_insamecol.setCheckState(int(self.default_insamecol))	
		self.connect(self.cb_insamecol, QtCore.SIGNAL('stateChanged(int)'), self.check_insamecol)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_ident,1,1)
		grid.addWidget(self.lbl_csvfile,2,1)
		grid.addWidget(self.lbl_delim,3,1)
		grid.addWidget(self.lbl_insamecol,4,1)
		grid.addWidget(self.lbl_idx_pre,5,1)
		grid.addWidget(self.lbl_noise,6,1)
		
		grid.addWidget(self.lbl_csvpath,2,2)
		grid.addWidget(self.btn_set_csvfile,2,3)
		grid.addWidget(self.qle_ident,1,2)
		grid.addWidget(self.qle_delim,3,2)
		grid.addWidget(self.cb_insamecol,4,2)
		grid.addWidget(self.qle_idx_pre,5,2)
		grid.addWidget(self.qle_noise,6,2)
		
		grid.addWidget(self.btn_done,8,1)
		
		grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.setLayout(grid)    
		
		self.setWindowTitle('Read csv file')    
		self.show()
	
	def sel_csvfile(self):
		
		fn_load=str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.parent.lastopen,"*.csv",))
		if fn_load=='':
			return
		
		self.fn_csv=fn_load
		self.parent.lastopen=os.path.dirname(fn_load)
	
		if len(self.fn_csv)>50:
			self.lbl_csvpath.setText("..."+self.fn_csv[-50:])
		else:
			self.lbl_csvpath.setText(self.fn_csv)
	
	def set_delim(self):
		self.default_delim=str(self.qle_delim.text())
	
	def set_ident(self):
		self.default_ident=self.str2list(str(self.qle_ident.text()),str)
		
	def set_idx_pre(self):
		
		s=str(self.qle_idx_pre.text())
		self.default_idx_pre=self.str2list(s,int)
		
	def str2list(self,s,dtype):
			
		s=s.replace("[","").replace("]","")
		s=s.split(",")
		l=[]
		for x in s:
			l.append(dtype(x))
		
		return l
		
	def set_noise(self):
		self.default_noise=float(str(self.qle_noise.text()))
		
	def check_insamecol(self,value):
		self.default_insamecol=bool(int(value))
	
	def done_pressed(self):
		
		if os.path.exists(self.fn_csv):
			self.embryos=pyfdap_misc.gen_embryos_from_csv(self.fn_csv,ident=self.default_ident,delim=self.default_delim,debug=True,mol=self.mol,
						 inSameCol=self.default_insamecol,idxPre=self.default_idx_pre,noiseVal=self.default_noise)
			self.done(1)
			return self.mol
			
		else:
			ret=QtGui.QMessageBox.warning(self,self.tr("Warning"),self.tr("CSV file does not exist. Exit?"),QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
					
			if ret == QtGui.QMessageBox.Yes:
				self.done(1)
				return self.mol
			else:	
				return
	
	def get_embryos(self):
		
		return self.embryos


#===================================================================================================================================
#Dialog for importing bkgds via csv sheet
#===================================================================================================================================
		

class csv_read_dialog_bkgds(QtGui.QDialog):
	def __init__(self,mol,parent):
		super(csv_read_dialog_bkgds,self).__init__(parent)
		
		self.mol=mol
		self.default_ident=["time","ext","int","slice"]
		self.parent=parent
		self.default_delim=','
		self.default_idx_pre=[0,1]
		self.default_insamecol=True
			
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btn_done=QtGui.QPushButton('Done')
		self.btn_set_csvfile=QtGui.QPushButton('Change')
		
		#Button Actions
		self.btn_set_csvfile.connect(self.btn_set_csvfile, QtCore.SIGNAL('clicked()'), self.sel_csvfile)
		
		self.btn_done.connect(self.btn_done, QtCore.SIGNAL('clicked()'), self.done_pressed)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Labels
		#-------------------------------------------------------------------------------------------------------------------
		
		self.lbl_ident = QtGui.QLabel("Identifiers:", self)
		self.lbl_csvfile = QtGui.QLabel("CSV File:", self)
		self.lbl_csvpath = QtGui.QLabel("", self)
		self.lbl_delim = QtGui.QLabel("Delimiter", self)
		self.lbl_idx_pre = QtGui.QLabel("Indices with Pre-values", self)
		self.lbl_insamecol = QtGui.QLabel("Pre-values in same column?", self)
		
		#-------------------------------------------------------------------------------------------------------------------
		#LineEdits
		#-------------------------------------------------------------------------------------------------------------------
		
		self.qle_ident = QtGui.QLineEdit(str(self.default_ident))
		self.qle_delim = QtGui.QLineEdit(str(self.default_delim))
		self.qle_idx_pre = QtGui.QLineEdit(str(self.default_idx_pre))
		
		self.qle_ident.editingFinished.connect(self.set_ident)
		self.qle_delim.editingFinished.connect(self.set_delim)
		self.qle_idx_pre.editingFinished.connect(self.set_idx_pre)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Checkboxes
		#-------------------------------------------------------------------------------------------------------------------
		
		self.cb_insamecol = QtGui.QCheckBox('', self)
		self.cb_insamecol.setCheckState(int(self.default_insamecol))	
		self.connect(self.cb_insamecol, QtCore.SIGNAL('stateChanged(int)'), self.check_insamecol)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		#Grid
		grid = QtGui.QGridLayout()
		
		grid.addWidget(self.lbl_ident,1,1)
		grid.addWidget(self.lbl_csvfile,2,1)
		grid.addWidget(self.lbl_delim,3,1)
		grid.addWidget(self.lbl_insamecol,4,1)
		grid.addWidget(self.lbl_idx_pre,5,1)
		
		grid.addWidget(self.lbl_csvpath,2,2)
		grid.addWidget(self.btn_set_csvfile,2,3)
		grid.addWidget(self.qle_ident,1,2)
		grid.addWidget(self.qle_delim,3,2)
		grid.addWidget(self.cb_insamecol,4,2)
		grid.addWidget(self.qle_idx_pre,5,2)
		
		grid.addWidget(self.btn_done,7,1)
		
		grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.setLayout(grid)    
		
		self.setWindowTitle('Read csv file')    
		self.show()
	
	def sel_csvfile(self):
		
		fn_load=str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.parent.lastopen,"*.csv",))
		if fn_load=='':
			return
		
		self.fn_csv=fn_load
	
		self.parent.lastopen=os.path.dirname(fn_load)
		if len(self.fn_csv)>50:
			self.lbl_csvpath.setText("..."+self.fn_csv[-50:])
		else:
			self.lbl_csvpath.setText(self.fn_csv)
	
	def set_delim(self):
		self.default_delim=str(self.qle_delim.text())
	
	def set_ident(self):
		self.default_ident=str(self.qle_ident.text())
		
	def set_idx_pre(self):
		self.default_idx_pre=str(self.qle_idx_pre.text())
		
	def check_insamecol(self,value):
		self.default_insamecol=bool(int(value))
	
	def done_pressed(self):
		
		if os.path.exists(self.fn_csv):
			self.bkgds=pyfdap_misc.gen_bkgds_from_csv(self.fn_csv,self.mol,ident=self.default_ident,delim=self.default_delim,debug=True,
						 inSameCol=self.default_insamecol,idxPre=self.default_idx_pre)
			self.done(1)
			return self.mol
			
		else:
			ret=QtGui.QMessageBox.warning(self,self.tr("Warning"),self.tr("CSV file does not exist. Exit?"),QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
					
			if ret == QtGui.QMessageBox.Yes:
				self.done(1)
				return self.mol
			else:	
				return
	
	def get_bkgds(self):
		
		return self.bkgds	
	
class listSelectorDialog(QtGui.QDialog):
	
	def __init__(self,parent,List,leftTitle="",rightTitle="",itemsRight=[],maxRight=2):
		super(listSelectorDialog,self).__init__(parent)
		#print type(self), type(parent)
		
		#QtGui.QDialog.__init__()
		
		
		self.itemsRight=itemsRight
		self.itemsLeft=list(List)
		self.maxRight=maxRight
		
		
		self.List=List
		
		#Buttons
		self.btnAdd=QtGui.QToolButton()
		self.btnAdd.connect(self.btnAdd, QtCore.SIGNAL('clicked()'), self.addItem)
		self.btnAdd.setArrowType(QtCore.Qt.RightArrow)
		
		self.btnRemove=QtGui.QToolButton()
		self.btnRemove.connect(self.btnRemove, QtCore.SIGNAL('clicked()'), self.removeItem)
		self.btnRemove.setArrowType(QtCore.Qt.LeftArrow)
		
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Left QtreeWidgetItem
		self.leftList=QtGui.QTreeWidget()
		self.leftList.setHeaderLabels([leftTitle])
		self.leftList.setColumnWidth(0,200)
		self.leftList.setColumnWidth(1,75)
		self.leftList.itemDoubleClicked.connect(self.addItem)
		
		#right QtreeWidgetItem
		self.rightList=QtGui.QTreeWidget()
		self.rightList.setHeaderLabels([rightTitle])
		self.rightList.setColumnWidth(0,200)
		self.rightList.setColumnWidth(1,75)
		self.rightList.itemDoubleClicked.connect(self.removeItem)
		
		#Layout
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.btnAdd)
		self.vbox.addWidget(self.btnRemove)
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addWidget(self.leftList)
		self.hbox.addLayout(self.vbox)
		self.hbox.addWidget(self.rightList)
		
		self.vbox2 = QtGui.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btnDone)
		
		#Init lists
		self.initLeftList()
		self.initRightList()
		
		self.resize(400,500)
		self.setLayout(self.vbox2)
		self.setWindowTitle("list Selector Dialog")
		
		self.show()
	
	def getListDifference(self):
		
		for item in self.itemsLeft:
			if item in self.itemsRight:
				self.itemsLeft.remove(item)
		
		return self.itemsLeft
				
	
	def initLeftList(self):
		
		self.getListDifference()
		
		for item in self.itemsLeft:
			QtGui.QTreeWidgetItem(self.leftList,[item])
			
	def initRightList(self):
		
		for item in self.itemsRight:
			QtGui.QTreeWidgetItem(self.rightList,[item])
				
	def addItem(self):
		
		#Check if there are already enough elements on the right hand side
		if len(self.itemsRight)>self.maxRight:
			return
		
		#Determine selected item
		self.currentItem=str(self.leftList.currentItem().data(0,0).toString())
		
		#Insert new node in right list
		newNode=QtGui.QTreeWidgetItem(self.rightList,[self.currentItem])
		
		#Remove node in left list
		self.currLeftInd=self.leftList.indexFromItem(self.leftList.currentItem()).row()
		self.leftList.takeTopLevelItem(self.currLeftInd)
		
		self.itemsRight.append(self.currentItem)
		self.itemsLeft.remove(self.currentItem)
			
	def removeItem(self):
	
		#Determine selected item
		self.currentItem=str(self.rightList.currentItem().data(0,0).toString())
		
		#Insert new node in left list
		newNode=QtGui.QTreeWidgetItem(self.leftList,[self.currentItem])
		
		#Remove node in right list
		self.currRightInd=self.rightList.indexFromItem(self.rightList.currentItem()).row()
		self.rightList.takeTopLevelItem(self.currRightInd)
		
		self.itemsRight.remove(self.currentItem)
		self.itemsLeft.append(self.currentItem)

	def getSelection(self):
		return self.itemsRight
	
	def donePressed(self):
		self.done(1)
		return

	