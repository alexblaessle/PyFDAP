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

#Module containing PyQT classes needed for PyFRAP GUI controlling plot behaviour:

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

import sys
from numpy import *
from PyQt4 import QtGui, QtCore
import pyfdap_img_module as pyfdap_img
import pyfdap_misc_module as pyfdap_misc

import matplotlib.pyplot as plt
import matplotlib

from matplotlib.figure import Figure


class modifyPlotDialog(QtGui.QDialog):
	
	def __init__(self,ax,parent):
		super(modifyPlotDialog,self).__init__(parent)
		
		#Initiate values
		self.ax=ax
		self.initUnit()
		
		#Buttons
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Label modification
		self.lbl_xlabel = QtGui.QLabel("Xlabel:", self)
		self.qle_xlabel = QtGui.QLineEdit(str(self.ax.get_xlabel()))
		self.qle_xlabel.editingFinished.connect(self.set_xlabel)
		
		self.lbl_ylabel = QtGui.QLabel("Ylabel:", self)
		self.qle_ylabel = QtGui.QLineEdit(str(self.ax.get_ylabel()))
		self.qle_ylabel.editingFinished.connect(self.set_ylabel)
		
		#xlim modification
		self.lbl_xlim = QtGui.QLabel("xlim:", self)
		self.qle_xlim_min = QtGui.QLineEdit(str(self.ax.get_xlim()[0]))
		self.qle_xlim_max = QtGui.QLineEdit(str(self.ax.get_xlim()[1]))
		self.qle_xlim_min.editingFinished.connect(self.set_xlim)
		self.qle_xlim_max.editingFinished.connect(self.set_xlim)
		
		self.lbl_ylim = QtGui.QLabel("ylim:", self)
		self.qle_ylim_min = QtGui.QLineEdit(str(self.ax.get_ylim()[0]))
		self.qle_ylim_max = QtGui.QLineEdit(str(self.ax.get_ylim()[1]))
		self.qle_ylim_min.editingFinished.connect(self.set_ylim)
		self.qle_ylim_max.editingFinished.connect(self.set_ylim)
		
		#Unit
		self.lbl_unit_xaxis = QtGui.QLabel("Unit x-axis:", self)
		
		self.combo_unit_xaxis = QtGui.QComboBox(self)
		self.combo_unit_xaxis.addItem("min")
		self.combo_unit_xaxis.addItem("sec")
		self.combo_unit_xaxis.addItem("None")
		self.initUnitCombo()
		self.combo_unit_xaxis.activated[str].connect(self.selectUnit)  
		
		#Xtick modification
		self.lbl_xticks = QtGui.QLabel("xticks", self)
		self.lbl_xticks.setAlignment(QtCore.Qt.AlignCenter)
		self.xtickList=QtGui.QTreeWidget()
		self.xtickList.setHeaderLabels(["lbl","pos"])
		self.xtickList.setColumnWidth(0,75)
		self.xtickList.setColumnWidth(1,75)
		self.xtickList.itemDoubleClicked.connect(self.editXTick)
		
		self.btnAddXTick=QtGui.QPushButton('Add xtick')
		self.btnAddXTick.connect(self.btnAddXTick, QtCore.SIGNAL('clicked()'), self.addXTick)
		
		self.btnRemoveXTick=QtGui.QPushButton('Remove xtick')
		self.btnRemoveXTick.connect(self.btnRemoveXTick, QtCore.SIGNAL('clicked()'), self.removeXTick)
		
		self.btnResetXTick=QtGui.QPushButton('Reset xtick')
		self.btnResetXTick.connect(self.btnResetXTick, QtCore.SIGNAL('clicked()'), self.resetXTick)
		
		#Ytick modification
		self.lbl_yticks = QtGui.QLabel("yticks", self)
		self.lbl_yticks.setAlignment(QtCore.Qt.AlignCenter)
		self.ytickList=QtGui.QTreeWidget()
		self.ytickList.setHeaderLabels(["lbl","pos"])
		self.ytickList.setColumnWidth(0,75)
		self.ytickList.setColumnWidth(1,75)
		self.ytickList.itemDoubleClicked.connect(self.editYTick)
		
		self.btnAddYTick=QtGui.QPushButton('Add ytick')
		self.btnAddYTick.connect(self.btnAddYTick, QtCore.SIGNAL('clicked()'), self.addYTick)
		
		self.btnRemoveYTick=QtGui.QPushButton('Remove ytick')
		self.btnRemoveYTick.connect(self.btnRemoveYTick, QtCore.SIGNAL('clicked()'), self.removeYTick)
		
		self.btnResetYTick=QtGui.QPushButton('Reset ytick')
		self.btnResetYTick.connect(self.btnResetYTick, QtCore.SIGNAL('clicked()'), self.resetYTick)
		
		#Lines modification
		self.lbl_lineList = QtGui.QLabel("Lines", self)
		self.lbl_lineList.setAlignment(QtCore.Qt.AlignCenter)
		self.linesList=QtGui.QTreeWidget()
		self.linesList.setHeaderLabels(["label","linestyle","marker"])
		self.linesList.setColumnWidth(0,75)
		self.linesList.setColumnWidth(1,75)
		self.linesList.setColumnWidth(2,75)
		self.linesList.itemDoubleClicked.connect(self.editLine)
		
		#Layout tick lists
		self.vboxXTick = QtGui.QVBoxLayout()
		self.vboxXTick.addWidget(self.lbl_xticks)
		self.vboxXTick.addWidget(self.xtickList)
		self.vboxXTick.addWidget(self.btnAddXTick)
		self.vboxXTick.addWidget(self.btnRemoveXTick)
		self.vboxXTick.addWidget(self.btnResetXTick)
		
		self.vboxYTick = QtGui.QVBoxLayout()
		self.vboxYTick.addWidget(self.lbl_yticks)
		self.vboxYTick.addWidget(self.ytickList)
		self.vboxYTick.addWidget(self.btnAddYTick)
		self.vboxYTick.addWidget(self.btnRemoveYTick)
		self.vboxYTick.addWidget(self.btnResetYTick)
		
		#Layout line list
		self.vboxLine = QtGui.QVBoxLayout()
		self.vboxLine.addWidget(self.lbl_lineList)
		self.vboxLine.addWidget(self.linesList)
		
		#Common layout for lists
		self.hboxTicks = QtGui.QHBoxLayout()
	
		self.hboxTicks.addLayout(self.vboxXTick)
		self.hboxTicks.addLayout(self.vboxYTick)
		self.hboxTicks.addLayout(self.vboxLine)
		
		#Layout for label column
		self.grid = QtGui.QGridLayout()
		
		self.grid.addWidget(self.lbl_xlabel,1,1)
		self.grid.addWidget(self.lbl_ylabel,2,1)
		self.grid.addWidget(self.lbl_xlim,3,1)
		self.grid.addWidget(self.lbl_ylim,4,1)
		self.grid.addWidget(self.lbl_unit_xaxis,5,1)
		
		self.grid.addWidget(self.qle_xlabel,1,2)
		self.grid.addWidget(self.qle_ylabel,2,2)
		self.grid.addWidget(self.qle_xlim_min,3,2)
		self.grid.addWidget(self.qle_ylim_min,4,2)
		self.grid.addWidget(self.combo_unit_xaxis,5,2)
		
		self.grid.addWidget(self.qle_xlim_max,3,3)
		self.grid.addWidget(self.qle_ylim_max,4,3)
		
		#Final layout
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(self.grid)
		self.hbox.addLayout(self.hboxTicks)
		
		self.vbox2 = QtGui.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btnDone, QtCore.Qt.AlignRight )
		
		#Init lists
		self.initXTickList()
		self.initYTickList()
		self.initLineList()
		
		self.resize(800,500)
		self.setLayout(self.vbox2)
		self.setWindowTitle("Modify Plot Dialog")
		
		self.show()
		
	def initLineList(self):
		
		"""Initiates list of lines."""
		
		self.linesList.clear()
		
		for l in self.ax.get_lines():
	
			item=QtGui.QTreeWidgetItem(self.linesList,[str(l.get_label()),str(l.get_linestyle()),str(l.get_marker())])
			
			b=QtGui.QBrush(self.color2QT(l.get_color(),l.get_alpha()))
			b2=QtGui.QBrush(self.color2QT(l.get_markerfacecolor(),1.))
			
			item.setForeground( 1 , b)
			item.setForeground( 2 , b2)

	def initXTickList(self):
		
		"""Initiates entries in xticklist."""
		
		self.xtickList.clear()
		
		for i,t in enumerate(self.ax.get_xticks()):
			QtGui.QTreeWidgetItem(self.xtickList,[str(self.ax.get_xticklabels()[i].get_text()),str(t)])
	
	def initYTickList(self):
		
		"""Initiates entries in yticklist."""
		
		self.ytickList.clear()
		
		for i,t in enumerate(self.ax.get_yticks()):
			QtGui.QTreeWidgetItem(self.ytickList,[str(self.ax.get_yticklabels()[i].get_text()),str(t)])
	
	def color2Hex(self,c):
		
		"""Converts color from MPL to hex code such that QT can use it."""
		
		colorRGB=self.color2RGB(c)
		colorHex=matplotlib.colors.rgb2hex(colorRGB) 
		
		return colorHex
	
	def set_xlabel(self):
		
		"""Sets xlabel when qle_xlabel was used."""
		
		text=str(self.qle_xlabel.text())
		self.ax.set_xlabel(text)
		
	def set_ylabel(self):
		
		"""Sets ylabel when qle_ylabel was used."""
		
		text=str(self.qle_ylabel.text())
		self.ax.set_ylabel(text)
	
	def initUnit(self):
	
		"""Initiates time unit."""
	
		if "min" in self.ax.get_xlabel():
			self.unit="min"
		elif "(s)" in self.ax.get_xlabel():
			self.unit="sec"
		elif "sec" in self.ax.get_xlabel():
			self.unit="sec"
		else:
			self.unit="None"
			
	def initUnitCombo(self):
		
		"""Initiates combo box for unit selection."""
		
		if self.unit=="sec":
			self.combo_unit_xaxis.setCurrentIndex(0)
		elif self.unit=="min":
			self.combo_unit_xaxis.setCurrentIndex(1)
		elif self.unit=="None":	
			self.combo_unit_xaxis.setCurrentIndex(2)
	
	def updateLabelQles(self):
		
		"""Updates Qles for labels."""
		
		self.qle_xlabel.setText(self.ax.get_xlabel())
		self.qle_ylabel.setText(self.ax.get_ylabel())
		
	def selectUnit(self,text):
		
		"""Updates xticklabels depending on unit."""
		
		oldUnit=str(self.unit)
		self.unit=str(text)
		
		if oldUnit==self.unit:
			return
		
		lbls=[]
		if oldUnit=="sec" and self.unit=="min":	
			for t in self.ax.get_xticklabels():
				lbls.append(round(float(t.get_text())/60.,2))
		if oldUnit=="min" and self.unit=="sec":		
			for t in self.ax.get_xticklabels():
				lbls.append(float(t.get_text())*60.)
		if len(lbls)>0:
			self.ax.set_xticklabels(lbls)
		
		self.ax.set_xlabel("Time ("+self.unit+")")
		self.updateLabelQles()
		
		self.redrawAxes()
		
	def redrawAxes(self):
		
		"""Redraws axes."""
		
		self.ax.get_figure().canvas.draw()
	
	def color2RGB(self,c,alpha):
		
		"""Converts color from MPL to RGB such that QT can use it."""
		
		colorRGB=list((array(matplotlib.colors.colorConverter.to_rgb(c),dtype=float)*255).astype(int))
		
		if not isinstance(alpha,float):
			alpha=255
		else:
			alpha=int(255.*alpha)
		
		colorRGB.append(alpha)
		
		return tuple(colorRGB)
	
	def color2QT(self,c,alpha):
		
		"""Converts color from MPL to QColor such that QT can use it."""
		
		C = self.color2RGB(c,alpha)
		C = QtGui.QColor(C[0],C[1],C[2],C[3])
		
		return C
	
	def getCurrentIndex(self,treeWidget):
		
		"""Gets current index of a treeWidget."""
		
		return int(treeWidget.indexFromItem(treeWidget.currentItem()).row())
	
	
	def getCurrentLine(self):
		
		"""Returns the current selected line."""
		
		idx=self.getCurrentIndex(self.linesList)
		
		return self.ax.get_lines()[idx]
	
	def editLine(self):
		
		"""Opens line editor."""
		
		line=self.getCurrentLine()
		
		lineEditor = lineDialog(line,self)
		if lineEditor.exec_():
			line = lineEditor.getLine()
		
		self.redrawAxes()
		self.initLineList()
	
	def getTickLabels(self,ticklabels):
		
		"""Returns the tick labels as list."""
		
		lbls=[]
		
		for t in ticklabels:
			lbls.append(t.get_text())
			
		return lbls
	
	def editXTick(self):
		
		"""Lets user specific xtick."""
		
		# Grab tick
		idx=self.getCurrentIndex(self.xtickList)
		pos=self.ax.get_xticks()[idx]
		lbl=self.ax.get_xticklabels()[idx].get_text()
		
		# Call dialog
		lbl,pos=self.editTick(lbl,pos)
		
		# Generate list and insert for xticklabels
		lbls=self.getTickLabels(self.ax.get_xticklabels())
		lbls[idx]=lbl
		self.ax.set_xticklabels(lbls)
		
		# Generate list and insert for xticks
		ticks=self.ax.get_xticks()
		ticks[idx]=pos
		self.ax.set_xticks(ticks)
		
		self.redrawAxes()
		self.initXTickList()
	
	def editYTick(self):
		
		"""Lets user specific ytick."""
		
		idx=self.getCurrentIndex(self.ytickList)
		pos=self.ax.get_yticks()[idx]
		lbl=self.ax.get_yticklabels()[idx].get_text()
		lbl,pos=self.editTick(lbl,pos)
		
		# Generate list and insert for xticklabels
		lbls=self.getTickLabels(self.ax.get_yticklabels())
		lbls[idx]=lbl
		self.ax.set_yticklabels(lbls)
		
		# Generate list and insert for xticks
		ticks=self.ax.get_yticks()
		ticks[idx]=pos
		self.ax.set_yticks(ticks)
	
		self.redrawAxes()
		self.initYTickList()
	
	def editTick(self,lbl,pos):
		
		"""Opens tick editor."""
		
		tickEditor = tickDialog(lbl,pos,self)
		if tickEditor.exec_():
			lbl,pos = tickEditor.getParms()
		
		return lbl,pos
	
	def addXTick(self):
		
		lbl,pos=self.editTick("",0.)
		
		# Generate list and insert for xticks
		ticks=list(self.ax.get_xticks())
		ticks.append(pos)
		ticks.sort()
		idx=ticks.index(pos)
		self.ax.set_xticks(ticks)
		
		# Generate list and insert for xticklabels
		lbls=self.getTickLabels(self.ax.get_xticklabels())
		lbls.insert(idx,lbl)
		self.ax.set_xticklabels(lbls)
		
		self.redrawAxes()
		self.initXTickList()
	
	def removeXTick(self):
		
		#Get index
		idx=self.getCurrentIndex(self.xtickList)
		
		#Get lbls and ticks
		lbls=self.getTickLabels(self.ax.get_xticklabels())
		ticks=list(self.ax.get_xticks())
		
		#Pop
		lbls.pop(idx)
		ticks.pop(idx)
		
		#Update
		self.ax.set_xticklabels(lbls)
		self.ax.set_xticks(ticks)
		
		self.redrawAxes()
		self.initXTickList()
		
	def addYTick(self):
		
		lbl,pos=self.editTick("",0.)
		
		# Generate list and insert for yticks
		ticks=list(self.ax.get_yticks())
		ticks.append(pos)
		ticks.sort()
		idx=ticks.index(pos)
		self.ax.set_yticks(ticks)
		
		# Generate list and insert for yticklabels
		lbls=self.getTickLabels(self.ax.get_yticklabels())
		lbls.insert(idx,lbl)
		self.ax.set_yticklabels(lbls)
		
		self.redrawAxes()
		self.initYTickList()
	
	def removeYTick(self):
		
		#Get index
		idx=self.getCurrentIndex(self.ytickList)
		
		#Get lbls and ticks
		lbls=self.getTickLabels(self.ax.get_yticklabels())
		ticks=list(self.ax.get_yticks())
		
		#Pop
		lbls.pop(idx)
		ticks.pop(idx)
		
		#Update
		self.ax.set_yticklabels(lbls)
		self.ax.set_yticks(ticks)
		
		self.redrawAxes()
		self.initYTickList()	
	
	def getXVec(self):
		
		"""Returns x-vector"""
		
		vec=linspace(min(self.ax.get_xticks()),max(self.ax.get_xticks()),len(self.ax.get_xticks()))
		
		return vec
			
	def getYVec(self):
		
		"""Returns x-vector"""
		
		vec=linspace(min(self.ax.get_yticks()),max(self.ax.get_yticks()),len(self.ax.get_yticks()))
		
		return vec
			
	def resetYTick(self):
		
		"""Edits y vector."""
		
		vec=self.editVec(self.getYVec())
		self.ax.set_yticks(vec)
		self.ax.set_yticklabels(vec)
		self.redrawAxes()
		self.initYTickList()
		self.unit=="sec"
		
	def resetXTick(self):
		
		"""Edits x vector."""
		
		vec=self.editVec(self.getXVec())
		self.ax.set_xticks(vec)
		self.ax.set_xticklabels(vec)
		self.redrawAxes()
		self.initXTickList()
		self.unit=="sec"
		
	def editVec(self,vec):
		
		"""Edits vector."""
		
		vecEditor = resetTickVecDialog(vec,self)
		if vecEditor.exec_():
			vec = vecEditor.getVec()
		
		return vec
	
	def set_xlim(self):
		
		self.ax.set_xlim([float(self.qle_xlim_min.text()),float(self.qle_xlim_max.text())])
		self.redrawAxes()
		
	def set_ylim(self):
		
		self.ax.set_ylim([float(self.qle_ylim_min.text()),float(self.qle_ylim_max.text())])
		self.redrawAxes()
		
	def donePressed(self):
		self.redrawAxes()
		self.done(1)
	
class lineDialog(QtGui.QDialog):
	
	def __init__(self,line,parent):
		super(lineDialog,self).__init__(parent)
		
		self.line=line
		
		self.linestyles = ['_', '-', '--', ':']
		self.initMarkers()
		
		#Linestyle
		self.lbl_linestyle = QtGui.QLabel("Linestyle:", self)
		self.initLinestyleCombo()
		
		#Marker
		self.lbl_marker = QtGui.QLabel("Marker:", self)
		self.initMarkerCombo()
		
		#Color
		self.lbl_color = QtGui.QLabel("Color:", self)
		self.initColorButton()
		
		#Marker Color
		self.lbl_markercolor = QtGui.QLabel("Marker Color:", self)
		self.initMarkerColorButton()
		
		#Linewidth
		self.lbl_linewidth = QtGui.QLabel("Linewidth:", self)
		self.qle_linewidth = QtGui.QLineEdit(str(self.line.get_linewidth()))
		self.qle_linewidth.editingFinished.connect(self.setLinewidth)
		
		#Alpha
		self.lbl_alpha = QtGui.QLabel("Alpha:", self)
		self.qle_alpha = QtGui.QLineEdit(str(self.line.get_alpha()))
		self.qle_alpha.editingFinished.connect(self.setAlpha)
		
		#Done button
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Layout
		self.grid = QtGui.QGridLayout()
		
		self.grid.addWidget(self.lbl_linestyle,1,1)
		self.grid.addWidget(self.lbl_color,2,1)
		self.grid.addWidget(self.lbl_marker,3,1)
		self.grid.addWidget(self.lbl_markercolor,4,1)
		self.grid.addWidget(self.lbl_linewidth,5,1)
		self.grid.addWidget(self.lbl_alpha,6,1)
		
		self.grid.addWidget(self.combo_linestyle,1,2)
		self.grid.addWidget(self.btnColor,2,2)
		self.grid.addWidget(self.combo_marker,3,2)
		self.grid.addWidget(self.btnMarkerColor,4,2)
		self.grid.addWidget(self.qle_linewidth,5,2)
		self.grid.addWidget(self.qle_alpha,6,2)
		self.grid.addWidget(self.btnDone,8,2)
		
		self.resize(300,500)
		self.setLayout(self.grid)
		self.setWindowTitle("Line editor")
		
		self.show()
		
	def donePressed(self):
		self.done(1)
		
	def getLine(self):
		return self.line
		
	def selectColor(self):
		
		"""Opens color selector and lets user select a color."""
		
		col = QtGui.QColorDialog.getColor(parent=self)
		col=tuple(asarray([col.red(),col.green(),col.blue()])/255.)
		
		self.line.set_color(col)
		self.updateColorButtonColor()
		
	def initColorButton(self):
		
		"""Initiates color button."""
		
		self.btnColor=QtGui.QPushButton('')
		self.btnColor.connect(self.btnColor, QtCore.SIGNAL('clicked()'), self.selectColor)
		
		self.updateColorButtonColor()
		
	def updateColorButtonColor(self):
		
		"""Updates the background color of the select color button."""
		
		colorRGB=matplotlib.colors.colorConverter.to_rgb(self.line.get_color())
		colorHex=matplotlib.colors.rgb2hex(colorRGB) 
		
		self.btnColor.setStyleSheet("background-color: "+colorHex)
		
	def selectMarkerColor(self):
		
		"""Opens color selector and lets user select a marker color."""
		
		col = QtGui.QColorDialog.getColor(parent=self)
		col=tuple(asarray([col.red(),col.green(),col.blue()])/255.)
		
		self.line.set_markerfacecolor(col)
		self.updateMarkerColorButtonColor()
		
	def initMarkerColorButton(self):
		
		"""Initiates marker color button."""
		
		self.btnMarkerColor=QtGui.QPushButton('')
		self.btnMarkerColor.connect(self.btnMarkerColor, QtCore.SIGNAL('clicked()'), self.selectMarkerColor)
		
		self.updateMarkerColorButtonColor()
		
	def updateMarkerColorButtonColor(self):
		
		"""Updates the background color of the select color button."""
		
		colorRGB=matplotlib.colors.colorConverter.to_rgb(self.line.get_markerfacecolor())
		colorHex=matplotlib.colors.rgb2hex(colorRGB) 
		
		self.btnMarkerColor.setStyleSheet("background-color: "+colorHex)	
		
	def initLinestyleCombo(self):
		
		"""Initiates linestyle combo box."""
		
		self.combo_linestyle = QtGui.QComboBox(self)
		
		for s in self.linestyles:
			self.combo_linestyle.addItem(s)
		
		self.combo_linestyle.activated[str].connect(self.selectLinestyle)  
		
		try:
			idx=self.linestyles.index(self.line.get_linestyle())
			self.combo_linestyle.setCurrentIndex(idx)
		except ValueError:
			pass
		
	def initMarkerCombo(self):
		
		"""Initiates marker combo box."""
		
		self.combo_marker = QtGui.QComboBox(self)
		
		for s in self.markers:
			self.combo_marker.addItem(s)
			
		self.combo_marker.activated[str].connect(self.selectMarker)  	
		
		try:
			idx=self.markers.index(self.line.get_marker())	
			self.combo_marker.setCurrentIndex(idx)
		except ValueError:
			pass
			
	def initMarkers(self):
		
		"""Initiates a list of possible markers."""
		
		markers = []
		for m in matplotlib.lines.Line2D.markers:
			try:
				if len(m) == 1 and m != ' ':
					markers.append(m)
			except TypeError:
				pass
		
		self.markers=markers
		
	def selectMarker(self,text):
		
		"""Sets marker according to text."""
		
		self.line.set_marker(str(text))
		
	def selectLinestyle(self,text):
		
		"""Sets marker according to text."""
		
		self.line.set_linestyle(str(text))
	
	def setLinewidth(self):
		
		"""Sets linewidth"""
		
		text=float(self.qle_linewidth.text())
		
		self.line.set_linewidth(text)
		
	def setAlpha(self):
		
		"""Sets alpha."""
		
		text=float(self.qle_alpha.text())
		
		self.line.set_alpha(text)
	
class tickDialog(QtGui.QDialog):
	
	def __init__(self,lbl,pos,parent):
		super(tickDialog,self).__init__(parent)		
		
		self.lbl=lbl
		self.pos=pos
		
		#Lbl
		self.lbl_lbl = QtGui.QLabel("Label:", self)
		self.qle_lbl = QtGui.QLineEdit(str(self.lbl))
		self.qle_lbl.editingFinished.connect(self.setLbl)
		
		#Pos
		self.lbl_pos = QtGui.QLabel("Position:", self)
		self.qle_pos = QtGui.QLineEdit(str(self.pos))
		self.qle_pos.editingFinished.connect(self.setPos)

		#Done button
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Layout
		self.grid = QtGui.QGridLayout()
		
		self.grid.addWidget(self.lbl_lbl,1,1)
		self.grid.addWidget(self.lbl_pos,2,1)
		
		self.grid.addWidget(self.qle_lbl,1,2)
		self.grid.addWidget(self.qle_pos,2,2)
		
		self.grid.addWidget(self.btnDone,4,2)
		
		self.resize(300,500)
		self.setLayout(self.grid)
		self.setWindowTitle("Tick editor")
		
		self.show()
		
	def donePressed(self):
		self.done(1)
	
	def setPos(self):	
		self.pos=float(self.qle_pos.text())
		
	def setLbl(self):	
		self.lbl=str(self.qle_lbl.text())	
		
	def getParms(self):
		return self.lbl,self.pos
	
class resetTickVecDialog(QtGui.QDialog):	
	
	def __init__(self,vec,parent):
		super(resetTickVecDialog,self).__init__(parent)		
		
		self.vec=vec
		
		#Minimum
		self.lbl_min = QtGui.QLabel("minimum:", self)
		self.qle_min = QtGui.QLineEdit(str(min(self.vec)))
		self.qle_min.editingFinished.connect(self.updateVec)
		
		#Maximum
		self.lbl_max = QtGui.QLabel("maximum:", self)
		self.qle_max = QtGui.QLineEdit(str(max(self.vec)))
		self.qle_max.editingFinished.connect(self.updateVec)
		
		#Steps
		self.lbl_steps = QtGui.QLabel("Number of steps:", self)
		self.qle_steps = QtGui.QLineEdit(str(len(self.vec)))
		self.qle_steps.editingFinished.connect(self.updateVec)
		
		#Scale
		self.lbl_scale = QtGui.QLabel("Scale:", self)
		self.combo_scale = QtGui.QComboBox(self)
		self.combo_scale.addItem("linear")
		self.combo_scale.addItem("logarithmic")
		self.combo_scale.activated[str].connect(self.selectScale)  
		
		#Done button
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Layout
		self.grid = QtGui.QGridLayout()
		
		self.grid.addWidget(self.lbl_min,1,1)
		self.grid.addWidget(self.lbl_max,2,1)
		self.grid.addWidget(self.lbl_steps,3,1)
		self.grid.addWidget(self.lbl_scale,4,1)
		
		self.grid.addWidget(self.qle_min,1,2)
		self.grid.addWidget(self.qle_max,2,2)
		self.grid.addWidget(self.qle_steps,3,2)
		self.grid.addWidget(self.combo_scale,4,2)
		
		self.grid.addWidget(self.btnDone,6,2)
		
		self.resize(300,500)
		self.setLayout(self.grid)
		self.setWindowTitle("Tick vector editor")
		
		self.show()
	
	def selectScale(self,text):
		self.updateVec()
	
	def updateVec(self):
		
		minv=float(self.qle_min.text())
		maxv=float(self.qle_max.text())
		steps=int(self.qle_steps.text())
		scale=str(self.combo_scale.currentText())
		
		if scale=="linear":
			self.vec=linspace(minv,maxv,steps)
		elif scale=="logarithmic":
			self.vec=logspace(minv,maxv,steps)
			
	def getVec(self):
		return self.vec
		
	def donePressed(self):
		self.done(1)	
		
		
		
		
		