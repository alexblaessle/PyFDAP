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

#! /usr/bin/python

#=====================================================================================================================================
#Description
#=====================================================================================================================================

#PyFRAP is a free Python-based software to analyze FDAP measurements. This file contains the main GUI helping to access the different 
#modules and classes:
#embryo.py: Class for FDAP datasets
#molecule.py: Class for collections of FDAP measurements
#pyfdap_fit_module.py: Module containing all necessary functions for fitting of FDAP data
#pyfdap_img_module.py: Module containing all necessary functions for FDAP data image analysis
#pyfdap_misc_module.py: Module containing some useful functions that keep PyFDAP going
#pyfdap_subwin.py: Module containing all pop-up, dialog and thread classes needed for PyFDAP
#pyfdap_term.py: Interactive Console class customized for PyFDAP needs

#=====================================================================================================================================
#Importing Packages and Modules
#=====================================================================================================================================

#Misc
import sys
from numpy import *
import code
import time
import os, os.path
import copy as cpy
import functools
import warnings
import Tkinter
import FileDialog

#PyFDP Modules
import pyfdap_img_module as pyfdap_img
import pyfdap_misc_module as pyfdap_misc
import pyfdap_stats_module as pyfdap_stats
import pyfdap_fit_module as pyfdap_fit
import pyfdap_subwin
import pyfdap_plot_dialogs
from pyfdap_term import *
from embryo import *
from molecule import *
from pyfdap_conf import *

#QT
from PyQt4 import QtGui, QtCore

#Matplotlib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
#from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as  NavigationToolbar
from matplotlib.figure import Figure

#=====================================================================================================================================
#Main Simulation window
#=====================================================================================================================================

class pyfdp(QtGui.QMainWindow):
	def __init__(self,ignWarnings=True,redirect=True, parent=None):
		QtGui.QWidget.__init__(self, parent)
	
		#-------------------------------------------
		#Initializing window
		#-------------------------------------------
		
		self.setWindowTitle('PyFDAP')
		self.setMinimumSize(400,300) 
		self.resize(1500,1000)
		self.dpi = 100
		self.version="1.1.2"
		self.website="http://people.tuebingen.mpg.de/mueller-lab"
		self.pyfdap_dir=os.getcwd()
		self.ignWarnings=ignWarnings
		self.redirect=redirect
		
		#print self.redirect, self.ignWarnings
		
		#-------------------------------------------
		#Statusbar
		#-------------------------------------------
		
		self.statusBar().showMessage("Idle")
		
		#-------------------------------------------
		#Some variables
		#-------------------------------------------
		
		self.molecules=[]
		self.tab_axes=[]
		self.tab_figs=[]
		self.curr_embr=None
		self.curr_fit=None
		self.curr_bkgd=None
		self.curr_embr_node=None
		self.curr_mol_node=None
		self.curr_fit_node=None
		self.lastopen=os.getcwd()
		self.copied_embr=None
		
		#-------------------------------------------
		#Menubar entries
		#-------------------------------------------
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#File
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		#Exit button
		newmolecule = QtGui.QAction('New Molecule', self)
		self.connect(newmolecule, QtCore.SIGNAL('triggered()'), self.add_molecule)
		
		loadmolecule = QtGui.QAction('Open Molecule', self)
		self.connect(loadmolecule, QtCore.SIGNAL('triggered()'), self.load_molecule)
		
		savemolecule = QtGui.QAction('Save Molecule', self)
		self.connect(savemolecule, QtCore.SIGNAL('triggered()'), self.save_molecule)
		
		exit = QtGui.QAction('Exit', self)
		exit.setShortcut('Ctrl+Q')	
		self.connect(exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Edit
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		removemolecule = QtGui.QAction('Remove Molecule', self)
		#removemolecule.setShortcut('Ctrl+Alt+R')
		QtGui.QShortcut(QtGui.QKeySequence("Ctrl+Alt+R"), self, self.delete_molecule)
		self.connect(removemolecule, QtCore.SIGNAL('triggered()'), self.delete_molecule)
		
		copymolecule = QtGui.QAction('Copy Molecule', self)
		self.connect(copymolecule, QtCore.SIGNAL('triggered()'), self.copy_molecule)
		
		editmolecule = QtGui.QAction('Edit Molecule', self)
		self.connect(editmolecule, QtCore.SIGNAL('triggered()'), self.edit_molecule)
		
		exportplot = QtGui.QAction('Export Plot', self)	
		self.connect(exportplot, QtCore.SIGNAL('triggered()'), self.export_plot)
		
		editplot = QtGui.QAction('Edit current Plot', self)	
		self.connect(editplot, QtCore.SIGNAL('triggered()'), self.edit_plot)
		
		exportplotseries = QtGui.QAction('Export Plot Series', self)	
		self.connect(exportplotseries, QtCore.SIGNAL('triggered()'), self.export_plot_series)
		
		exportmovie = QtGui.QAction('Export Movie', self)	
		self.connect(exportmovie, QtCore.SIGNAL('triggered()'), self.export_movie)
		
		exportembryo = QtGui.QAction('Export Embryo to csv-file', self)	
		self.connect(exportembryo, QtCore.SIGNAL('triggered()'), self.export_embryo_csv)
		
		exportmolecule = QtGui.QAction('Export Molecule to csv-file', self)	
		self.connect(exportmolecule, QtCore.SIGNAL('triggered()'), self.export_molecule_csv)
		
		exportfit = QtGui.QAction('Export Fit to csv-file', self)	
		self.connect(exportfit, QtCore.SIGNAL('triggered()'), self.export_fit_to_csv)
		
		exporterror = QtGui.QAction('Export errorbar plot to csv-file', self)	
		self.connect(exporterror, QtCore.SIGNAL('triggered()'), self.export_errorbar_to_csv)
		
		exportselobjtocsv = QtGui.QAction('Export selected object to csv-file', self)	
		self.connect(exportselobjtocsv, QtCore.SIGNAL('triggered()'), self.export_sel_obj_to_csv)
		
		exportselfitstocsv = QtGui.QAction('Export selected fits to csv-file', self)	
		self.connect(exportselfitstocsv, QtCore.SIGNAL('triggered()'), self.export_sel_fits_to_csv)
		
		memusage = QtGui.QAction('Print object memory usage', self)
		self.connect(memusage, QtCore.SIGNAL('triggered()'), self.print_mem_usage)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#View
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		showconsole = QtGui.QAction('Show Console', self)	
		self.connect(showconsole, QtCore.SIGNAL('triggered()'), self.show_console)
		
		hideconsole = QtGui.QAction('Hide Console', self)	
		self.connect(hideconsole, QtCore.SIGNAL('triggered()'), self.hide_console)
	
		showproplist = QtGui.QAction('Show Property Column', self)	
		self.connect(showproplist, QtCore.SIGNAL('triggered()'), self.show_proplist)
		
		hideproplist = QtGui.QAction('Hide Property Column', self)	
		self.connect(hideproplist, QtCore.SIGNAL('triggered()'), self.hide_proplist)
		
		showplottab = QtGui.QAction('Show Plot-tab', self)	
		self.connect(showplottab, QtCore.SIGNAL('triggered()'), self.show_plottab)
		
		hideplottab = QtGui.QAction('Hide Plot-tab', self)	
		self.connect(hideplottab, QtCore.SIGNAL('triggered()'), self.hide_plottab)
	
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Data analysis
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		analyze = QtGui.QAction('Analyze Embryo', self)
		self.connect(analyze, QtCore.SIGNAL('triggered()'), self.analyze_dataset)
		
		analyzebkgds = QtGui.QAction('Analyze Background Data Sets', self)
		self.connect(analyzebkgds, QtCore.SIGNAL('triggered()'), self.analyze_bkgds)
		
		analyzeall = QtGui.QAction('Analyze Molecule', self)
		self.connect(analyzeall, QtCore.SIGNAL('triggered()'), self.analyze_all)
		
		newembryo = QtGui.QAction('New Embryo', self)
		self.connect(newembryo, QtCore.SIGNAL('triggered()'), self.add_embryo)
		
		loadembryo = QtGui.QAction('Load Embryo', self)
		self.connect(loadembryo, QtCore.SIGNAL('triggered()'), self.load_embryo)
		
		loadembryocsv = QtGui.QAction('Load Embryos from CSV', self)
		self.connect(loadembryocsv, QtCore.SIGNAL('triggered()'), self.load_embryos_from_csv)
		
		saveembryo = QtGui.QAction('Save Embryo', self)
		self.connect(saveembryo, QtCore.SIGNAL('triggered()'), self.save_embryo)
		
		removeembryo = QtGui.QAction('Remove Embryo', self)
		self.connect(removeembryo, QtCore.SIGNAL('triggered()'), self.delete_embryo)
	
		copyany = QtGui.QAction('Copy Any', self)
		QtGui.QShortcut(QtGui.QKeySequence("Ctrl+C"), self, self.copy_any)
		self.connect(copyany, QtCore.SIGNAL('triggered()'), self.copy_any)
		
		pasteany = QtGui.QAction('Paste Any', self)
		QtGui.QShortcut(QtGui.QKeySequence("Ctrl+V"), self, self.paste_any)
		self.connect(pasteany, QtCore.SIGNAL('triggered()'), self.paste_any)
		
		copyembryo = QtGui.QAction('Copy Embryo', self)
		self.connect(copyembryo, QtCore.SIGNAL('triggered()'), self.copy_embryo)
		
		pasteembryo = QtGui.QAction('Paste Embryo', self)
		self.connect(pasteembryo, QtCore.SIGNAL('triggered()'), self.paste_embryo)
		
		editdataset = QtGui.QAction('Edit Embryo', self)
		self.connect(editdataset, QtCore.SIGNAL('triggered()'), self.edit_dataset)
		
		editpre = QtGui.QAction('Edit Pre', self)
		self.connect(editpre, QtCore.SIGNAL('triggered()'), self.edit_pre)
		
		editnoise = QtGui.QAction('Edit Noise', self)
		self.connect(editnoise, QtCore.SIGNAL('triggered()'), self.edit_noise)
		
		editignored = QtGui.QAction('Edit Ignored Frames', self)
		self.connect(editignored, QtCore.SIGNAL('triggered()'), self.select_ignored_frames)
		
		editthresh = QtGui.QAction('Edit Threshhold', self)
		self.connect(editthresh, QtCore.SIGNAL('triggered()'), self.edit_thresh)
		
		plotdata = QtGui.QAction('Plot analyzed time series', self)
		self.connect(plotdata, QtCore.SIGNAL('triggered()'), self.plot_data_timeseries)
		
		plotdatabkgds = QtGui.QAction('Plot analyzed time series and background values', self)
		self.connect(plotdatabkgds, QtCore.SIGNAL('triggered()'), self.plot_data_bkgd)
		
		plotallextdata = QtGui.QAction('Plot all extracelluar time series and background values', self)
		self.connect(plotallextdata, QtCore.SIGNAL('triggered()'), self.plot_all_ext_data)
		
		plotallintdata = QtGui.QAction('Plot all intracellular time series and background values', self)
		self.connect(plotallintdata, QtCore.SIGNAL('triggered()'), self.plot_all_int_data)
		
		plotallslicedata = QtGui.QAction('Plot all slice time series and background values', self)
		self.connect(plotallslicedata, QtCore.SIGNAL('triggered()'), self.plot_all_slice_data)
		
		plotallextdataign = QtGui.QAction('Plot all extracelluar time series and background values (including ignored tps)', self)
		self.connect(plotallextdataign, QtCore.SIGNAL('triggered()'), self.plot_all_ext_data_ign)
		
		plotallintdataign = QtGui.QAction('Plot all intracellular time series and background values (including ignored tps)', self)
		self.connect(plotallintdataign, QtCore.SIGNAL('triggered()'), self.plot_all_int_data_ign)
		
		plotallslicedataign = QtGui.QAction('Plot all slice time series and background values (including ignored tps)', self)
		self.connect(plotallslicedataign, QtCore.SIGNAL('triggered()'), self.plot_all_slice_data_ign)
		
		plotembryosliceimgs = QtGui.QAction('Plot slice images', self)
		self.connect(plotembryosliceimgs, QtCore.SIGNAL('triggered()'), self.plot_embryo_slice_imgs)
		
		plotembryoextimgs = QtGui.QAction('Plot extracelluar images', self)
		self.connect(plotembryoextimgs, QtCore.SIGNAL('triggered()'), self.plot_embryo_ext_imgs)
		
		plotembryointimgs = QtGui.QAction('Plot intracellular images', self)
		self.connect(plotembryointimgs, QtCore.SIGNAL('triggered()'), self.plot_embryo_int_imgs)
		
		plotembryoslicemasks = QtGui.QAction('Plot slice mask', self)
		self.connect(plotembryoslicemasks, QtCore.SIGNAL('triggered()'), self.plot_embryo_masks_embryo)
		
		plotembryoextmasks = QtGui.QAction('Plot extracelluar mask', self)
		self.connect(plotembryoextmasks, QtCore.SIGNAL('triggered()'), self.plot_embryo_masks_ext)
		
		plotembryointmasks = QtGui.QAction('Plot intracellular mask', self)
		self.connect(plotembryointmasks, QtCore.SIGNAL('triggered()'), self.plot_embryo_masks_int)
		
		plotbkgdtimeseries = QtGui.QAction('Plot background timeseries', self)
		self.connect(plotbkgdtimeseries, QtCore.SIGNAL('triggered()'), self.plot_bkgd_timeseries)
		
		plotbkgdsliceimgs = QtGui.QAction('Plot slice images', self)
		self.connect(plotbkgdsliceimgs, QtCore.SIGNAL('triggered()'), self.plot_bkgd_slice_imgs)
		
		plotbkgdextimgs = QtGui.QAction('Plot extracelluar images', self)
		self.connect(plotbkgdextimgs, QtCore.SIGNAL('triggered()'), self.plot_bkgd_ext_imgs)
		
		plotbkgdintimgs = QtGui.QAction('Plot intracellular images', self)
		self.connect(plotbkgdintimgs, QtCore.SIGNAL('triggered()'), self.plot_bkgd_int_imgs)
		
		plotbkgdslicemasks = QtGui.QAction('Plot slice mask', self)
		self.connect(plotbkgdslicemasks, QtCore.SIGNAL('triggered()'), self.plot_bkgd_masks_embryo)
		
		plotbkgdextmasks = QtGui.QAction('Plot extracelluar mask', self)
		self.connect(plotbkgdextmasks, QtCore.SIGNAL('triggered()'), self.plot_bkgd_masks_ext)
		
		plotbkgdintmasks = QtGui.QAction('Plot intracellular mask', self)
		self.connect(plotbkgdintmasks, QtCore.SIGNAL('triggered()'), self.plot_bkgd_masks_int)
			
		addbkgd = QtGui.QAction('Add background dataset', self)
		self.connect(addbkgd, QtCore.SIGNAL('triggered()'), self.add_bkgd)
		
		removebkgd = QtGui.QAction('Remove background dataset', self)
		self.connect(removebkgd, QtCore.SIGNAL('triggered()'), self.delete_bkgd)
		
		#Edit Background dataset
		editbkgd = QtGui.QAction('Edit background dataset', self)
		self.connect(editbkgd, QtCore.SIGNAL('triggered()'), self.edit_bkgd)
		
		editbkgdpre = QtGui.QAction('Edit Pre', self)
		self.connect(editbkgdpre, QtCore.SIGNAL('triggered()'), self.edit_bkgd_pre)
		
		copybkgd = QtGui.QAction('Copy Background', self)
		self.connect(copybkgd, QtCore.SIGNAL('triggered()'), self.copy_bkgd)
		
		pastebkgd = QtGui.QAction('Paste Background', self)
		self.connect(pastebkgd, QtCore.SIGNAL('triggered()'), self.paste_bkgd)
		
		editignoredbkgd = QtGui.QAction('Edit Ignored Frames', self)
		self.connect(editignoredbkgd, QtCore.SIGNAL('triggered()'), self.select_ignored_frames)
		
		loadbkgdcsv = QtGui.QAction('Load Bkgds from CSV', self)
		self.connect(loadbkgdcsv, QtCore.SIGNAL('triggered()'), self.load_bkgds_from_csv)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Fitting
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		addfit = QtGui.QAction('New fit', self)
		self.connect(addfit, QtCore.SIGNAL('triggered()'), self.add_fit)
		
		removefit = QtGui.QAction('Remove fit', self)
		self.connect(removefit, QtCore.SIGNAL('triggered()'), self.delete_fit)
		
		editfit = QtGui.QAction('Edit fit', self)
		self.connect(editfit, QtCore.SIGNAL('triggered()'), self.edit_fit)
		
		editmultfit = QtGui.QAction('Edit multiple fits', self)
		self.connect(editmultfit, QtCore.SIGNAL('triggered()'), self.edit_mult_fit)
		
		copyfit = QtGui.QAction('Copy fit', self)
		self.connect(copyfit, QtCore.SIGNAL('triggered()'), self.copy_fit)
		
		copyfitforall = QtGui.QAction('Copy fit into all embryos', self)
		self.connect(copyfitforall, QtCore.SIGNAL('triggered()'), self.copy_fit_to_all)
		
		performfit = QtGui.QAction('Perform fit', self)
		self.connect(performfit, QtCore.SIGNAL('triggered()'), self.perform_fit)
		
		fitallseries = QtGui.QAction('Fit all data series', self)
		self.connect(fitallseries, QtCore.SIGNAL('triggered()'), self.fit_all_series)
		
		fitall = QtGui.QAction('Perform all fits in molecule', self)
		self.connect(fitall, QtCore.SIGNAL('triggered()'), self.perform_fits_molecule)
		
		plotfit = QtGui.QAction('Plot fit', self)
		self.connect(plotfit, QtCore.SIGNAL('triggered()'), self.plot_fit)
		
		plottrackfit = QtGui.QAction('Plot fitting progress', self)
		self.connect(plottrackfit, QtCore.SIGNAL('triggered()'), self.plot_track_fit)
		
		plotextrapfit = QtGui.QAction('Plot extrapolated fit', self)
		self.connect(plotextrapfit, QtCore.SIGNAL('triggered()'), self.plot_extrap_fit)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Stats
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		sumupmol = QtGui.QAction('Average Molecule', self)
		self.connect(sumupmol, QtCore.SIGNAL('triggered()'), self.sumup_molecule)
		
		normplot = QtGui.QAction('Plot normed average fit', self)
		self.connect(normplot, QtCore.SIGNAL('triggered()'), self.norm_data_fit_plot)
		
		unpinplot = QtGui.QAction('Plot average fit', self)
		self.connect(unpinplot, QtCore.SIGNAL('triggered()'), self.unpin_data_fit_plot)
		
		barks = QtGui.QAction('Plot ks by fit', self)
		self.connect(barks, QtCore.SIGNAL('triggered()'), self.plot_ks_by_fit)
	
		barks = QtGui.QAction('Plot ks by fit', self)
		self.connect(barks, QtCore.SIGNAL('triggered()'), self.plot_ks_by_fit)
		
		bary0s = QtGui.QAction('Plot y0s by fit', self)
		self.connect(bary0s, QtCore.SIGNAL('triggered()'), self.plot_ynaughts_by_fit)
		
		barc0s = QtGui.QAction('Plot c0s by fit', self)
		self.connect(barc0s, QtCore.SIGNAL('triggered()'), self.plot_cnaughts_by_fit)
		
		barall = QtGui.QAction('Plot all parameters by fit', self)
		self.connect(barall, QtCore.SIGNAL('triggered()'), self.plot_all_by_fit)
		
		histk = QtGui.QAction('Histogram ks', self)
		self.connect(histk, QtCore.SIGNAL('triggered()'), self.hist_k)
		
		histtaumin = QtGui.QAction('Histogram halflives', self)
		self.connect(histtaumin, QtCore.SIGNAL('triggered()'), self.hist_taumin)
		
		avgFs  = QtGui.QAction('Compute average correction factors', self)
		self.connect(avgFs, QtCore.SIGNAL('triggered()'), self.compute_av_corr_Fs)
		
		ttest  = QtGui.QAction('Perform standard t-test', self)
		self.connect(ttest, QtCore.SIGNAL('triggered()'), self.perform_ttest)
		
		ttest_welch  = QtGui.QAction('Perform Welchs t-test', self)
		self.connect(ttest_welch, QtCore.SIGNAL('triggered()'), self.perform_welch)
		
		wilcoxontest  = QtGui.QAction('Perform Wilcoxon test', self)
		self.connect(wilcoxontest, QtCore.SIGNAL('triggered()'), self.perform_wilcoxon)
		
		mannwhitneytest  = QtGui.QAction('Perform Mann-Whitney test', self)
		self.connect(mannwhitneytest, QtCore.SIGNAL('triggered()'), self.perform_mann_whitney)
		
		shirapotest  = QtGui.QAction('Perform Shirapo test', self)
		self.connect(shirapotest, QtCore.SIGNAL('triggered()'), self.perform_sharipo)
		
		compnormplot = QtGui.QAction('Plot normed average fit', self)
		self.connect(compnormplot, QtCore.SIGNAL('triggered()'), self.compare_molecule_norm)
		
		compunpinplot = QtGui.QAction('Plot average fit', self)
		self.connect(compunpinplot, QtCore.SIGNAL('triggered()'), self.compare_molecule_unpin)
		
		compbarks = QtGui.QAction('Plot ks by molecule', self)
		self.connect(compbarks, QtCore.SIGNAL('triggered()'), self.compare_molecule_ks)
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#Help
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		about  = QtGui.QAction('About PyFDAP', self)
		self.connect(about, QtCore.SIGNAL('triggered()'), self.show_about)
		
		self.menubar = self.menuBar()
		
		self.file_mb = self.menubar.addMenu('&File')
		self.file_mb.addAction(newmolecule)
		self.file_mb.addAction(loadmolecule)
		self.file_recent_mb=self.file_mb.addMenu('&Open recent')
		self.file_mb.addAction(savemolecule)
		self.file_mb.addAction(exit)
		
		self.edit_mb = self.menubar.addMenu('&Edit')
		self.edit_mb.addAction(copymolecule)
		self.edit_mb.addAction(editmolecule)
		self.edit_mb.addAction(removemolecule)
		self.edit_plot_mb=self.edit_mb.addMenu('&Plotting')
		self.edit_plot_mb.addAction(editplot)
		self.edit_export_mb=self.edit_mb.addMenu('&Export')
		self.edit_export_mb.addAction(exportplot)
		self.edit_export_mb.addAction(exportplotseries)
		self.edit_export_mb.addAction(exportmovie)
		self.edit_export_mb.addAction(exportembryo)
		self.edit_export_mb.addAction(exportmolecule)
		self.edit_export_mb.addAction(exportfit)
		self.edit_export_mb.addAction(exporterror)
		self.edit_export_mb.addAction(exportselobjtocsv)
		self.edit_export_mb.addAction(exportselfitstocsv)
		self.edit_mb.addAction(memusage)
		
		self.view_mb = self.menubar.addMenu('&View')
		self.view_console_mb = self.view_mb.addMenu('&Console')
		self.view_console_mb.addAction(showconsole)
		self.view_console_mb.addAction(hideconsole)
		self.view_proplist_mb = self.view_mb.addMenu('&Property Column')
		self.view_proplist_mb.addAction(showproplist)
		self.view_proplist_mb.addAction(hideproplist)
		self.view_plottab_mb = self.view_mb.addMenu('&Plot Tabs')
		self.view_plottab_mb.addAction(showplottab)
		self.view_plottab_mb.addAction(hideplottab)
		
		self.data_mb = self.menubar.addMenu('&Data Analyses')
		self.data_data_mb=self.data_mb.addMenu('&Embryo')
		self.data_data_mb.addAction(newembryo)
		self.data_data_mb.addAction(loadembryo)
		self.data_data_mb.addAction(loadembryocsv)
		self.data_data_mb.addAction(saveembryo)
		self.data_data_mb.addAction(removeembryo)
		self.data_data_mb.addAction(copyembryo)
		self.data_data_mb.addAction(pasteembryo)
		self.data_data_mb.addAction(editdataset)
		self.data_data_mb.addAction(editpre)
		self.data_data_mb.addAction(editnoise)
		self.data_data_mb.addAction(editignored)
		self.data_data_mb.addAction(editthresh)
				
		self.data_bkgd_mb=self.data_mb.addMenu('&Background Datasets')
		self.data_bkgd_mb.addAction(addbkgd)
		self.data_bkgd_mb.addAction(removebkgd)
		self.data_bkgd_mb.addAction(editbkgd)
		self.data_bkgd_mb.addAction(editbkgdpre)
		self.data_bkgd_mb.addAction(copybkgd)
		self.data_bkgd_mb.addAction(pastebkgd)
		self.data_bkgd_mb.addAction(editignored)
		self.data_bkgd_mb.addAction(loadbkgdcsv)
		
		self.data_analysis_mb=self.data_mb.addMenu('&Analysis')
		
		self.data_analysis_mb.addAction(analyzeall)
		self.data_analysis_mb.addAction(analyze)
		self.data_analysis_mb.addAction(analyzebkgds)
		
		self.data_plotting_mb=self.data_mb.addMenu('&Plotting')
		self.data_plotting_main_mb=self.data_plotting_mb.addMenu('&Main Dataset')
		self.data_plotting_main_mb.addAction(plotdata)
		self.data_plotting_main_mb.addAction(plotdatabkgds)
		self.data_plotting_main_mb.addAction(plotallextdata)
		self.data_plotting_main_mb.addAction(plotallintdata)
		self.data_plotting_main_mb.addAction(plotallslicedata)
		self.data_plotting_main_mb.addAction(plotallextdataign)
		self.data_plotting_main_mb.addAction(plotallintdataign)
		self.data_plotting_main_mb.addAction(plotallslicedataign)
		self.data_plotting_main_mb.addAction(plotembryosliceimgs)
		self.data_plotting_main_mb.addAction(plotembryoextimgs)
		self.data_plotting_main_mb.addAction(plotembryointimgs)
		self.data_plotting_main_mb.addAction(plotembryoslicemasks)
		self.data_plotting_main_mb.addAction(plotembryoextmasks)
		self.data_plotting_main_mb.addAction(plotembryointmasks)
		
		self.data_plotting_bkgd_mb=self.data_plotting_mb.addMenu('&Background Datasets')
		self.data_plotting_bkgd_mb.addAction(plotbkgdtimeseries)
		self.data_plotting_bkgd_mb.addAction(plotbkgdsliceimgs)
		self.data_plotting_bkgd_mb.addAction(plotbkgdextimgs)
		self.data_plotting_bkgd_mb.addAction(plotbkgdintimgs)
		self.data_plotting_bkgd_mb.addAction(plotbkgdslicemasks)
		self.data_plotting_bkgd_mb.addAction(plotbkgdextmasks)
		self.data_plotting_bkgd_mb.addAction(plotbkgdintmasks)
		
		self.fit_mb = self.menubar.addMenu('&Fitting')
		self.fit_fits_mb=self.fit_mb.addMenu('&Fits')
		self.fit_fitting_mb=self.fit_mb.addMenu('&Perform Fits')
		
		self.fit_fits_mb.addAction(addfit)
		self.fit_fits_mb.addAction(removefit)
		self.fit_fits_mb.addAction(editfit)
		self.fit_fits_mb.addAction(editmultfit)
		self.fit_fits_mb.addAction(copyfit)
		self.fit_fits_mb.addAction(copyfitforall)
		self.fit_fitting_mb.addAction(performfit)
		self.fit_fitting_mb.addAction(fitallseries)
		self.fit_fitting_mb.addAction(fitall)
		self.fit_plot_mb=self.fit_mb.addMenu('&Plotting')
		self.fit_plot_mb.addAction(plotfit)
		self.fit_plot_mb.addAction(plottrackfit)
		self.fit_plot_mb.addAction(plotextrapfit)
	
		
		self.stats_mb = self.menubar.addMenu('&Statistics')
		self.stats_mb.addAction(sumupmol)
		self.stats_mb.addAction(avgFs)
		self.stats_test_mb = self.stats_mb.addMenu('&Tests')
		self.stats_test_mb.addAction(ttest)
		self.stats_test_mb.addAction(ttest_welch)
		self.stats_test_mb.addAction(wilcoxontest)
		self.stats_test_mb.addAction(mannwhitneytest)
		self.stats_test_mb.addAction(shirapotest)
		self.stats_plot_mb = self.stats_mb.addMenu('&Plotting')
		self.stats_plot_single_mb = self.stats_plot_mb.addMenu('&Single Molecule')
		self.stats_plot_comp_mb = self.stats_plot_mb.addMenu('&Molecule comparison')
		self.stats_plot_single_mb.addAction(unpinplot)
		self.stats_plot_single_mb.addAction(normplot)
		self.stats_plot_single_mb.addAction(barks)
		self.stats_plot_single_mb.addAction(bary0s)
		self.stats_plot_single_mb.addAction(barc0s)
		self.stats_plot_single_mb.addAction(barall)
		self.stats_plot_single_mb.addAction(histk)
		self.stats_plot_single_mb.addAction(histtaumin)
		
		self.stats_plot_comp_mb.addAction(compnormplot)
		self.stats_plot_comp_mb.addAction(compunpinplot)
		self.stats_plot_comp_mb.addAction(compbarks)
		
		self.help_mb = self.menubar.addMenu('&Help')
		self.help_mb.addAction(about)
				
		#-------------------------------------------
		#Embryo list
		#-------------------------------------------

		self.embryos_list=QtGui.QTreeWidget()
		self.embryos_list.setHeaderLabels(["Embryo","Analyzed","Fitted"])
		self.embryos_list.setColumnWidth(0,200)
		self.embryos_list.setColumnWidth(1,75)
		self.embryos_list.setColumnWidth(2,75)
		self.embryos_list.itemActivated.connect(self.show_embryo_props) 
		
		#-------------------------------------------
		#Property bar
		#-------------------------------------------
		
		self.prop_list=QtGui.QListWidget()
		
		#-------------------------------------------
		#Console
		#-------------------------------------------
			
		self.console = PyInterp(self,ignWarnings=self.ignWarnings,redirect=self.redirect)
		self.console.initInterpreter(locals())
		
		#-------------------------------------------
		#Splitter
		#-------------------------------------------
		
		#Creating splitters
		self.splitter_hor = QtGui.QSplitter(QtCore.Qt.Horizontal)
		self.splitter_ver = QtGui.QSplitter(QtCore.Qt.Vertical)
		
		#Creating widgets
		self.embryos_frame   = QtGui.QFrame()
		self.embryos_frame   = self.create_frame(self.embryos_frame)
		
		self.prop_frame   = QtGui.QFrame()
		self.prop_frame   = self.create_frame(self.prop_frame)
		
		self.term_frame   = QtGui.QFrame()
		self.term_frame   = self.create_frame(self.term_frame)
		
		#-------------------------------------------
		#Plotting tabs
		#-------------------------------------------
		
		self.plot_tabs = QtGui.QTabWidget()
		self.plot_tabs.setTabsClosable(False)
		self.plot_tabs.currentChanged.connect(self.curr_tab_changed)
		self.plot_tabs.tabCloseRequested.connect(self.curr_tab_closed)
		
		#-------------------------------------------
		#Final Layout
		#-------------------------------------------
		
		#Adding widgets to splitters	
		self.splitter_hor.addWidget(self.embryos_list)
		self.splitter_hor.addWidget(self.plot_tabs)
		self.splitter_hor.addWidget(self.prop_list)
		self.splitter_ver.addWidget(self.splitter_hor)
		self.splitter_ver.addWidget(self.console)
		
		#Setting default sizes of splitters
		self.splitter_ver.setSizes([750,250])
		self.splitter_hor.setSizes([350,900,250])
		
		#Connecting splitter movement to figure size adjustment
		self.splitter_hor.splitterMoved.connect(self.adjust_canvas)
		self.splitter_ver.splitterMoved.connect(self.adjust_canvas)

		self.curr_tab=QtGui.QWidget()
		self.first_tab=self.plot_tabs.addTab(self.curr_tab,"PlotTab")
		
		#Load config file
		self.init_conf()
		
		self.bin_width_halfilfe_min=10
		
		self.setCentralWidget(self.splitter_ver)
		QtGui.QApplication.setStyle(QtGui.QStyleFactory.create('Cleanlooks'))
		self.show()
				
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Functions
	#----------------------------------------------------------------------------------------------------------------------------------------
	
	def closeEvent(self, event):
			
		reply = QtGui.QMessageBox.question(self, 'Message',"Are you sure you want to quit?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
	
		if reply == QtGui.QMessageBox.Yes:
			fn=self.pyfdap_dir+"/pyfdap_conf.pk"
			self.curr_conf.save_conf(fn)
			event.accept()
		else:
			event.ignore()
	
	def show_about(self):
		
		ret=pyfdap_subwin.about_dialog(self).exec_()
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Create frame
	
	def create_frame(self,frame):
		
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		frame.setBackgroundRole(QtGui.QPalette.Light)
		frame.setAutoFillBackground(True)        
		frame.setLineWidth (1)
		frame.setFrameShadow (frame.Sunken)
		
		return frame
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Recently opened
	
	def init_conf(self):
		
		fn=self.pyfdap_dir+"/pyfdap_conf.pk"
		self.curr_conf=pyfdap_conf()
		if os.path.isfile(fn):
			self.curr_conf=self.curr_conf.load_conf(fn)
		
			self.console.setHidden(self.curr_conf.term_hidden)
			self.prop_list.setHidden(self.curr_conf.prop_hidden)
			self.plot_tabs.setHidden(self.curr_conf.plot_hidden)
			self.splitter_ver.refresh()
			self.splitter_hor.refresh()
			self.adjust_canvas()
			
			self.add_recent_mbs()
		
	def append_recent(self,fn_new):
	
		if fn_new in self.curr_conf.recent:
			ind=self.curr_conf.recent.index(fn_new)
			self.curr_conf.recent.pop(ind)
			
		self.curr_conf.recent.insert(0,str(fn_new))
		self.file_recent_mb.clear()
		self.add_recent_mbs()
	
	def add_recent_mbs(self):
		self.recent_actions=[]
		for i in range(shape(self.curr_conf.recent)[0]): 
			if i>5:
				self.curr_conf.recent.pop(i)
			else:
				self.recent_actions.append(QtGui.QAction(self.curr_conf.recent[i], self))
				item=self.recent_actions[i]
				item.triggered.connect(functools.partial(self.open_molecule,self.curr_conf.recent[i]))
			
				self.file_recent_mb.addAction(item)	
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Add molecule
	
	def add_molecule(self):
		
		#Create new embryo object
		curr_name="newmolecule_1"
		mol=molecule(curr_name)
		
		#Check if there is already an embryo with the name newembryo_1
		in_list=1
		
		while in_list==1:
			#Go through all embryos
			for item in self.molecules:
				if item.name==curr_name:
					#If exists, split name in nmbr and name and add +1 to number
					nme,nmbr=curr_name.split("_")
					nmbr=str(int(nmbr)+1)
					curr_name=nme+"_"+nmbr
			
			#Actually in list
			in_list=0
		
		#Set new name
		mol.name=curr_name
		
		#Add embryo to list of embryos
		self.molecules.append(mol)
		
		#Add to embryo bar
		self.curr_node=QtGui.QTreeWidgetItem(self.embryos_list,[curr_name])
		self.curr_embryos=QtGui.QTreeWidgetItem(self.curr_node,["Embryos",'',''])
		self.curr_bkgds=QtGui.QTreeWidgetItem(self.curr_node,["Backgrounds",'',''])
		self.parent_node=None	
		
		self.curr_mol_node=self.curr_node
		
		self.curr_mol=mol
	
		self.embryos_list.setCurrentItem(self.curr_node)
		self.show_embryo_props()
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Add embryo
	
	def add_embryo(self):
		
		if self.curr_mol_node==None:
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
			
		#Create new embryo object
		curr_name="newembryo_1"
		emb=embryo(curr_name,"fdap")
		
		#Check if there is already an embryo with the name newembryo_1
		in_list=1
		
		while in_list==1:
			#Go through all embryos
			for item in self.curr_mol.embryos:
				if item.name==curr_name:
					#If exists, split name in nmbr and name and add +1 to number
					nme,nmbr=curr_name.split("_")
					nmbr=str(int(nmbr)+1)
					curr_name=nme+"_"+nmbr
			
			#Actually in list
			in_list=0
		
		#Set new name
		emb.name=curr_name
		
		#Add embryo to list of embryos
		self.curr_mol.embryos.append(emb)
		
		#Add to embryo bar
		self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[curr_name,"0","0"])
		self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'',''])
		self.curr_pre_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Pre",'0',''])
		self.curr_noise_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Noise",'0',''])
		
		self.parent_node=self.curr_mol	
		
		self.curr_embr=emb
	
		self.embryos_list.setCurrentItem(self.curr_embr_node)
		self.show_embryo_props()
		self.edit_dataset()
		self.edit_pre()
		self.edit_noise()
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Save molecule
	
	def save_molecule(self):
		
		#Check if anything is highlighted
		if self.embryos_list.currentItem()==None:
		
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return

		else:
			
			reply = QtGui.QMessageBox.question(self, 'Message',"Do you also want to save the image data into the molecule file?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

			if reply == QtGui.QMessageBox.Yes:
				
				#Grab filename
				fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen+"/"+self.curr_mol.name+".pk","*.pk",)
				fn_save=str(fn_save)
				self.lastopen=os.path.dirname(str(fn_save))
				
				#Save molecule
				if fn_save=='':
					pass
				else:
					self.curr_mol.save_molecule(fn_save)
		
			else:
				
				#Moving everything to temporary molecule
				if self.curr_conf.backup_to_file:
					self.fn_backup=self.lastopen+"/"+self.curr_mol.name+"_backup.pk"
					self.curr_mol.save_molecule(self.fn_backup)
					temp_mol=self.curr_mol
				if self.curr_conf.backup_to_mem:
					temp_mol=cpy.deepcopy(self.curr_mol)
				else:
					temp_mol=self.curr_mol
				
				#Make some temporary lists
				masks_embryo_list=[]
				masks_ext_list=[]
				masks_int_list=[]
				
				fits_list=[]
				
				bkgd_masks_embryo_list=[]
				bkgd_masks_ext_list=[]
				bkgd_masks_int_list=[]
				bkgd_vals_slice_list=[]
				
				for temp_emb in temp_mol.embryos:
					
					#Dumping mask and image data into temporay lists
					if temp_emb.masks_ext!=None:
						masks_embryo_list.append(list(temp_emb.masks_embryo))
						masks_ext_list.append(list(temp_emb.masks_ext))
						masks_int_list.append(list(temp_emb.masks_int))
					
						#Deleting all image data to reduce file size
						temp_emb.masks_embryo=[]
						temp_emb.masks_ext=None
						temp_emb.masks_int=None
						temp_emb.vals_slice=None
					
					#Deleting all track fit data
					fit_list=[]
					for fit in temp_emb.fits:
						
						#Dumping fit data in list
						fit_list.append([fit.save_track,fit.track_fit,fit.track_parms])	
						
						#Clearing all costly fit data
						fit.save_track=0
						fit.track_fit=[]
						fit.track_parms=[]
					
					fits_list.append(fit_list)	
						
				for bkgd in temp_mol.bkgds:
					
					#Dumping fit data in list
					if bkgd.masks_ext!=None:
						bkgd_masks_embryo_list.append(list(bkgd.masks_embryo))
						bkgd_masks_ext_list.append(list(bkgd.masks_ext))
						bkgd_masks_int_list.append(list(bkgd.masks_int))
						bkgd_vals_slice_list.append(list(bkgd.bkgd_vals_slice))
					
					#Clearing all costly background data
					bkgd.masks_embryo=[]
					bkgd.masks_ext=None
					bkgd.masks_int=None
					bkgd.bkgd_vals_slice=[]
				
				#Grab filename
				fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen+"/"+self.curr_mol.name+".pk","*.pk",)
				fn_save=str(fn_save)	
				self.lastopen=os.path.dirname(str(fn_save))
				
				#Save molecule
				if fn_save=='':
					pass
				else:
					temp_mol.save_molecule(fn_save)
				
				#Recovering data
				if self.curr_conf.backup_to_file:
					self.curr_mol=self.curr_mol.load_molecule(self.fn_backup)
					os.remove(self.fn_backup)
				if self.curr_conf.backup_to_mem:
					temp_mol=None
				
				#Mapping all data back
				if len(masks_embryo_list)>0:
					for i,emb in enumerate(self.curr_mol.embryos):
						
						emb.masks_embryo=masks_embryo_list[i]
						emb.masks_ext=masks_ext_list[i]
						emb.masks_int=masks_int_list[i]
						
						for j,fit in enumerate(emb.fits):
							fit.save_track=fits_list[i][j][0]
							fit.track_fit=fits_list[i][j][1]
							fit.track_parms=fits_list[i][j][2]
				
				if len(bkgd_masks_embryo_list)>0:
					for i,bkgd in enumerate(self.curr_mol.bkgds):
						bkgd.masks_embryo=bkgd_masks_embryo_list[i]
						bkgd.masks_ext=bkgd_masks_ext_list[i]
						bkgd.masks_int=bkgd_masks_int_list[i]
						bkgd.bkgd_vals_slice=bkgd_vals_slice_list[i]
						
		if fn_save!='':
			self.append_recent(fn_save)
		
		return True
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Save embryo
	
	def save_embryo(self):
		
		#Check if anything is highlighted
		if self.embryos_list.currentItem()==None:
		
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		else:
			if self.curr_node.parent()!=self.curr_embryos:
				#If not give error msg
				QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				return
			
			else:
				reply = QtGui.QMessageBox.question(self, 'Message',"Do you also want to save the image data into the embryo file?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
	
				if reply == QtGui.QMessageBox.Yes:
					
					fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,".pk",".pk")
					
					fn_save=str(fn_save)
					self.lastopen=os.path.dirname(str(fn_save))
					if fn_save=='':
						pass
					else:
					
						self.curr_embr.save_embryo(fn_save)
				
				else:
					
					#Temporarily saving the embryo
					
					###NOTE: This seems really slow for big objects, it is probably better to just backup all image data in numpy arrays and then dump the whole thing
					#temp_emb=cpy.deepcopy(self.curr_embr)
								
					#Dumping mask and image data into temporay lists
					masks_embryo=list(self.curr_embr.masks_embryo)
					masks_ext=list(self.curr_embr.masks_ext)
					masks_int=list(self.curr_embr.masks_int)
					
					#Deleting all image data to reduce file size
					self.curr_embr.masks_embryo=[]
					self.curr_embr.masks_ext=None
					self.curr_embr.masks_int=None
					
					fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,".pk",".pk")
					fn_save=str(fn_save)	
					self.lastopen=os.path.dirname(str(fn_save))
					if fn_save=='':
						pass
					else:
						self.curr_embr.save_embryo(fn_save)
	
					#Putting everything back
					self.curr_embr.masks_embryo=masks_embryo
					self.curr_embr.masks_ext=masks_ext
					self.curr_embr.masks_int=masks_int
						
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Load molecule
	
	def load_molecule(self):
		
		fn_load=QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.lastopen)
		if fn_load=='':
			return
		
		
		self.open_molecule(fn_load)
		
	def open_molecule(self,fn_load):
		self.lastopen=os.path.dirname(str(fn_load))
		self.append_recent(fn_load)
		
		#Create new molecule object
		curr_name="newmolecule_1"
		mol=molecule(curr_name)
		
		#Load into new molecule object
		mol=mol.load_molecule(fn_load)
		curr_name=mol.name
		
		#Check if there is already a molecule with the name newmolecule_1
		in_list=1
		
		while in_list==1:
			#Go through all embryos
			for item in self.molecules:
				if item.name==curr_name:
					ret=QtGui.QMessageBox.warning(self,self.tr("Warning"),self.tr("Molecule already exists, do you want to overwrite it?"),QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
					
					if ret == QtGui.QMessageBox.Yes:
						self.delete_specific_molecule(molecule=item)
						
					else:
						if "_" in curr_name:
							nme,nmbr=curr_name.split("_")
							nmbr=str(int(nmbr)+1)
						else:
							nme=curr_name
							nmbr=str(1)
						
						curr_name=curr_name+"_"+nmbr
			#Actually in list
			in_list=0
		
		#Set new name
		mol.name=curr_name
		
		#Add embryo to list of embryos
		self.molecules.append(mol)
		
		#Adding molecule to sidebar
		self.curr_node=QtGui.QTreeWidgetItem(self.embryos_list,[curr_name,"",""])
		self.curr_embryos=QtGui.QTreeWidgetItem(self.curr_node,["Embryos","",""])
		self.curr_bkgds=QtGui.QTreeWidgetItem(self.curr_node,["Backgrounds","",""])
		
		self.curr_mol_node=self.curr_node
		
		for emb in mol.embryos:
		
			#Add to embryo bar
			analyzed=str(0)
			fitted=str(0)
			if shape(emb.ext_av_data_d)[0]>1:
				analyzed=str(1)
			if shape(emb.fits)>0:
				for fit in emb.fits:
					if shape(fit.fit_av_d)[0]>1:	
						fitted=str(1)
				
			#Adding embryo to sidebar
			self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[emb.name,analyzed,fitted])
			self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'',''])
			self.curr_pre_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Pre",analyzed,''])
			self.curr_noise_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Noise",analyzed,''])
			
			#Add fits if they exist
			for fit in emb.fits:
				if shape(fit.fit_av_d)[0]>0:		
					QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','1'])
				else:	
					QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','0'])
						
			self.curr_embr=emb
			self.embryos_list.expandItem(self.curr_fits)
			
		#Add bkgds if they exist
		for bkgd in mol.bkgds:	
			self.add_bkgd_to_sidebar(bkgd=bkgd)
							
		self.embryos_list.expandItem(self.curr_mol_node)
		self.embryos_list.expandItem(self.curr_embryos)
		self.embryos_list.expandItem(self.curr_bkgds)
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Load embryo
	
	def load_embryo(self):
		
		fn_load=QtGui.QFileDialog.getOpenFileName(self, 'Open file', self.lastopen)
		if fn_load=='':
			return
		self.lastopen=os.path.dirname(str(fn_load))

		#Create new embryo object
		curr_name="newembryo_1"
		emb=embryo(curr_name,"default")
		
		#Load into new embryo object
		emb=emb.load_embryo(fn_load)
		curr_name=emb.name
		
		#Check if there is already an embryo with the name newembryo_1
		in_list=1
		
		while in_list==1:
			#Go through all embryos
			for item in self.curr_mol.embryos:
				if item.name==curr_name:
					ret=QtGui.QMessageBox.warning(self,self.tr("Warning"),self.tr("Embryo already exists, do you want to overwrite it?"),QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
					
					if ret == QtGui.QMessageBox.Yes:
						self.delete_specific_embryo(embryo=item)
						
					else:
						if "_" in curr_name:
							nme,nmbr=curr_name.split("_")
							nmbr=str(int(nmbr)+1)
						else:
							nme=curr_name
							nmbr=str(1)
						
						curr_name=curr_name+"_"+nmbr
			#Actually in list
			in_list=0
		
		#Set new name
		emb.name=curr_name
		
		#Add embryo to list of embryos
		self.curr_mol.embryos.append(emb)
		
		self.add_embryo_to_sidebar(emb=emb)
						
		self.curr_embr=emb
		self.embryos_list.expandItem(self.curr_embr_node)
		self.embryos_list.expandItem(self.curr_fits)
	
	def add_embryo_to_sidebar(self,emb=None):
		
		if emb==None:
			emb=self.curr_embr
		
		#Add to embryo bar
		analyzed=str(0)
		fitted=str(0)
		if shape(emb.ext_av_data_d)[0]>1:
			analyzed=str(1)
		if shape(emb.fits)>0:
			for fit in emb.fits:
				if shape(fit.fit_av_d)[0]>1:	
					fitted=str(1)
		
		#Adding embryo to sidebar
		self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[emb.name,analyzed,fitted])
		self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'',''])
		self.curr_pre_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Pre",analyzed,''])
		self.curr_noise_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Noise",analyzed,''])
		
		#Add fits if they exist
		for fit in emb.fits:
			if shape(fit.fit_av_d)[0]>0:		
				QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','1'])
			else:	
				QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','0'])
		
		return
		
	def add_bkgd_to_sidebar(self,bkgd=None):
		
		if bkgd==None:
			bkgd=self.curr_bkgd
		
		if shape(bkgd.bkgd_slice)!=None:
			self.curr_bkgd_node=QtGui.QTreeWidgetItem(self.curr_bkgds,[bkgd.name,'1',''])
			self.curr_bkgd_pre_node=QtGui.QTreeWidgetItem(self.curr_bkgd_node,["Pre",'1',''])
		else:
			self.curr_bkgd_node=QtGui.QTreeWidgetItem(self.curr_bkgds,[bkgd.name,'0',''])
			self.curr_bkgd_pre_node=QtGui.QTreeWidgetItem(self.curr_bkgd_node,["Pre",'0',''])
				
	def load_embryos_from_csv(self):
		
		csv_dialog=pyfdap_subwin.csv_read_dialog(self.curr_mol,self)
		
		if csv_dialog.exec_():
			embryos = csv_dialog.get_embryos()
		
		for emb in embryos:
			self.add_embryo_to_sidebar(emb=emb)
			
	def load_bkgds_from_csv(self):
		
		csv_dialog=pyfdap_subwin.csv_read_dialog_bkgds(self.curr_mol,self)
		
		if csv_dialog.exec_():
			bkgds = csv_dialog.get_bkgds()
		
		for bkgd in bkgds:
			self.add_bkgd_to_sidebar(bkgd=bkgd)		
				
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Show embryo or fit properties in property bar
	
	def show_embryo_props(self):
		
		if self.ignWarnings:
			warnings.catch_warnings()
			
			warnings.filterwarnings("ignore",category=DeprecationWarning)
			warnings.filterwarnings("ignore",category=FutureWarning)	
		
		#Clearing property list
		self.prop_list.clear()
			
		#Getting selected embryo
		self.selected_index=self.embryos_list.currentIndex()
		self.curr_node=self.embryos_list.currentItem()
		
		#Getting curr_mol_node, curr_embr_node and curr_fits and curr_bkgds and curr_embryos	
		if self.curr_node.parent()==None:
			#This is a molecule node
			self.curr_mol_node=self.curr_node
			self.curr_embr_node=None
			self.curr_bkgd_node=None
			
		else:		
			if self.curr_node.data(0,0).toString()=="Embryos" or self.curr_node.data(0,0).toString()=="Backgrounds":
				#This is an embryos/bkgds node
				self.curr_mol_node=self.curr_node.parent()
				self.curr_embr_node=None
				self.curr_bkgd_node=None
				
			elif self.curr_node.parent().data(0,0).toString()=="Embryos":
				#This is an embryo node
				self.curr_mol_node=self.curr_node.parent().parent()
				self.curr_embr_node=self.curr_node
				self.curr_bkgd_node=None
				
			elif self.curr_node.parent().parent().data(0,0).toString()=="Backgrounds":
				#This is a bkgd_pre node
				self.curr_mol_node=self.curr_node.parent().parent().parent()
				self.curr_embr_node=None
				self.curr_bkgd_node=self.curr_node.parent()
				
			elif self.curr_node.parent().data(0,0).toString()=="Backgrounds":
				#This is a bkgd node
				self.curr_mol_node=self.curr_node.parent().parent()
				self.curr_embr_node=None
				self.curr_bkgd_node=self.curr_node
			
			elif self.curr_node.data(0,0).toString()=="Fits" or self.curr_node.data(0,0).toString()=="Noise":
				#This is fits/noise node
				self.curr_mol_node=self.curr_node.parent().parent().parent()
				self.curr_embr_node=self.curr_node.parent()
				self.curr_bkgd_node=None
				
			elif self.curr_node.data(0,0).toString()=="Pre" and self.curr_node.parent().parent().data(0,0).toString()=="Embryos":
				#This is a pre node
				self.curr_mol_node=self.curr_node.parent().parent().parent()
				self.curr_embr_node=self.curr_node.parent()
				self.curr_bkgd_node=None
					
			elif self.curr_node.parent().data(0,0).toString()=="Fits":
				#This is a fit node
				self.curr_mol_node=self.curr_node.parent().parent().parent().parent()
				self.curr_embr_node=self.curr_node.parent().parent()
				self.curr_bkgd_node=None
		
		#By default set curr_obj=None
		self.curr_obj=None
		
		#Now define curr_fits and curr_bkgds
		if self.curr_embr_node!=None:
			self.curr_fits=self.curr_embr_node.child(0)
			self.curr_pre_node=self.curr_embr_node.child(1)
			self.curr_noise_node=self.curr_embr_node.child(2)
		else:
			self.curr_fits=None
			self.curr_pre_node=None
			self.curr_noise_node=None
		self.curr_fit_node=None
		
		if self.curr_mol_node!=None:
			self.curr_embryos=self.curr_mol_node.child(0)
			self.curr_bkgds=self.curr_mol_node.child(1)	
		else:
			self.curr_embryos=None
			self.curr_bkgds=None
			
		if self.curr_bkgd_node!=None:
			self.curr_bkgd_pre_node=self.curr_bkgd_node.child(0)
		else:
			self.curr_bkgd_pre_node=None
		
		#Find out which embryo,molecule,bkgd is selected
		for mol in self.molecules:
			if self.curr_mol_node.data(0,0).toString()==mol.name:
				self.curr_mol=mol
		if self.curr_embr_node!=None:
			for emb in self.curr_mol.embryos:
				if self.curr_embr_node.data(0,0).toString()==emb.name:
					self.curr_embr=emb
		if self.curr_bkgd_node!=None:
			for bkgd in self.curr_mol.bkgds:
				if self.curr_bkgd_node.data(0,0).toString()==bkgd.name:
					self.curr_bkgd=bkgd
		
		#If selected Fits or Bkgds or Embryos, just do nothing
		if self.curr_node.data(0,0)=="Fits" or self.curr_node.data(0,0)=="Backgrounds" or self.curr_node.data(0,0)=="Embryos":
			return
		
		#Check if selected node is molecule
		if self.curr_mol_node==self.curr_node:
			
			self.parent_node=None
			self.curr_fit=None
			self.curr_bkgd=None
			self.curr_obj=self.curr_mol
			
			#Adding all properties of curr_embr to prop_list	
			for item in vars(self.curr_mol):
	
				#Don't print out large arrays
				if isinstance(vars(self.curr_mol)[str(item)],(int,float,str)) or vars(self.curr_mol)[str(item)]==None:
					curr_str=item+"="+str(vars(self.curr_mol)[str(item)])
					self.prop_list.addItem(curr_str)
		
		#Check if selected node is embryo
		elif self.curr_embr_node==self.curr_node:
			
			self.parent_node=self.curr_mol_node
			self.curr_fit=None
			self.curr_bkgd=None
			self.curr_obj=self.curr_embr
			
			#Adding all properties of curr_embr to prop_list	
			for item in vars(self.curr_embr):
	
				#Don't print out large arrays
				if isinstance(vars(self.curr_embr)[str(item)],(int,float,str)) or vars(self.curr_embr)[str(item)]==None:
					curr_str=item+"="+str(vars(self.curr_embr)[str(item)])
					self.prop_list.addItem(curr_str)
		else:
	
			self.parent_node=self.curr_node.parent()
							
			#Get fit index if fit is selected
			if self.parent_node==self.curr_fits:
				
				fit_ind=self.parent_node.indexOfChild(self.curr_node)
				self.curr_fit=self.curr_embr.fits[fit_ind]
			
				self.curr_fit_node=self.curr_node
				self.curr_obj=self.curr_fit
				
				#Adding all properties of curr_fit to prop_list	
				for item in vars(self.curr_fit):
		
					#Don't print out large arrays
					if isinstance(vars(self.curr_fit)[str(item)],(int,float,str)) or vars(self.curr_fit)[str(item)]==None:
						curr_str=item+"="+str(vars(self.curr_fit)[str(item)])
						self.prop_list.addItem(curr_str)
						
			#Pre selected
			if self.parent_node==self.curr_embr_node and self.curr_node.data(0,0).toString()=="Pre":
				
				self.curr_pre=self.curr_embr.pre
				self.curr_pre_node=self.curr_node
				self.curr_obj=self.curr_pre
				
				#Adding all properties of curr_fit to prop_list	
				for item in vars(self.curr_pre):
		
					#Don't print out large arrays
					if isinstance(vars(self.curr_pre)[str(item)],(int,float,str)) or vars(self.curr_pre)[str(item)]==None:
						curr_str=item+"="+str(vars(self.curr_pre)[str(item)])
						self.prop_list.addItem(curr_str)	
						
			#Bkgd Pre selected
			if self.parent_node==self.curr_bkgd_node and self.curr_node.data(0,0).toString()=="Pre":
				
				self.curr_bkgd_pre=self.curr_bkgd.pre
				self.curr_bkgd_pre_node=self.curr_node
				self.curr_obj=self.curr_bkgd_pre
				
				#Adding all properties of curr_fit to prop_list	
				for item in vars(self.curr_bkgd_pre):
		
					#Don't print out large arrays
					if isinstance(vars(self.curr_bkgd_pre)[str(item)],(int,float,str)) or vars(self.curr_bkgd_pre)[str(item)]==None:
						curr_str=item+"="+str(vars(self.curr_bkgd_pre)[str(item)])
						self.prop_list.addItem(curr_str)						
			
			#Noise selected
			if self.parent_node==self.curr_embr_node and self.curr_node.data(0,0).toString()=="Noise":
				
				self.curr_noise=self.curr_embr.noise
				self.curr_noise_node=self.curr_node
				self.curr_obj=self.curr_noise
				
				#Adding all properties of curr_fit to prop_list	
				for item in vars(self.curr_noise):
		
					#Don't print out large arrays
					if isinstance(vars(self.curr_noise)[str(item)],(int,float,str)) or vars(self.curr_noise)[str(item)]==None:
						curr_str=item+"="+str(vars(self.curr_noise)[str(item)])
						self.prop_list.addItem(curr_str)			
						
			#Get bkgd index if bkgd is selected
			if self.parent_node==self.curr_bkgds:
		
				self.bkgd_ind=self.parent_node.indexOfChild(self.curr_node)
				self.curr_bkgd=self.curr_mol.bkgds[self.bkgd_ind]
				self.curr_bkgd_node=self.curr_node
				self.curr_obj=self.curr_bkgd
				
				#Adding all properties of curr_fit to prop_list	
				for item in vars(self.curr_bkgd):
		
					#Don't print out large arrays
					if isinstance(vars(self.curr_bkgd)[str(item)],(int,float,str)) or vars(self.curr_bkgd)[str(item)]==None:
						curr_str=item+"="+str(vars(self.curr_bkgd)[str(item)])
						self.prop_list.addItem(curr_str)		
						
		#Sort prop_list
		self.prop_list.sortItems()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Print out hierachy tree of embryo_list
	
	def print_embryo_list(self):
		
		for i in range(self.embryos_list.topLevelItemCount()):
				print "-Top level item:", self.embryos_list.topLevelItem(i)," with name :", self.embryos_list.topLevelItem(i).data(0,0).toString(),"and parent:",self.embryos_list.topLevelItem(i).parent()
				
				for j in range(self.embryos_list.topLevelItem(i).childCount()):
					print "----" , self.embryos_list.topLevelItem(i).child(j)," with name :", self.embryos_list.topLevelItem(i).child(j).data(0,0).toString(), "and parent:",self.embryos_list.topLevelItem(i).child(j).parent()
					
					for k in range(self.embryos_list.topLevelItem(i).child(j).childCount()):
						print "   --- " , self.embryos_list.topLevelItem(i).child(j).child(k)," with name :", self.embryos_list.topLevelItem(i).child(j).child(k).data(0,0).toString(),"and parent:",self.embryos_list.topLevelItem(i).child(j).child(k).parent()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Delete embryo
		
	def delete_embryo(self):
		
		#Check if anything is highlighted
		if  self.curr_embr_node==None:
		
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return		
			
		#Remove from list of embryo objects
		self.curr_mol.embryos.remove(self.curr_embr)
		
		#Remove from sidebar
		ind=self.curr_embryos.indexOfChild(self.curr_embr_node)
		self.curr_embryos.takeChild(ind)
		
		self.curr_node=self.embryos_list.currentItem()
		
		if self.curr_node!=None:
			self.show_embryo_props()
		else:
			self.prop_list.clear()
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Delete molecule
		
	def delete_molecule(self):
		
		#Check if anything is highlighted
		if self.embryos_list.currentItem()==None:
		
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
			
		#Remove from list of embryo objects
		self.molecules.remove(self.curr_mol)
		
		#Remove from sidebar
		ind=self.embryos_list.indexOfTopLevelItem(self.curr_mol_node)
		self.embryos_list.takeTopLevelItem(ind)
		
		self.curr_node=self.embryos_list.currentItem()
		
		if self.curr_node!=None:
			self.show_embryo_props()
		else:
			self.prop_list.clear()		
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Delete specific molecule
	
	def delete_specific_molecule(self,molecule=None):
		
		#Find out which node molecule is in embryo column
		for i in range(self.embryos_list.topLevelItemCount()):
			if self.embryos_list.topLevelItem(i).data(0,0).toString()==molecule.name:
				self.embryos_list.takeTopLevelItem(i)
				break
		
		#Remove from list of embryo objects
		self.molecules.remove(self.curr_mol)
		
		self.curr_node=self.embryos_list.currentItem()
		if self.curr_node!=None:
			self.show_embryo_props()
		else:
			self.prop_list.clear()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Delete specific embryo
	
	def delete_specific_embryo(self,embryo=None):
		
		#Find out which node embryo is in embryos column
		for i in range(self.curr_embryos.childCount()):
			if self.curr_embryos.child(i).data(0,0).toString()==embryo.name:
				self.curr_embryos.takeChild(i)
				break
		
		#Remove from list of embryo objects
		self.curr_mol.embryos.remove(self.curr_embr)
		
		self.curr_node=self.embryos_list.currentItem()
		if self.curr_node!=None:
			self.show_embryo_props()
		else:
			self.prop_list.clear()
					
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Add fit
	
	def add_fit(self):
		
		#Check if top level node is a embryo
		if self.curr_embr_node==None:
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			curr_name="fit_0"
			
			#Check if there is already a fit with the name fit_0
			in_list=1
			
			while in_list==1:
				#Go through all embryos
				for i in range(self.curr_fits.childCount()):
					
					if self.curr_fits.child(i).data(0,0).toString()==curr_name:
						#If exists, split name in nmbr and name and add +1 to number
						nme,nmbr=curr_name.split("_")
						nmbr=str(int(nmbr)+1)
						curr_name=nme+"_"+nmbr
				
				#Actually in list
				in_list=0
			
			#Add to embryo bar
			self.curr_fit_node=QtGui.QTreeWidgetItem(self.curr_fits,[curr_name,'','0'])
			
			#Create fit object
			self.curr_embr.add_fit(shape(self.curr_embr.fits)[0],curr_name,"default")
			self.curr_fit=self.curr_embr.fits[-1]
			
			#If molecule has been analyzed, apply optimal values
			if shape(self.curr_fit.embryo.ext_av_data_d)[0]>0 and shape(self.curr_mol.bkgds)[0]>0:
				if shape(self.curr_mol.bkgds[0].bkgd_ext_vec)[0]>0:
					self.curr_fit=pyfdap_fit.opt_parms_fit(self.curr_fit,self.curr_mol)
			
			ret=pyfdap_subwin.fit_dialog(self.curr_fit,self.curr_mol,self.curr_embr,self).exec_()
			
			#Update fit name in embryos column
			self.curr_fit_node.setText(0,self.curr_fit.name)
			self.curr_node=self.curr_fit_node
			
			self.embryos_list.setCurrentItem(self.curr_fit_node)
			
			#Perform the new fit
			self.perform_fit()
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Copy fit
	
	def copy_fit(self):
		
		if self.curr_node.parent()==self.curr_fits:
			
			curr_name=self.curr_fit.name
			curr_fit_mumber=self.curr_fit.fit_number
			self.curr_fit=cpy.deepcopy(self.curr_fit)
			
			
			if "copy" in curr_name:
				new_name=curr_name
			else:
				new_name=curr_name+"_copy1"
			
			in_list=1
			while in_list==1:
				for fit in self.curr_embr.fits:
					
					if fit.name==new_name:
						
						empty,nmbr=new_name.split("_copy")
						nmbr=str(int(nmbr)+1)
						new_name=empty+"_copy"+nmbr
						break
				in_list=0
			
			self.curr_fit.name=new_name
			
			#Get maximum fit number
			tmp=[]
			for fit in self.curr_embr.fits:
				tmp.append(fit.fit_number)
			
			self.curr_fit.fit_number=max(tmp)+1
			self.curr_embr.fit_number=self.curr_embr.fit_number+1
			self.curr_embr.fits.append(self.curr_fit)
			self.curr_fit_node=QtGui.QTreeWidgetItem(self.curr_fits,[new_name,'','0'])		
	
	def copy_fit_for_other_embryo(self,fit,embryo):
		
		new_fit=cpy.deepcopy(fit)
		new_fit.embryo=embryo
		
		for fit in embryo.fits:
			if fit.name==new_fit.name:
				new_fit.name=new_fit.name+"_copied_from"+embryo.name
		
		embryo.fit_number=embryo.fit_number+1
		new_fit.fit_number=embryo.fit_number-1
		embryo.fits.append(new_fit)
		
		for i in range(self.curr_embryos.childCount()):
		
			if self.curr_embryos.child(i).data(0,0).toString()==embryo.name:
				curr_embr_node=self.curr_embryos.child(i)
				curr_fits=curr_embr_node.child(0)
				QtGui.QTreeWidgetItem(curr_fits,[new_fit.name,'','0'])
		
	def copy_fit_to_all(self):
		
		if self.curr_fit==None:
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		for emb in self.curr_mol.embryos:
			
			if emb!=self.curr_embr:
			
				self.copy_fit_for_other_embryo(self.curr_fit,emb)
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Delete fit
	
	def delete_fit(self):
		
		#Check what is highlighted
		if self.embryos_list.currentItem()==None or self.embryos_list.currentItem()==self.curr_embr_node or self.embryos_list.currentItem().parent()==self.curr_bkgds or self.embryos_list.currentItem().parent()==self.curr_embr_node:
		
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
									
		#Remove from list of embryo objects		
		self.curr_embr.fits.remove(self.curr_fit)
		self.curr_embr.fit_number=self.curr_embr.fit_number-1
		
		#Remove from embryo bar
		self.parent_node.removeChild(self.curr_node)
		self.prop_list.clear()
		
		self.curr_node=self.embryos_list.currentItem()
		if self.curr_node!=None:
			self.show_embryo_props()
		else:
			self.prop_list.clear()
					
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit fit
	
	def edit_fit(self):
		
		#Check if highlighted node is child of embryo
		if self.embryos_list.currentItem().parent()!=self.curr_fits:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Open fit dialog
			ret=pyfdap_subwin.fit_dialog(self.curr_fit,self.curr_mol,self.curr_embr,self).exec_()
			
			#Update property column
			self.show_embryo_props()
			
			#Update fit name in embryos column
			self.curr_fit_node=self.embryos_list.currentItem()
			self.curr_fit_node.setText(0,self.curr_fit.name)
				
	def edit_mult_fit(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_mol==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Backup saving sel_fits
			sel_fits_backup=list(self.curr_mol.sel_fits)
			
			#Open fit selection dialog
			ret=pyfdap_subwin.select_fits(self.curr_mol,0,self).exec_()
			
			if self.curr_mol.sel_fits!=[]:
				#Open mult fit editing dialog
				ret=pyfdap_subwin.mult_fit_dialog(self.curr_mol,self.curr_mol.sel_fits,self).exec_()
				
			#Mapping back sel_fits_backup
			self.curr_mol.sel_fits=list(sel_fits_backup)
			
			#Update property column
			self.show_embryo_props()
				
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Add bkgd
	
	def add_bkgd(self):
		
		#Check if highlighted node is a embryo
		if self.curr_mol_node==None:
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			curr_name="bkgd_0"
			
			#Check if there is already a bkgd with the name bkgd_0
			in_list=1
			
			while in_list==1:
				#Go through all embryos
				for i in range(self.curr_bkgds.childCount()):
					
					if self.curr_bkgds.child(i).data(0,0).toString()==curr_name:
						#If exists, split name in nmbr and name and add +1 to number
						nme,nmbr=curr_name.split("_")
						nmbr=str(int(nmbr)+1)
						curr_name=nme+"_"+nmbr
				
				#Actually in list
				in_list=0
			
			#Add to embryo bar
			self.curr_bkgd_node=QtGui.QTreeWidgetItem(self.curr_bkgds,[curr_name,'0',''])
			self.curr_bkgd_pre_node=QtGui.QTreeWidgetItem(self.curr_bkgd_node,["Pre",'0',''])
			
			#Create new bkgd object
			self.curr_mol.add_bkgd(shape(self.curr_mol.bkgds)[0],curr_name,"default")
			self.curr_bkgd=self.curr_mol.bkgds[-1]
			
			#Show new bkgd props in props column
			self.embryos_list.setCurrentItem(self.curr_bkgd_node)
			self.show_embryo_props()
			
			#Open Bkgd dialog for bkgd dataset
			ret=pyfdap_subwin.bkgd_dataset_dialog(self.curr_bkgd,self).exec_()
			
			self.edit_bkgd_pre()
			
			#Update property column
			self.show_embryo_props()
			
			#Update bkgd name in embryos column
			self.curr_bkgd_node.setText(0,self.curr_bkgd.name)
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Delete bkgd
	
	def delete_bkgd(self):
			
		#Check what is highlighted
		if self.curr_node.parent()!=self.curr_bkgds:
		
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
							
		#Remove from list of embryo objects		
		self.curr_mol.bkgds.remove(self.curr_bkgd)
		
		#Remove from embryo bar
		self.curr_bkgds.removeChild(self.curr_node)
		self.prop_list.clear()
		
		self.curr_node=self.embryos_list.currentItem()
		if self.curr_node!=None:
			self.show_embryo_props()
		else:
			self.prop_list.clear()		
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit bkgd
	
	def edit_bkgd(self):
		
		#Check if highlighted node is child of embryo
		if self.embryos_list.currentItem().parent()!=self.curr_bkgds:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Open Bkgd dialog for bkgd data set
			ret=pyfdap_subwin.bkgd_dataset_dialog(self.curr_bkgd,self).exec_()
			
			#Update property column
			self.show_embryo_props()
			
			#Update bkgd name in embryos column
			self.curr_bkgd_node=self.embryos_list.currentItem()
			self.curr_bkgd_node.setText(0,self.curr_bkgd.name)
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit pre dataset
	
	def edit_bkgd_pre(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_bkgd_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Open Bkgd dialog for bkgd dataset
			ret=pyfdap_subwin.pre_dataset_dialog(self.curr_bkgd.pre,self).exec_()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit main dataset
	
	def edit_dataset(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_embr_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Open Bkgd dialog for bkgd dataset
			ret=pyfdap_subwin.dataset_dialog(self.curr_embr,self).exec_()
			
			#Update property column
			self.embryos_list.setCurrentItem(self.curr_embr_node)
			self.show_embryo_props()
			
			#Update embr name in embryos column
			self.curr_embr_node=self.embryos_list.currentItem()
			self.curr_embr_node.setText(0,self.curr_embr.name)
	
	def select_ignored_frames(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_embr_node==None and self.curr_bkgd_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo or background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Open selection tool
		if self.curr_embr_node!=None:
			ret=pyfdap_subwin.select_ignored_frames(self.curr_embr,self).exec_()
		elif self.curr_bkgd_node!=None:
			ret=pyfdap_subwin.select_ignored_frames(self.curr_bkgd,self).exec_()
	
	def edit_thresh(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_embr_node==None and self.curr_bkgd_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo or background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Open selection tool
		if self.curr_embr_node!=None:
			ret=pyfdap_subwin.select_threshhold(self.curr_embr,self).exec_()
		elif self.curr_bkgd_node!=None:
			ret=pyfdap_subwin.select_threshhold(self.curr_bkgd,self).exec_()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit Molecule
	
	def edit_molecule(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_mol_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			#Open Bkgd dialog for bkgd data set
			ret=pyfdap_subwin.molecule_dialog(self.curr_mol,self).exec_()
			
			#Update mol name in embryos column
			self.curr_mol_node.setText(0,self.curr_mol.name)		
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit pre data set
	
	def edit_pre(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_embr_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Open Bkgd dialog for bkgd data set
			ret=pyfdap_subwin.pre_dataset_dialog(self.curr_embr.pre,self).exec_()
			
	
	#----------------------------------------------------------------------------------------------------------------------------------------
        #Edit noise data set
	
	def edit_noise(self):
		
		#Check if highlighted node is child of embryo
		if self.curr_embr_node==None:
			
			#If not give error msg
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		else:
			
			#Open Bkgd dialog for bkgd data set
			ret=pyfdap_subwin.noise_dataset_dialog(self.curr_embr.noise,self).exec_()		
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Copy any
	
	def copy_any(self):
	
		if self.curr_embr_node!=None:
			self.copy_embryo()
		elif self.curr_bkgd_node!=None:
			self.copy_bkgd()
		elif self.curr_embr_node==None and self.curr_bkgd_node==None:
			self.copy_molecule()
		
	def paste_any(self):	
		
		if self.copied_embr!=None:
			self.paste_embryo()
		elif self.copied_bkgd!=None:
			self.paste_bkgd()
		else:
			QtGui.QMessageBox.critical(None, "Error","Nothing in clipboard.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return	
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Copy embryo
	
	def copy_embryo(self):
		
		#Check if an embryo is selected
		if self.curr_embr_node!=None:
			self.copied_embr=cpy.deepcopy(self.curr_embr)
			self.copied_bkgd=None
		else:
			QtGui.QMessageBox.critical(None, "Error","No background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return	
		
	def paste_embryo(self):
		if self.curr_mol!=None:
			if self.copied_embr!=None:
				emb=self.copied_embr
				
				name=str(emb.name)
				
				for embr in self.curr_mol.embryos:
					if name==str(embr.name):
						curr_name=name+"_copy"
					else:
						curr_name=name
						
				#Add to embryo bar
				analyzed=str(0)
				fitted=str(0)
				if shape(emb.ext_av_data_d)[0]>1:
					analyzed=str(1)
				if shape(emb.fits)>0:
					for fit in emb.fits:
						if shape(fit.fit_av_d)[0]>1:	
							fitted=str(1)
				
				#Adding embryo to sidebar
				self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[emb.name,analyzed,fitted])
				self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'',''])
				self.curr_pre_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Pre",analyzed,''])
				
				#Add fits if they exist
				for fit in emb.fits:
					if shape(fit.fit_av_d)[0]>0:		
						QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','1'])
					else:	
						QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','0'])
									
				self.curr_embr=emb
				self.curr_mol.embryos.append(emb)
				self.embryos_list.expandItem(self.curr_embr_node)
				self.embryos_list.expandItem(self.curr_fits)
				self.copied_embr=None
				self.copied_bkgd=None
				
			else:
				QtGui.QMessageBox.critical(None, "Error","No embryo copied yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				return
		
		else:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Copy bkgd
	
	def copy_bkgd(self):
		
		#Check if an embryo is selected
		if self.curr_bkgd_node!=None:	
			self.copied_bkgd=cpy.deepcopy(self.curr_bkgd)
			self.copied_embr=None
		else:
			QtGui.QMessageBox.critical(None, "Error","No background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		
	def paste_bkgd(self):
		if self.curr_mol!=None:
			if self.copied_bkgd!=None:
				bkg=self.copied_bkgd
				
				name=str(bkg.name)
				
				for b in self.curr_mol.bkgds:
					if name==str(b.name):
						curr_name=name+"_copy"
					else:
						curr_name=name
						
				#Add to embryo bar
				analyzed=str(0)
				if shape(bkg.bkgd_ext_vec)[0]>1:
					analyzed=str(1)
					
				#Adding bkgd to sidebar
				self.curr_bkgd_node=QtGui.QTreeWidgetItem(self.curr_bkgds,[bkg.name,analyzed,''])
				self.curr_bkgd_pre_node=QtGui.QTreeWidgetItem(self.curr_bkgd_node,["Pre",analyzed,''])
				
				#Add to molecule
				self.curr_mol.bkgds.append(bkg)
				
				self.curr_bkgd=bkg
				self.embryos_list.expandItem(self.curr_bkgd_node)	
				self.copied_embr=None
				self.copied_bkgd=None
				
			else:
				QtGui.QMessageBox.critical(None, "Error","No embryo copied yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				return
		else:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Copy molecule
	
	def copy_molecule(self):
		
		#Check if a molecule is selected
		if self.curr_mol_node!=None:
			
			name=str(self.curr_mol.name)
			curr_name=name+"_copy"
			
			self.curr_mol=cpy.deepcopy(self.curr_mol)
			self.curr_mol.name=curr_name
			self.molecules.append(self.curr_mol)
			
			#Adding molecule to sidebar
			self.curr_node=QtGui.QTreeWidgetItem(self.embryos_list,[curr_name,"",""])
			self.curr_embryos=QtGui.QTreeWidgetItem(self.curr_node,["Embryos","",""])
			self.curr_bkgds=QtGui.QTreeWidgetItem(self.curr_node,["Backgrounds","",""])
			self.curr_mol_node=self.curr_node
			
			for emb in self.curr_mol.embryos:
		
				#Add to embryo bar
				analyzed=str(0)
				fitted=str(0)
				if shape(emb.ext_av_data_d)[0]>1:
					analyzed=str(1)
				if shape(emb.fits)>0:
					for fit in emb.fits:
						if shape(fit.fit_av_d)[0]>1:	
							fitted=str(1)
					
				#Adding embryo to sidebar
				self.curr_embr_node=QtGui.QTreeWidgetItem(self.curr_embryos,[emb.name,analyzed,fitted])
				self.curr_fits=QtGui.QTreeWidgetItem(self.curr_embr_node,["Fits",'',''])
				self.curr_pre_node=QtGui.QTreeWidgetItem(self.curr_embr_node,["Pre",analyzed,''])
				
				#Add fits if they exist
				for fit in emb.fits:
					if shape(fit.fit_av_d)[0]>0:		
						QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','1'])
					else:	
						QtGui.QTreeWidgetItem(self.curr_fits,[fit.name,'','0'])
				
				
				self.curr_embr=emb
				self.embryos_list.expandItem(self.curr_fits)
				
			#Add bkgds if they exist
			for bkgd in self.curr_mol.bkgds:	
				if shape(bkgd.bkgd_vals_slice)[0]>0:
					QtGui.QTreeWidgetItem(self.curr_bkgds,[bkgd.name,'1',''])
				else:
					QtGui.QTreeWidgetItem(self.curr_bkgds,[bkgd.name,'0',''])
								
			self.embryos_list.expandItem(self.curr_embryos)
			self.embryos_list.expandItem(self.curr_bkgds)
			
		else:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Edit properties via doubleclick on property
	
	def edit_prop(self):
		
		#Getting name of property
		varname,varvalue=self.prop_list.currentItem().text().split("=")
		
		#Getting property
		for item in vars(self.curr_embr):
			if varname==item:
				
				#For general purpose, define what the current property is
				self.curr_prop=vars(self.curr_embr)[str(item)]
				
				#Check what dialogbox to use
				if isinstance(self.curr_prop,(int)):
					varvalue, ok=QtGui.QInputDialog.getInt(self, "Set"+varname, varname+"=")
					varvalue=int(varvalue)
				elif isinstance(self.curr_prop,(float)):
					varvalue, ok=QtGui.QInputDialog.getDouble(self, "Set"+varname, varname+"=")
					varvalue=float(varvalue)
				elif isinstance(self.curr_prop,(str)):
					varvalue, ok=QtGui.QInputDialog.getText(self, "Set"+varname, varname+"=")
					varvalue=str(varvalue)
					
				#Enter the new value into curr_embr
				vars(self.curr_embr)[str(item)]=varvalue
				break
		
		if varname=="name":
			self.curr_node.setText(0,varvalue)
			
		#Update property bar
		self.show_embryo_props()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Analyze complete molecule
	
	def analyze_all(self):
		
		if self.curr_mol_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		#Make backup copy of molecule
		if self.curr_conf.backup_to_file:
			self.fn_backup=self.lastopen+"/"+self.curr_mol.name+"_backup.pk"
			self.curr_mol.save_molecule(self.fn_backup)
		if self.curr_conf.backup_to_mem:
			self.temp_mol=cpy.deepcopy(self.curr_mol)
		
		#Generate wait popup
		self.wait_popup=pyfdap_subwin.analyze_all_prog(None)
		self.wait_popup.accepted.connect(self.analyze_all_canceled)
		self.statusBar().showMessage("Analyzing molecule " + self.curr_mol.name)
		self.setDisabled(True)
		
		#Generate Qthread and pass analyze there
		self.analyze_all_task=pyfdap_subwin.analyze_all_thread(molecule=self.curr_mol)
		self.analyze_all_task.taskFinished.connect(self.analyze_all_finished)
		self.analyze_all_task.start()
				
	def analyze_all_finished(self):
		
		self.wait_popup.close()
		self.statusBar().showMessage("Idle")
		#Setting analyzed=1
		for i in range(self.curr_embryos.childCount()):
			self.curr_embryos.child(i).setText(1,"1")
			for j in range(self.curr_embryos.child(i).childCount()):
				if self.curr_embryos.child(i).child(j).data(0,0).toString()!="Fits":
					self.curr_embryos.child(i).child(j).setText(1,"1")
			
		for i in range(self.curr_bkgds.childCount()):
			self.curr_bkgds.child(i).setText(1,"1")
			for j in range(self.curr_bkgds.child(i).childCount()):
				self.curr_bkgds.child(i).child(j).setText(1,"1")
		
		self.setEnabled(True)	
		
		if self.curr_conf.backup_to_mem:
			self.temp_mol=None
		if self.curr_conf.backup_to_file:
			os.remove(self.fn_backup)
			
		return
	
	def analyze_all_canceled(self):
		self.statusBar().showMessage("Idle")
		self.setEnabled(True)
		
		self.analyze_all_task.terminate()
		if self.curr_conf.backup_to_file:
			self.curr_mol=self.curr_mol.load_molecule(self.fn_backup)
			os.remove(self.fn_backup)
		if self.curr_conf.backup_to_mem:
			self.curr_mol=cpy.deepcopy(self.temp_mol)
			self.temp_mol=None
			
		self.wait_popup.close()	
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Analyze bkgd data sets
	
	def analyze_bkgds(self):
		
		if self.curr_mol_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return

		#Make backup copy of molecule
		if self.curr_conf.backup_to_file:
			self.fn_backup=self.lastopen+"/"+self.curr_mol.name+"_backup.pk"
			self.curr_mol.save_molecule(self.fn_backup)
		if self.curr_conf.backup_to_mem:
			self.temp_mol=cpy.deepcopy(self.curr_mol)
		
		#Generate wait popup
		self.wait_popup=pyfdap_subwin.analyze_bkgd_prog(None)
		self.wait_popup.accepted.connect(self.analyze_bkgds_canceled)
		self.statusBar().showMessage("Analyzing Dataset")
		self.setDisabled(True)
		
		#Generate Qthread and pass analyze there
		self.analyze_bkgd_task=pyfdap_subwin.analyze_bkgd_thread(molecule=self.curr_mol)
		self.analyze_bkgd_task.taskFinished.connect(self.analyze_bkgds_finished)
		self.analyze_bkgd_task.start()
				
	def analyze_bkgds_finished(self):
		
		self.wait_popup.close()
		self.statusBar().showMessage("Idle")
		#Setting analyzed=1
		for i in range(self.curr_bkgds.childCount()):
			self.curr_bkgds.child(i).setText(1,"1")
		
		self.setEnabled(True)	
		if self.curr_conf.backup_to_mem:
			self.temp_mol=None
		if self.curr_conf.backup_to_file:
			os.remove(self.fn_backup)
		return
	
	def analyze_bkgds_canceled(self):
		self.statusBar().showMessage("Idle")
		self.setEnabled(True)
		
		self.analyze_bkgd_task.terminate()
		if self.curr_conf.backup_to_file:
			self.curr_mol=self.curr_mol.load_molecule(self.fn_backup)
			os.remove(self.fn_backup)
		if self.curr_conf.backup_to_mem:
			self.curr_mol=cpy.deepcopy(self.temp_mol)
			self.temp_mol=None
		
		self.wait_popup.close()	
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Analyze data set
	
	def analyze_dataset(self):
		
		if self.curr_embr_node==None:
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Make backup copy of embryo
		self.backup_emb=cpy.deepcopy(self.curr_embr)
		
		#Generate wait popup
		self.wait_popup=pyfdap_subwin.analyze_prog(None)
		self.wait_popup.accepted.connect(self.analyze_canceled)
		self.statusBar().showMessage("Analyzing Dataset")
		self.setDisabled(True)
		
		#Generate Qthread and pass analyze there
		self.analyze_task=pyfdap_subwin.analyze_thread(embryo=self.curr_embr)
		self.analyze_task.taskFinished.connect(self.analyze_finished)
		self.analyze_task.start()
				
	def analyze_finished(self):
		
		self.wait_popup.close()
		self.statusBar().showMessage("Idle")
		#Setting analyzed=1
		self.curr_embr_node.setText(1,"1")
		self.curr_pre_node.setText(1,"1")
		self.curr_noise_node.setText(1,"1")
		
		self.setEnabled(True)
			
		return
	
	def analyze_canceled(self):
		self.statusBar().showMessage("Idle")
		self.setEnabled(True)
		
		self.analyze_task.terminate()
		self.curr_embr=cpy.deepcopy(self.backup_emb)
		self.backup_emb=None
		self.wait_popup.close()
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Plot data 
	
	def plot_data_timeseries(self):
		
		if self.curr_embr_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if shape(self.curr_embr.ext_av_data_d)[0]>0:
			pass
		else:
			QtGui.QMessageBox.critical(None, "Error","No analyzed data.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		self.create_plot_tab("data")
		if len(self.curr_embr.ignored)>0:
			self.try_plot_data(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b','--','int_data without ign',self.ax)
			self.try_plot_data(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r','--','ext_data without ign',self.ax)
			self.try_plot_data(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g','--','slice_data without ign',self.ax)
		else:	
			self.try_plot_data(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b','-','int_data',self.ax)
			self.try_plot_data(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r','-','ext_data',self.ax)
			self.try_plot_data(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g','-','slice_data',self.ax)
			
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		
		self.adjust_canvas()
	
	def try_plot_data(self,tvec,data,color,ls,label,ax):
		
		try:
			ax.plot(tvec,data,c=color,ls=ls,label=label)
		except:
			print "Was not able to plot ", label
			
		return ax
		
	
	def plot_data_bkgd(self):
		
		if self.curr_embr_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if shape(self.curr_embr.ext_av_data_d)[0]>0:
			pass
		else:
			QtGui.QMessageBox.critical(None, "Error","No analyzed data.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if self.curr_mol.bkgd_slice==None:
			QtGui.QMessageBox.critical(None, "Error","No analyzed background.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
			
		self.create_plot_tab("dataandbkgd")
		
		self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b-',label="int_data")
		self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r-',label="ext_data")
		self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-',label="slice_data")
		
		bkgd_slice_vec=ones(shape(self.curr_embr.ext_av_data_d))*self.curr_mol.bkgd_slice
		bkgd_ext_vec=ones(shape(self.curr_embr.ext_av_data_d))*self.curr_mol.bkgd_ext
		bkgd_int_vec=ones(shape(self.curr_embr.ext_av_data_d))*self.curr_mol.bkgd_int
		
		self.ax.plot(self.curr_embr.tvec_data,bkgd_int_vec,'b--',label="int_bkgd")
		self.ax.plot(self.curr_embr.tvec_data,bkgd_ext_vec,'r--',label="ext_bkgd")
		self.ax.plot(self.curr_embr.tvec_data,bkgd_slice_vec,'g--',label="slice_bkgd")
		
		pre_slice_vec=ones(shape(self.curr_embr.ext_av_data_d))*self.curr_embr.pre.pre_slice
		pre_ext_vec=ones(shape(self.curr_embr.ext_av_data_d))*self.curr_embr.pre.pre_ext
		pre_int_vec=ones(shape(self.curr_embr.ext_av_data_d))*self.curr_embr.pre.pre_int
		
		self.ax.plot(self.curr_embr.tvec_data,pre_int_vec,'b-.',label="int_pre")
		self.ax.plot(self.curr_embr.tvec_data,pre_ext_vec,'r-.',label="ext_pre")
		self.ax.plot(self.curr_embr.tvec_data,pre_slice_vec,'g-.',label="slice_pre")
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()		
	
	def plot_embryo_slice_imgs(self):
		self.plot_embryo_img_series("slice")
		
	def plot_embryo_ext_imgs(self):
		self.plot_embryo_img_series("ext")	
	
	def plot_embryo_int_imgs(self):
		self.plot_embryo_img_series("int")
	
	def plot_embryo_masks_embryo(self):
		self.plot_embryo_img_series("masks_embryo")
		
	def plot_embryo_masks_ext(self):
		self.plot_embryo_img_series("masks_ext")	
	
	def plot_embryo_masks_int(self):
		self.plot_embryo_img_series("masks_int")
	
	def plot_embryo_img_series(self,imgtype):
		
		if self.curr_embr_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Check if embryo is already analyzed
		if shape(self.curr_embr.slice_av_data_d)[0]==0:
			reply = QtGui.QMessageBox.question(self, 'Message',"Data set has not been analyzed yet?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
			if reply == QtGui.QMessageBox.Yes:
				self.analyze_dataset()
				return
			else:
				return
		elif shape(self.curr_embr.masks_embryo)[0]==0 and shape(self.curr_embr.slice_av_data_d)[0]>0:
			reply = QtGui.QMessageBox.question(self, 'Message',"Image data is currently not loaded, do you want to analyze the data again?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
			if reply == QtGui.QMessageBox.Yes:
				self.analyze_dataset()
				return
			else:
				return
			
		self.create_slider_plot_tab(imgtype)
		self.img_axes=[]
		
		self.curr_tab.imgs=[]
		
		for i in range(shape(self.curr_embr.vals_slice)[0]):
			if imgtype=="masks_embryo":
				img=self.curr_embr.masks_embryo[i]
			elif imgtype=="masks_ext":
				img=self.curr_embr.masks_ext[i]
			elif imgtype=="masks_int":
				img=self.curr_embr.masks_int[i]
			elif imgtype=="slice":
				img=self.curr_embr.vals_slice[i]
			elif imgtype=="ext":
				img=self.curr_embr.vals_slice[i]*self.curr_embr.masks_ext[i]
			elif imgtype=="int":
				img=self.curr_embr.vals_slice[i]*self.curr_embr.masks_int[i]
				
			self.img_axes.append(self.ax.imshow(img,cmap='jet'))
			self.curr_tab.imgs.append(img)
			
		self.curr_tab.img_axes=self.img_axes
		
		self.adjust_canvas()
		self.ax.draw_artist(self.curr_tab.img_axes[0])
	
	def plot_bkgd_slice_imgs(self):
		self.plot_bkgd_img_series("bkgd_slice")
		
	def plot_bkgd_ext_imgs(self):
		self.plot_bkgd_img_series("bkgd_ext")	
	
	def plot_bkgd_int_imgs(self):
		self.plot_bkgd_img_series("bkgd_int")
	
	def plot_bkgd_masks_embryo(self):
		self.plot_bkgd_img_series("bkgd_masks_embryo")
		
	def plot_bkgd_masks_ext(self):
		self.plot_bkgd_img_series("bkgd_masks_ext")	
	
	def plot_bkgd_masks_int(self):
		self.plot_bkgd_img_series("bkgd_masks_int")
		
	def plot_bkgd_img_series(self,imgtype):
		
		#Check if bkgd node is selected
		if self.embryos_list.currentItem().parent()!=self.curr_bkgds:
			
			QtGui.QMessageBox.critical(None, "Error","No background selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Check if background is already analyzed		
		if shape(self.curr_bkgd.bkgd_vals_slice)[0]==0 and shape(self.curr_embr.slice_av_data_d)[0]==0:
			
			reply = QtGui.QMessageBox.question(self, 'Message',"This background dataset has not been analyzed yet. Do you want to analyze the background datasets now?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

			if reply == QtGui.QMessageBox.Yes:
				self.analyze_bkgds()
				return
			else:
				return
			
		elif shape(self.curr_bkgd.bkgd_vals_slice)[0]==0 and shape(self.curr_embr.slice_av_data_d)[0]>0:
		
			#Do analysis for bkgd data set
			reply = QtGui.QMessageBox.question(self, 'Message',"This background dataset does not contain the image data anymore. Do you want to reanalyze all background datasets now?", QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

			if reply == QtGui.QMessageBox.Yes:
				self.analyze_bkgds()
				return
			else:
				return
			
		self.create_slider_plot_tab(imgtype)
		self.img_axes=[]
		self.curr_tab.imgs=[]
		for i in range(shape(self.curr_bkgd.bkgd_vals_slice)[0]):
			if imgtype=="bkgd_masks_embryo":
				img=self.curr_bkgd.masks_embryo[i]
			elif imgtype=="bkgd_masks_ext":
				img=self.curr_bkgd.masks_ext[i]
			elif imgtype=="bkgd_masks_int":
				img=self.curr_bkgd.masks_int[i]
			elif imgtype=="bkgd_slice":
				img=self.curr_bkgd.bkgd_vals_slice[i]
			elif imgtype=="bkgd_ext":
				img=self.curr_bkgd.bkgd_vals_slice[i]*self.curr_bkgd.masks_ext[i]
			elif imgtype=="bkgd_int":
				img=self.curr_bkgd.bkgd_vals_slice[i]*self.curr_bkgd.masks_int[i]
				
			self.img_axes.append(self.ax.imshow(img,cmap='jet'))
			self.curr_tab.imgs.append(img)
			
		self.curr_tab.img_axes=self.img_axes
		
		self.adjust_canvas()
		self.ax.draw_artist(self.curr_tab.img_axes[0])
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Plot fits
	
	def plot_fit(self):
		
		if self.curr_node.parent().data(0,0).toString()!="Fits":
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if shape(self.curr_fit.fit_av_d)[0]>0:
			pass
		else:
			QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		self.create_plot_tab("fit")
		
		if shape(self.curr_embr.ignored)[0]>0:
			if self.curr_fit.fit_ext==1:
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r*',label='ext_data')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.fit_av_d,'r--',label='fit')
			if self.curr_fit.fit_slice==1:
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g*',label='slice_data')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.fit_av_d,'g--',label='fit')
			if self.curr_fit.fit_int==1:
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b*',label='int_data')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.fit_av_d,'b--',label='fit')
		else:
			if self.curr_fit.fit_ext==1:
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r*',label='ext_data')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.fit_av_d,'r--',label='fit')
			if self.curr_fit.fit_slice==1:
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g*',label='slice_data')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.fit_av_d,'g--',label='fit')
			if self.curr_fit.fit_int==1:
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b*',label='int_data')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.fit_av_d,'b--',label='fit')
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		
		self.adjust_canvas()
		
	def plot_track_fit(self):
		
		if self.embryos_list.currentItem()!=self.curr_fit_node:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if shape(self.curr_fit.fit_av_d)[0]>0:
			pass
		else:
			QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if self.curr_fit.save_track==1:
			self.create_slider_plot_tab("track_fit")
			
			if shape(self.curr_embr.ignored)[0]>0:
									
				if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
					
					self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g*')
					self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[-1],'g--')
					
				elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
					self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r*')
					self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[-1],'r--')
					
				elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:

					self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b*')
					self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[-1],'b--')	
			
			else:
				if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
					
					self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g*')
					self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[-1],'g--')
					
				elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
					self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r*')
					self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[-1],'r--')
					
				elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:

					self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b*')
					self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[-1],'b--')	
		
			self.curr_tab.lbl_track_k.setText(str(self.curr_fit.track_parms[-1][0]))
			self.curr_tab.lbl_track_y0.setText(str(self.curr_fit.track_parms[-1][2]))
			self.curr_tab.lbl_track_c0.setText(str(self.curr_fit.track_parms[-1][1]))
			
			self.curr_tab.curr_slider.setSliderPosition(shape(self.curr_fit.track_fit)[0]-1)
			
			self.canvas.draw()
		
		else: 
			QtGui.QMessageBox.critical(None, "Error","Did not save fitting tracks.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
	def plot_extrap_fit(self):
		
		if self.curr_node.parent().data(0,0).toString()!="Fits":
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if shape(self.curr_fit.fit_av_d)[0]>0:
			pass
		else:
			QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		self.create_plot_tab("fit")
		
		tvec_extrap,f_extrap=pyfdap_fit.extrap_fit_to_end(self.curr_fit,0,6*self.curr_fit.embryo.tvec_data[-1])
		
		if shape(self.curr_embr.ignored)[0]>0:
			if self.curr_fit.fit_ext==1:
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r*',label='ext_data')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.fit_av_d,'r--',label='fit')
			if self.curr_fit.fit_slice==1:
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g*',label='slice_data')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.fit_av_d,'g--',label='fit')
			if self.curr_fit.fit_int==1:
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b*',label='int_data')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.fit_av_d,'b--',label='fit')
		else:
			if self.curr_fit.fit_ext==1:
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r*',label='ext_data')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.fit_av_d,'r--',label='fit')
			if self.curr_fit.fit_slice==1:
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g*',label='slice_data')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.fit_av_d,'g--',label='fit')
			if self.curr_fit.fit_int==1:
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b*',label='int_data')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.fit_av_d,'b--',label='fit')
		
		if self.curr_fit.fit_ext==1:
			self.ax.plot(tvec_extrap,f_extrap,'r:',label='extrap fit')
		if self.curr_fit.fit_slice==1:
			self.ax.plot(tvec_extrap,f_extrap,'g:',label='extrap fit')
		if self.curr_fit.fit_int==1:
			self.ax.plot(tvec_extrap,f_extrap,'b:',label='extrap fit')
		
		#Show ynaught
		ynaught_plot=self.curr_fit.ynaught_opt*ones(shape(tvec_extrap))
		self.ax.plot(tvec_extrap,ynaught_plot,'k-',label='y0_opt')
		
		#Show cnaught
		cnaught_plot=(self.curr_fit.cnaught_opt+self.curr_fit.ynaught_opt)*ones(shape(tvec_extrap))
		self.ax.plot(tvec_extrap,cnaught_plot,'k--',label='c0_opt+y0_opt')
		
		#Show halflife in plot
		f_tau=self.curr_fit.cnaught_opt*exp(-self.curr_fit.k_opt*self.curr_fit.halflife_s)+self.curr_fit.ynaught_opt
		self.ax.plot([self.curr_fit.halflife_s,self.curr_fit.halflife_s],[self.curr_fit.ynaught_opt,f_tau],'k-')
		self.ax.plot([0,self.curr_fit.halflife_s],[f_tau,f_tau],'k-')
		
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		
		self.adjust_canvas()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Plot background values
	
	def plot_bkgd_timeseries(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		if self.curr_bkgd==None:
			
			QtGui.QMessageBox.critical(None, "Error","No bkgd selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return	
		
		else:
			
			for embr in self.curr_mol.embryos:
				
				#Find fitting tvec for background (bkgds normally don't have tvec)
				if len(embr.tvec_data)==len(self.curr_bkgd.bkgd_ext_vec):
					tvec_data=embr.tvec_data
					break
				
			#Plot
			self.ax.plot(tvec_data,self.curr_bkgd.bkgd_ext_vec,'r-',label="ext")		
			self.ax.plot(tvec_data,self.curr_bkgd.bkgd_int_vec,'b-',label="int")		
			self.ax.plot(tvec_data,self.curr_bkgd.bkgd_slice_vec,'g-',label="slice")		
				
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()		
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Plot all dataseries
	
	
	def plot_all_ext_data_ign(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		#Loop through embryos
		for i,embr in enumerate(self.curr_mol.embryos):
		
			if len(embr.ext_av_data_d)>0:
				
				#Grab color
				#c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				
				#Plot
				self.ax.plot(embr.tvec_data,embr.ext_av_data_d,'-',color=c,label=embr.name)
				
		#Loop through backgrounds
		for i,bkgd in enumerate(self.curr_mol.bkgds):
			
			if len(embr.ext_av_data_d)>0:
				
				#Grab color
				c = cm.gray(float(i)/float(len(self.curr_mol.bkgds)),1)
				#summer
				
				tvec_data=None
				#Find fitting tvec for background (bkgds normally don't have tvec)
				for j,embr in enumerate(self.curr_mol.embryos):
					if len(embr.tvec_data)==len(bkgd.bkgd_ext_vec):
						tvec_data=embr.tvec_data
						break
					elif len(embr.tvec_data)>len(bkgd.bkgd_ext_vec):
						tvec_data=embr.tvec_data[0:len(bkgd.bkgd_ext_vec)]
				
				#Plot
				self.ax.plot(tvec_data,bkgd.bkgd_ext_vec,'-',color=c,label=bkgd.name)		
		
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()		
	
	def plot_all_slice_data_ign(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		#Loop through embryos
		for i,embr in enumerate(self.curr_mol.embryos):
			
			if len(embr.slice_av_data_d)>0:
				
				#Grab color
				c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				
				#Plot
				self.ax.plot(embr.tvec_data,embr.slice_av_data_d,'-',color=c,label=embr.name)
				
		#Loop through backgrounds
		for i,bkgd in enumerate(self.curr_mol.bkgds):
			
			if len(embr.slice_av_data_d)>0:
				
				#Grab color
				c = cm.gray(float(i)/float(len(self.curr_mol.bkgds)),1)
				
				#Find fitting tvec for background (bkgds normally don't have tvec)
				for j,embr in enumerate(self.curr_mol.embryos):
					if len(embr.tvec_data)==len(bkgd.bkgd_slice_vec):
						tvec_data=embr.tvec_data
						break
					elif len(embr.tvec_data)>len(bkgd.bkgd_slice_vec):
						tvec_data=embr.tvec_data[0:len(bkgd.bkgd_slice_vec)]
				
				#Plot
				self.ax.plot(tvec_data,bkgd.bkgd_slice_vec,'-',color=c,label=bkgd.name)		
		
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()		
		
	def plot_all_int_data_ign(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		#Loop through embryos
		for i,embr in enumerate(self.curr_mol.embryos):
			
			if len(embr.int_av_data_d)>0:
				
				#Grab color
				c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				
				#Plot
				self.ax.plot(embr.tvec_data,embr.int_av_data_d,'-',color=c,label=embr.name)
				
		#Loop through backgrounds
		for i,bkgd in enumerate(self.curr_mol.bkgds):
			
			if len(embr.int_av_data_d)>0:
				
				#Grab color
				c = cm.gray(float(i)/float(len(self.curr_mol.bkgds)),1)
				
				#Find fitting tvec for background (bkgds normally don't have tvec)
				for j,embr in enumerate(self.curr_mol.embryos):
					if len(embr.tvec_data)==len(bkgd.bkgd_int_vec):
						tvec_data=embr.tvec_data
						break
					elif len(embr.tvec_data)>len(bkgd.bkgd_int_vec):
						tvec_data=embr.tvec_data[0:len(bkgd.bkgd_int_vec)]
						
				#Plot
				self.ax.plot(tvec_data,bkgd.bkgd_int_vec,'-',color=c,label=bkgd.name)
				
		
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()			
	
	def plot_all_ext_data(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		#Loop through embryos
		for i,embr in enumerate(self.curr_mol.embryos):
			
			if len(embr.ext_av_data_d)>0:
				
				#Grab color
				c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				
				#Plot
				if shape(embr.ignored)[0]>0:	
					self.ax.plot(embr.tvec_ignored,embr.ext_av_data_ign,'-',color=c,label=embr.name)
				else:
					self.ax.plot(embr.tvec_data,embr.ext_av_data_d,'-',color=c,label=embr.name)
				
		#Loop through backgrounds
		for i,bkgd in enumerate(self.curr_mol.bkgds):
			
			if len(embr.ext_av_data_d)>0:
				
				#Grab color
				c = cm.gray(float(i)/float(len(self.curr_mol.bkgds)),1)
				
				#Find fitting tvec for background (bkgds normally don't have tvec)
				for j,embr in enumerate(self.curr_mol.embryos):
					#print len(embr.tvec_data),len(bkgd.bkgd_ext_vec)
					if len(embr.tvec_data)==len(bkgd.bkgd_ext_vec):
						tvec_data=embr.tvec_data
						break
					elif len(embr.tvec_data)>len(bkgd.bkgd_ext_vec):
						tvec_data=embr.tvec_data[0:len(bkgd.bkgd_ext_vec)]
					elif len(embr.tvec_data)<len(bkgd.bkgd_ext_vec):
						print len(embr.tvec_data)
						tvec_data=pyfdap_misc.equ_extend_vec(embr.tvec_data,len(bkgd.bkgd_ext_vec)-len(embr.tvec_data))
						print len(embr.tvec_data)
						
				#Plot
				self.ax.plot(tvec_data,bkgd.bkgd_ext_vec,'-',color=c,label=bkgd.name)		
		
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()			
	
	def plot_all_slice_data(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		#Loop through embryos
		for i,embr in enumerate(self.curr_mol.embryos):
			
			if len(embr.slice_av_data_d)>0:
				
				#Grab color
				c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				
				#Plot
				if shape(embr.ignored)[0]>0:	
					self.ax.plot(embr.tvec_ignored,embr.slice_av_data_ign,'-',color=c,label=embr.name)
				else:
					self.ax.plot(embr.tvec_data,embr.slice_av_data_d,'-',color=c,label=embr.name)
				
		#Loop through backgrounds
		for i,bkgd in enumerate(self.curr_mol.bkgds):
			
			if len(embr.slice_av_data_d)>0:
				
				#Grab color
				c = cm.gray(float(i)/float(len(self.curr_mol.bkgds)),1)
				
				#Find fitting tvec for background (bkgds normally don't have tvec)
				for j,embr in enumerate(self.curr_mol.embryos):
					if len(embr.tvec_data)==len(bkgd.bkgd_slice_vec):
						tvec_data=embr.tvec_data
						break
					elif len(embr.tvec_data)>len(bkgd.bkgd_slice_vec):
						tvec_data=embr.tvec_data[0:len(bkgd_slice_vec)]	
				
				#Plot
				self.ax.plot(tvec_data,bkgd.bkgd_slice_vec,'-',color=c,label=bkgd.name)		
		
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()		
		
	def plot_all_int_data(self):
		
		#Create Plot tab
		self.create_plot_tab("data")
		
		#Loop through embryos
		for i,embr in enumerate(self.curr_mol.embryos):
			
			if len(embr.int_av_data_d)>0:
				
				#Grab color
				c = cm.hsv(float(i)/float(len(self.curr_mol.embryos)),1)
				
				#Plot
				if shape(embr.ignored)[0]>0:	
					self.ax.plot(embr.tvec_ignored,embr.int_av_data_ign,'-',color=c,label=embr.name)
				else:
					self.ax.plot(embr.tvec_data,embr.int_av_data_d,'-',color=c,label=embr.name)
				
		#Loop through backgrounds
		for i,bkgd in enumerate(self.curr_mol.bkgds):
			
			if len(embr.int_av_data_d)>0:
				
				#Grab color
				c = cm.gray(float(i)/float(len(self.curr_mol.bkgds)),1)
				
				#Find fitting tvec for background (bkgds normally don't have tvec)
				for j,embr in enumerate(self.curr_mol.embryos):
					if len(embr.tvec_data)==len(bkgd.bkgd_int_vec):
						tvec_data=embr.tvec_data
						break
					elif len(embr.tvec_data)>len(bkgd.bkgd_int_vec):
						tvec_data=embr.tvec_data[0:len(bkgd_int_vec)]
						
				#Plot
				self.ax.plot(tvec_data,bkgd.bkgd_int_vec,'-',color=c,label=bkgd.name)
				
		
		#Legend and draw
		self.ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		self.adjust_canvas()			
					
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Create plot tab
	
	def create_plot_tab(self,plottype,tabname=None):
		
		if not isinstance(tabname,str):
			
			if plottype in ["data","dataandbkgd"]:
				tabname=self.curr_mol.name+"/"+self.curr_embr.name+"/"+plottype+"#1"
			elif plottype=="fit":
				tabname=self.curr_mol.name+"/"+self.curr_embr.name+"/"+self.curr_fit.name+"/"+plottype+"#1"
			elif plottype in ["err","bar"]:
				tabname=self.curr_mol.name+"/"+plottype+"#1"
			
			else:
				tabname="newtab#1"
		
			for i in range(self.plot_tabs.count()):
				if self.plot_tabs.tabText(i)==tabname:
					nme,nmbr=tabname.split("#")
					nmbr=str(int(nmbr)+1)
					tabname=nme+"#"+nmbr
		
		self.curr_tab=QtGui.QWidget()	
		self.plot_tabs.addTab(self.curr_tab,tabname)
			
		h=10000/self.dpi
		v=10000/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.curr_tab)
		
		self.curr_tab.typ=plottype
		
		self.ax = self.fig.add_subplot(111)
		if plottype=="bar":
			self.fig.subplots_adjust(bottom=0.3)
		else:	
			self.fig.subplots_adjust(right=0.75)
			self.ax.set_xlabel("Time (s)")
			self.ax.set_ylabel("Intensity (AU)")
			
		self.tab_axes.append(self.ax)
		self.tab_figs.append(self.fig)
		self.adjust_canvas()
		if self.plot_tabs.tabText(0)=="PlotTab":
		
			self.plot_tabs.removeTab(self.plot_tabs.currentIndex())
			self.plot_tabs.setTabsClosable(True)
		
		self.plot_tabs.setCurrentWidget(self.curr_tab)
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#What happens when tab is changed
	
	def curr_tab_changed(self,value):
			
		self.curr_tab=self.plot_tabs.widget(value)
		if shape(self.tab_axes)[0]>0:
			self.ax=self.tab_axes[value]
			self.fig=self.tab_figs[value]	
		
		if self.curr_tab!=None:
			self.splitter_ver.refresh()
			self.splitter_hor.refresh()
			self.curr_tab.setHidden(True)
			self.curr_tab.setVisible(True)
			self.adjust_canvas()	
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#What happens when tab is closed
	
	def curr_tab_closed(self,value):
			
		self.curr_tab=self.plot_tabs.widget(value)
		if shape(self.tab_axes)[0]>0:
			self.tab_axes.pop(value)
			self.tab_figs.pop(value)
			
		self.plot_tabs.removeTab(value)
		
		if self.plot_tabs.count()==0:
			self.curr_tab=QtGui.QWidget()
			self.first_tab=self.plot_tabs.addTab(self.curr_tab,"PlotTab")
			self.plot_tabs.setTabsClosable(False)
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Adjust canvas if slider changes
	
	def adjust_canvas(self):
		
		if hasattr(self,'fig'):
		
			h=float(self.splitter_hor.sizes()[1])/float(self.dpi)
			v=float(self.splitter_ver.sizes()[0])/float(self.dpi)
			
			if hasattr(self.curr_tab,'curr_slider'):
				
				h_slider=float(self.curr_tab.curr_slider.size().width())/float(self.dpi)
				v_slider=float(self.curr_tab.curr_slider.size().height())/float(self.dpi)
				
				v=v-5*v_slider
			
			self.fig.set_size_inches(h,v,forward=False)
			
			self.canvas.draw()
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Make slider plot
	
	def create_slider_plot_tab(self,plottype):
		
		#Create new tab
		if plottype in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int"]:
			tabname=self.curr_mol.name+"/"+self.curr_bkgd.name+"/"+plottype+"#1"
		elif plottype in ["ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			tabname=self.curr_embr.name+"/"+plottype+"#1"	
		elif plottype=="track_fit":
			tabname=self.curr_embr.name+"/"+self.curr_fit.name+"/"+"track#1"
			
		else:	
			tabname="newtab#1"

		for i in range(self.plot_tabs.count()):
			if self.plot_tabs.tabText(i)==tabname:
				nme,nmbr=tabname.split("#")
				nmbr=str(int(nmbr)+1)
				tabname=nme+"#"+nmbr
		
		self.curr_tab=QtGui.QWidget()	
		self.plot_tabs.addTab(self.curr_tab,tabname)
		
		h=10000/self.dpi
		v=10000/self.dpi
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.curr_tab)
		
		self.curr_tab.curr_slider=QtGui.QSlider(parent=self.curr_tab)
		self.curr_tab.curr_slider.setOrientation(Qt.Horizontal)
		
		self.curr_tab.typ=plottype
		
		self.hbox_arrows = QtGui.QHBoxLayout()
		self.vbox_slider = QtGui.QVBoxLayout()
		self.vbox_slider.addWidget(self.canvas,stretch=1)
		self.vbox_slider.addLayout(self.hbox_arrows,stretch=1)
	
		if plottype in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int"]:
			self.curr_tab.curr_slider.setRange(0,shape(self.curr_bkgd.bkgd_vals_slice)[0]-1)
			self.curr_tab.curr_slider.setSingleStep(1)
			self.connect(self.curr_tab.curr_slider, QtCore.SIGNAL('sliderReleased()'), self.update_slider_bkgd)
			self.curr_tab.lbl_time_slider = QtGui.QLabel("img = 0", self)
			self.curr_tab.lbl_time_slider.setAlignment(Qt.AlignHCenter)
			
			self.vbox_slider_small=QtGui.QVBoxLayout()
			self.vbox_slider_small.addWidget(self.curr_tab.curr_slider,stretch=0)
			self.vbox_slider_small.addWidget(self.curr_tab.lbl_time_slider)
			
			self.curr_tab.btn_right=QtGui.QToolButton()
			self.curr_tab.btn_right.connect(self.curr_tab.btn_right, QtCore.SIGNAL('clicked()'), self.right_img)
			self.curr_tab.btn_right.setArrowType(QtCore.Qt.RightArrow)
			
			self.curr_tab.btn_left=QtGui.QToolButton()
			self.curr_tab.btn_left.connect(self.curr_tab.btn_left, QtCore.SIGNAL('clicked()'), self.left_img)
			self.curr_tab.btn_left.setArrowType(QtCore.Qt.LeftArrow)
			
			self.hbox_arrows.addWidget(self.curr_tab.btn_left)
			self.hbox_arrows.addLayout(self.vbox_slider_small)
			self.hbox_arrows.addWidget(self.curr_tab.btn_right)
			
		if plottype in ["ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			self.curr_tab.curr_slider.setRange(0,shape(self.curr_embr.vals_slice)[0]-1)
			self.curr_tab.curr_slider.setSingleStep(1)
			self.connect(self.curr_tab.curr_slider, QtCore.SIGNAL('sliderReleased()'), self.update_slider_bkgd)
			self.curr_tab.lbl_time_slider = QtGui.QLabel("img = 0", self)
			self.curr_tab.lbl_time_slider.setAlignment(Qt.AlignHCenter)
			
			self.vbox_slider_small=QtGui.QVBoxLayout()
			self.vbox_slider_small.addWidget(self.curr_tab.curr_slider,stretch=0)
			self.vbox_slider_small.addWidget(self.curr_tab.lbl_time_slider)
			
			self.curr_tab.btn_right=QtGui.QToolButton()
			self.curr_tab.btn_right.connect(self.curr_tab.btn_right, QtCore.SIGNAL('clicked()'), self.right_img)
			self.curr_tab.btn_right.setArrowType(QtCore.Qt.RightArrow)
			
			self.curr_tab.btn_left=QtGui.QToolButton()
			self.curr_tab.btn_left.connect(self.curr_tab.btn_left, QtCore.SIGNAL('clicked()'), self.left_img)
			self.curr_tab.btn_left.setArrowType(QtCore.Qt.LeftArrow)
			
			self.hbox_arrows.addWidget(self.curr_tab.btn_left)
			self.hbox_arrows.addLayout(self.vbox_slider_small)
			self.hbox_arrows.addWidget(self.curr_tab.btn_right)
		
		elif plottype=="track_fit":
			self.curr_tab.curr_slider.setRange(0,shape(self.curr_fit.track_fit)[0]-1)
			self.curr_tab.curr_slider.setSingleStep(1)
			self.connect(self.curr_tab.curr_slider, QtCore.SIGNAL('valueChanged(int)'), self.update_slider_track)
			
			#Labels under Slider
			self.curr_tab.lbl_track_name_k = QtGui.QLabel("k = ", self)
			self.curr_tab.lbl_track_name_y0 = QtGui.QLabel("y0 = ", self)
			self.curr_tab.lbl_track_name_c0 = QtGui.QLabel("c0 = ", self)
			self.curr_tab.lbl_track_name_iter = QtGui.QLabel("iter = ", self)
			self.curr_tab.lbl_track_k = QtGui.QLabel("", self)
			self.curr_tab.lbl_track_y0 = QtGui.QLabel("", self)
			self.curr_tab.lbl_track_c0 = QtGui.QLabel("", self)
			self.curr_tab.lbl_track_iter = QtGui.QLabel("", self)
			
			#Add labels to grid sublayout
			self.curr_tab.grid_track=QtGui.QGridLayout()
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_name_k,0,0)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_k,0,1)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_name_c0,0,2)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_c0,0,3)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_name_y0,0,4)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_y0,0,5)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_name_iter,0,6)
			self.curr_tab.grid_track.addWidget(self.curr_tab.lbl_track_iter,0,7)
			
			#Add left/right buttons
			self.curr_tab.btn_right=QtGui.QToolButton()
			self.curr_tab.btn_right.connect(self.curr_tab.btn_right, QtCore.SIGNAL('clicked()'), self.right_track)
			self.curr_tab.btn_right.setArrowType(QtCore.Qt.RightArrow)
			
			self.curr_tab.btn_left=QtGui.QToolButton()
			self.curr_tab.btn_left.connect(self.curr_tab.btn_left, QtCore.SIGNAL('clicked()'), self.left_track)
			self.curr_tab.btn_left.setArrowType(QtCore.Qt.LeftArrow)
			
			#Add Grid + Slider to vbox layout
			self.vbox_slider_small=QtGui.QVBoxLayout()
			self.vbox_slider_small.addWidget(self.curr_tab.curr_slider,stretch=0)
			self.vbox_slider_small.addLayout(self.curr_tab.grid_track)
			
			#Add vbox layout + buttons to hbox layout
			self.hbox_arrows.addWidget(self.curr_tab.btn_left)
			self.hbox_arrows.addLayout(self.vbox_slider_small)
			self.hbox_arrows.addWidget(self.curr_tab.btn_right)
			
			
		self.curr_tab.setLayout(self.vbox_slider)    
		
		self.ax = self.fig.add_subplot(111)
		self.tab_axes.append(self.ax)
		self.tab_figs.append(self.fig)
		self.adjust_canvas()
		
		if self.plot_tabs.tabText(0)=="PlotTab":
		
			self.plot_tabs.removeTab(self.plot_tabs.currentIndex())
			self.plot_tabs.setTabsClosable(True)
		
		self.plot_tabs.setCurrentWidget(self.curr_tab)
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Update plot through slider for bkgd images
	
	def update_slider_bkgd(self):
		value=self.curr_tab.curr_slider.value()
		self.ax.clear()
		self.ax.draw_artist(self.curr_tab.img_axes[value])
		self.curr_tab.lbl_time_slider.setText("#img = "+str(value+1))
		self.curr_tab.setHidden(True)
		self.curr_tab.setVisible(True)
		
		
	def right_img(self):
		value=self.curr_tab.curr_slider.value()
		self.ax.clear()
		self.ax.draw_artist(self.curr_tab.img_axes[value+1])
		self.curr_tab.lbl_time_slider.setText("#img = "+str(value+2))
		self.curr_tab.setHidden(True)
		self.curr_tab.setVisible(True)
		self.curr_tab.curr_slider.setValue(value+1)
	
	def left_img(self):
		value=self.curr_tab.curr_slider.value()
		self.ax.clear()
		self.ax.draw_artist(self.curr_tab.img_axes[value-1])
		self.curr_tab.lbl_time_slider.setText("#img = "+str(value))
		self.curr_tab.setHidden(True)
		self.curr_tab.setVisible(True)
		self.curr_tab.curr_slider.setValue(value-1)
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Update plot through slider for tracking images
	
	def update_slider_track(self,value):
		self.ax.clear()
		
		if shape(self.curr_embr.ignored)[0]>0:
			if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
						
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'g--')
					
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'r--')
				
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
				
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'b--')	
		
		else:
			if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
					
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'g--')
					
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'r--')
				
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
				
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'b--')	
		
		self.curr_tab.lbl_track_k.setText(str(self.curr_fit.track_parms[value][0]))
		self.curr_tab.lbl_track_y0.setText(str(self.curr_fit.track_parms[value][1]))
		self.curr_tab.lbl_track_c0.setText(str(self.curr_fit.track_parms[value][2]))
		self.curr_tab.lbl_track_iter.setText(str(value))
		self.canvas.draw()
	
	def right_track(self):
		#Retrieve value
		value=self.curr_tab.curr_slider.value()
		
		#Check if there is a next frame
		if value+1>=len(self.curr_fit.track_fit):
			return
		
		#Increase value and plot
		value=value+1
		self.ax.clear()
		
		if shape(self.curr_embr.ignored)[0]>0:
			if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
						
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'g--')
					
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'r--')
				
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
				
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'b--')	
		
		else:
			if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
					
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'g--')
					
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'r--')
				
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
				
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'b--')	
		
		#Update values
		self.curr_tab.lbl_track_k.setText(str(self.curr_fit.track_parms[value][0]))
		self.curr_tab.lbl_track_y0.setText(str(self.curr_fit.track_parms[value][1]))
		self.curr_tab.lbl_track_c0.setText(str(self.curr_fit.track_parms[value][2]))
		self.curr_tab.lbl_track_iter.setText(str(value))
		
		#Update slider
		self.curr_tab.curr_slider.setValue(value)
		
		self.canvas.draw()
	
	def left_track(self):
		#Retrieve value
		value=self.curr_tab.curr_slider.value()
		
		#Check if there is a previous frame
		if value-1<0:
			return
		
		#Decrease value and plot
		value=value-11
		self.ax.clear()
		
		if shape(self.curr_embr.ignored)[0]>0:
			if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
						
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.slice_av_data_ign,'g-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'g--')
					
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.ext_av_data_ign,'r-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'r--')
				
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
				
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_embr.int_av_data_ign,'b-')
				self.ax.plot(self.curr_embr.tvec_ignored,self.curr_fit.track_fit[value],'b--')	
		
		else:
			if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
					
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'g--')
					
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
				
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'r--')
				
			elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
				
				self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b-')
				self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[value],'b--')	
		
		#Update values
		self.curr_tab.lbl_track_k.setText(str(self.curr_fit.track_parms[value][0]))
		self.curr_tab.lbl_track_y0.setText(str(self.curr_fit.track_parms[value][1]))
		self.curr_tab.lbl_track_c0.setText(str(self.curr_fit.track_parms[value][2]))
		self.curr_tab.lbl_track_iter.setText(str(value))
		
		#Update slider
		self.curr_tab.curr_slider.setValue(value)
		
		self.canvas.draw()
	
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Perform current fit
	
	def perform_fit(self):
		
		if self.curr_fit_node==None:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return		
		
		#Generate wait popup
		self.wait_popup=pyfdap_subwin.fitting_prog(None)
		self.wait_popup.accepted.connect(self.fitting_canceled)
		self.statusBar().showMessage("Fitting")
		self.setDisabled(True)
		
		#Generate Qthread and pass fitting there
		self.fitting_task=pyfdap_subwin.fitting_thread(embryo=self.curr_embr,fit=self.curr_fit,gui=self)
		self.fitting_task.taskFinished.connect(self.fitting_finished)
		self.fitting_task.start()
				
	def fitting_finished(self):
	
		self.wait_popup.close()
		self.statusBar().showMessage("Idle")
		#Setting fitted=1
		self.curr_embr_node.setText(2,"1")
		self.curr_fit_node.setText(2,"1")
		
		#Plotting
		if self.curr_fit.save_track==1:
			self.plot_track_fit()
		else:
			self.plot_fit()
		
		self.setEnabled(True)
		return
	
	def fitting_canceled(self):
		
		self.setEnabled(True)
		self.statusBar().showMessage("Idle")
		self.fitting_task.terminate()
		
		self.wait_popup.close()
		
	def fit_all_series(self):
		fits_to_fit=[]
		
		#Check if fit is selected
		if self.curr_node.parent().data(0,0).toString()=="Fits":
			
			if hasattr(self,"fits_to_fit"):
				pass
			else:
				#Generate all necessary fits
				fits_to_fit.append(self.curr_fit)
				self.copy_fit()
				if fits_to_fit[-1].fit_ext==1:
					self.curr_fit.fit_ext=0
					self.curr_fit.fit_int=1
				elif fits_to_fit[-1].fit_int==1:
					self.curr_fit.fit_ext=1
					self.curr_fit.fit_int=0
				elif fits_to_fit[-1].fit_slice==1:
					self.curr_fit.fit_ext=1
					self.curr_fit.fit_slice=0
				
				fits_to_fit.append(self.curr_fit)
				self.copy_fit()
				
				if fits_to_fit[-2].fit_ext==1:
					self.curr_fit.fit_ext=0
					self.curr_fit.fit_int=0
					self.curr_fit.fit_slice=1
				elif fits_to_fit[-2].fit_int==1:
					self.curr_fit.fit_ext=0
					self.curr_fit.fit_int=0
					self.curr_fit.fit_slice=1
				elif fits_to_fit[-2].fit_slice==1:
					self.curr_fit.fit_ext=0
					self.curr_fit.fit_slice=0
					self.curr_fit.fit_int=1
				fits_to_fit.append(self.curr_fit)
				
				self.fits_to_fit=fits_to_fit
				print fits_to_fit
			
			#Generate wait popup
			self.wait_popup=pyfdap_subwin.fitting_prog(None)
			self.wait_popup.accepted.connect(self.fitting_canceled)
			self.statusBar().showMessage("Fitting")
			self.setDisabled(True)
			
			#Generate Qthread and pass fitting there
			self.fitting_task=pyfdap_subwin.fitting_all_thread(embryo=self.curr_embr,fits=fits_to_fit,gui=self)
			self.fitting_task.taskFinished.connect(self.fitting_series_finished)
			self.fitting_task.start()
			
		else:
			
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
	def perform_fits_molecule(self):
		
		#Generate wait popup
		self.wait_popup=pyfdap_subwin.fitting_prog(None)
		self.wait_popup.accepted.connect(self.fitting_canceled)
		self.statusBar().showMessage("Fitting")
		self.setDisabled(True)
		
		#Generate Qthread and pass fitting there
		self.fitting_task=pyfdap_subwin.fitting_mol_thread(molecule=self.curr_mol,gui=self)
		self.fitting_task.taskFinished.connect(self.fitting_all_finished)
		self.fitting_task.start()
	
	def fitting_all_finished(self):
		
		self.wait_popup.close()
		self.statusBar().showMessage("Idle")
	
		#Setting fitted=1
		for j in range(self.curr_embryos.childCount()):
			curr_embr_node=self.curr_embryos.child(j)
			curr_embr_node.setText(2,"1")
			curr_fits=curr_embr_node.child(0)
			
			for i in range(curr_fits.childCount()):
				curr_fits.child(i).setText(2,"1")
			
		self.setEnabled(True)
		
		return	
	
	def fitting_series_finished(self):
		
		self.wait_popup.close()
		self.statusBar().showMessage("Idle")
		#Setting fitted=1
		self.curr_embr_node.setText(2,"1")
		for i in range(self.curr_fits.childCount()):
			self.curr_fits.child(i).setText(2,"1")
		
		self.setEnabled(True)

		return	
	
	def compute_av_corr_Fs(self):
		
		if self.curr_mol_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		self.curr_mol=pyfdap_fit.comp_av_corr_F(self.curr_mol)
			
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Show/Hide Console
	
	def show_console(self):	
		self.console.setVisible(True)
		self.splitter_ver.refresh()
		self.splitter_hor.refresh()
		self.adjust_canvas()
		self.curr_conf.term_hidden=False
		
	def hide_console(self):
		self.console.setVisible(False)
		self.splitter_ver.refresh()
		self.splitter_hor.refresh()
		self.adjust_canvas()
		self.curr_conf.term_hidden=True
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Show/Hide Prop List
	
	def show_proplist(self):	
		self.prop_list.setVisible(True)
		self.splitter_ver.refresh()
		self.splitter_hor.refresh()
		self.adjust_canvas()
		self.curr_conf.prop_hidden=False
		
	def hide_proplist(self):
		self.prop_list.setVisible(False)
		self.splitter_ver.refresh()
		self.splitter_hor.refresh()
		self.adjust_canvas()
		self.curr_conf.prop_hidden=True
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Show/Hide Plotting Tab
	
	def show_plottab(self):	
		self.plot_tabs.setVisible(True)
		self.splitter_ver.refresh()
		self.splitter_hor.refresh()
		self.adjust_canvas()
		self.curr_conf.plot_hidden=False
		
	def hide_plottab(self):
		self.plot_tabs.setVisible(False)
		self.splitter_ver.refresh()
		self.splitter_hor.refresh()
		self.curr_conf.plot_hidden=True
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Print current memory usage
	
	def print_mem_usage(self):
		
		if self.curr_embr_node!=None:
			
			if self.curr_fit_node==self.embryos_list.currentItem():
				pyfdap_misc.print_mem_usage(self.curr_fit)
				return
			elif self.curr_noise_node==self.embryos_list.currentItem():
				pyfdap_misc.print_mem_usage(self.curr_noise)
				return
			elif self.curr_pre_node==self.embryos_list.currentItem():
				pyfdap_misc.print_mem_usage(self.curr_pre)
				return
			else:
				pyfdap_misc.print_mem_usage(self.curr_embr)
				return
		
		elif self.curr_bkgd!=None:
		
			if self.curr_bkgd_pre_node==self.embryos_list.currentItem():
				pyfdap_misc.print_mem_usage(self.curr_bkgd_pre)
				return
			else:
				pyfdap_misc.print_mem_usage(self.curr_bkgd)
				return
		
		else:
			pyfdap_misc.print_mem_usage(self.curr_mol)
			return
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Export plots/movies	
	
	def export_plot(self):
		
		#Check if plot is selected
		ind=self.plot_tabs.currentIndex()
		if self.plot_tabs.tabText(ind)=="PlotTab":
			QtGui.QMessageBox.critical(None, "Error","No plot selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Building filter
		file_choices = "IMAGES (*pdf *.png *.tif *eps)"
		
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,self.tr(file_choices),"*.png")
		fn_save=str(fn_save)
		if fn_save=='':
			return
		self.lastopen=os.path.dirname(str(fn_save))
		
		if  self.curr_tab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			ind=self.curr_tab.curr_slider.value()
			matplotlib.image.imsave(fn_save,self.curr_tab.imgs[ind])
		else:
			self.fig.savefig(fn_save)
	
	def export_plot_series(self):
		
		ind=self.plot_tabs.currentIndex()
		if self.plot_tabs.tabText(ind)=="PlotTab":
			QtGui.QMessageBox.critical(None, "Error","No plot selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		if self.curr_tab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int"]:
			
			fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen)
			fn_save=str(fn_save)
			if fn_save=='':
				return
			self.lastopen=os.path.dirname(str(fn_save))
			
			if "." in fn_save:
				fn_save,ending=fn_save.split(".")
			
			files=[]
			
			#Make folder to put images in
			os.mkdir(fn_save)
			
			fn_file=os.path.basename(fn_save)
			
			for i in range(shape(self.curr_tab.img_axes)[0]):
			
				curr_fn_save=fn_save+"/"+fn_file+str(i)+'.tif'
				matplotlib.image.imsave(curr_fn_save,self.curr_tab.imgs[i])
				files.append(curr_fn_save)
			
		else:
			QtGui.QMessageBox.critical(None, "Error","No timeseries selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
	def export_movie(self):
		
		ind=self.plot_tabs.currentIndex()
		if self.plot_tabs.tabText(ind)=="PlotTab":
			QtGui.QMessageBox.critical(None, "Error","No plot selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
			
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.mpg *.avi","*.mpg")
		fn_save=str(fn_save)
		if fn_save=='':
			return
		self.lastopen=os.path.dirname(str(fn_save))
		
		if "." in fn_save:
			fn_save_temp,ending=fn_save.split(".")
		else:
			fn_save_temp=fn_save
			ending="mpg"
			fn_save=fn_save+"."+ending
		
		if ending == "avi" or ending=="mpg":
			encstr="mencoder 'mf://*_tmp.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=mpeg2video  -o "+ fn_save
		else:
			encstr="mencoder 'mf://*_tmp.png' -mf type=png:fps=5 -ovc lavc -lavcopts vcodec=mpeg2video  -o "+ fn_save+'.avi'
		
		files=[]
		if self.curr_tab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int","track_fit"]:	
			if self.curr_tab.typ in ["bkgd_ext","bkgd_slice","bkgd_int","bkgd_masks_embryo","bkgd_masks_ext","bkgd_masks_int","ext","slice","int","masks_embryo","masks_ext","masks_int"]:	
				for i in range(shape(self.curr_tab.img_axes)[0]):
				
					curr_fn_save=fn_save_temp+str(i)+'_tmp'+'.png'
					matplotlib.image.imsave(curr_fn_save,self.curr_tab.imgs[i])
					files.append(curr_fn_save)
			
			if self.curr_tab.typ in ["track_fit"]:
				
				for i in range(shape(self.curr_fit.track_fit)[0]):
					
					self.ax.clear()
				
					if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
							
						self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.slice_av_data_d,'g-')
						self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[i],'g--')
							
					elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
						
						self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.ext_av_data_d,'r-')
						self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[i],'r--')
						
					elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
						
						self.ax.plot(self.curr_embr.tvec_data,self.curr_embr.int_av_data_d,'b-')
						self.ax.plot(self.curr_embr.tvec_data,self.curr_fit.track_fit[i],'b--')	
					
					curr_fn_save=fn_save_temp+str(i)+'_tmp'+'.png'
					self.fig.savefig(curr_fn_save)
					files.append(curr_fn_save)
		
			os.chdir(self.lastopen)
			os.system(encstr)
				
			for fn in files:
				os.remove(fn)

		else:
			QtGui.QMessageBox.critical(None, "Error","No timeseries selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Plot editing
	
	def edit_plot(self):
		
		if not hasattr(self,'ax'):
			QtGui.QMessageBox.critical(None, "Error","No plot tab selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		pyfdap_plot_dialogs.modifyPlotDialog(self.ax,self).exec_()
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Export Embryo object to csv file	
	
	
	def export_embryo_csv(self):
		
		if self.curr_embr_node==None:
			QtGui.QMessageBox.critical(None, "Error","No embryo selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.cvv","*.csv")
		fn_save=str(fn_save)
		if fn_save=='':
			return
		self.lastopen=os.path.dirname(str(fn_save))
		
		pyfdap_misc.write_csv_embryo(fn_save,self.curr_embr)
		
	def export_molecule_csv(self):
		
		if self.curr_mol_node==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.csv","*.csv")
		fn_save=str(fn_save)
		if fn_save=='':
			return
		self.lastopen=os.path.dirname(str(fn_save))
		
		pyfdap_misc.write_csv_molecule(fn_save,self.curr_mol)	
	
	def export_fit_to_csv(self):
		
		if self.curr_fit_node==None:
			QtGui.QMessageBox.critical(None, "Error","No fit selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		if self.curr_fit.k_opt==None: 
			QtGui.QMessageBox.critical(None, "Error","Fit not performed yet.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.csv","*.csv")
		fn_save=str(fn_save)
		if fn_save=='':
			return
		
		if self.curr_fit.fit_slice==1 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==0:
							
			pyfdap_misc.write_csv_timeseries(fn_save,[self.curr_embr.slice_av_data_d,self.curr_fit.fit_av_d],["slice_data","fit"])	
			
			
		elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==1 and self.curr_fit.fit_int==0:
			
			pyfdap_misc.write_csv_timeseries(fn_save,[self.curr_embr.ext_av_data_d,self.curr_fit.fit_av_d],["ext_data","fit"])	
			
		elif self.curr_fit.fit_slice==0 and self.curr_fit.fit_ext==0 and self.curr_fit.fit_int==1:
			
			pyfdap_misc.write_csv_timeseries(fn_save,[self.curr_embr.int_av_data_d,self.curr_fit.fit_av_d],["int_data","fit"])	
			
	def export_errorbar_to_csv(self):
		
		if self.curr_mol==None:
			QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
	
		if self.curr_mol.tvec_errors==None: 
			QtGui.QMessageBox.critical(None, "Error","You haven't made an error bar plot yet. Go to Statistics -> Plotting -> Plot normed fit first.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen,"*.csv","*.csv")
		fn_save=str(fn_save)
		if fn_save=='':
			return
			
		pyfdap_misc.write_csv_timeseries(fn_save,[self.curr_mol.tvec_avg,self.curr_mol.tvec_errors,self.curr_mol.data_av,self.curr_mol.data_errors,self.curr_mol.fit_av],["tvec_avg","tvec_errors","data_av","data_errors","fit_av"])	
		
		return
	
	def export_sel_obj_to_csv(self):
		
		if self.curr_obj==None:
			QtGui.QMessageBox.critical(None, "Error","Nothing selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Property selector
		ret=pyfdap_subwin.select_obj_props(self.curr_obj,self).exec_()
		
		#Filename
		if hasattr(self.curr_obj,'name'):
			fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen+"/"+self.curr_obj.name+".csv","*.csv")
		else:
			fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen+"/"+'saved_obj'+".csv","*.csv")	
		fn_save=str(fn_save)
		self.lastopen=os.path.dirname(str(fn_save))
		
		#Write csv
		pyfdap_misc.write_csv_sel_props([self.curr_obj],fn_save)
		
	def export_sel_fits_to_csv(self):
		
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		if self.curr_mol.sel_fits==[]:
			return
		
		#Property selector
		ret=pyfdap_subwin.select_obj_props(self.curr_mol.sel_fits[0],self).exec_()
		
		#Write same properties in all fit.selected_props
		names=[]
		for fit in self.curr_mol.sel_fits:
			fit.selected_properties=list(self.curr_mol.sel_fits[0].selected_props)
			names.append(fit.embryo.name+"_"+fit.name)
			
		#Filename
		fn_save=QtGui.QFileDialog.getSaveFileName(self, 'Save file', self.lastopen+"/"+self.curr_mol.name+".csv","*.csv")
		fn_save=str(fn_save)
		self.lastopen=os.path.dirname(str(fn_save))
		
		#Write csv
		pyfdap_misc.write_csv_sel_props(self.curr_mol.sel_fits,fn_save,names=names)
		
	#----------------------------------------------------------------------------------------------------------------------------------------
	#Statistics
	
	def sumup_molecule(self,mol=None):
		
		if mol==None:
			if self.curr_mol_node==None:
				QtGui.QMessageBox.critical(None, "Error","No molecule selected.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
				return
			mol=self.curr_mol
		
		if shape(mol.embryos)[0]==0:
			QtGui.QMessageBox.critical(None, "Error","Molecule does not have any embryos to average.",QtGui.QMessageBox.Ok | QtGui.QMessageBox.Default)
			return
		
		#Open Bkgd dialog for bkgd data set
		ret=pyfdap_subwin.select_fits(mol,1,self).exec_()
		
		mol.sumup_results()
	
	def plot_ks_by_fit(self):
		
		#Get some fits and embryos if not already existent
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		ks=[]
		names=[]
		for fit	in self.curr_mol.sel_fits:
			
			ks.append(fit.k_opt)
			names.append(fit.embryo.name)
		
		N=shape(self.curr_mol.sel_fits)[0]
		ind= arange(N)
		width=0.5
		
		self.create_plot_tab("bar")
		self.ax.bar(ind+width, ks, width, color='m')
		self.ax.set_xticks(ind+1.5*width)
		self.ax.set_xticklabels(names,rotation=90)
		self.adjust_canvas()
		
	def plot_ynaughts_by_fit(self):
		
		#Get some fits and embryos if not already existent
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		ynaughts=[]
		names=[]
		for fit	in self.curr_mol.sel_fits:
			
			ynaughts.append(fit.ynaught_opt)
			names.append(fit.embryo.name)
		
		N=shape(self.curr_mol.sel_fits)[0]
		ind= arange(N)
		width=0.5
		
		self.create_plot_tab("bar")
		self.ax.bar(ind+width, ynaughts, width, color='y')
		self.ax.set_xticks(ind+1.5*width)
		self.ax.set_xticklabels(names,rotation=90)
		self.adjust_canvas()	
		
	def plot_cnaughts_by_fit(self):
		
		#Get some fits and embryos if not already existent
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		cnaughts=[]
		names=[]
		for fit	in self.curr_mol.sel_fits:
			
			cnaughts.append(fit.cnaught_opt)
			names.append(fit.embryo.name)
		
		N=shape(self.curr_mol.sel_fits)[0]
		ind= arange(N)
		width=0.5
		
		self.create_plot_tab("bar")
		self.ax.bar(ind+width, cnaughts, width, color='c')
		self.ax.set_xticks(ind+1.5*width)
		self.ax.set_xticklabels(names,rotation=90)
		self.adjust_canvas()		
			
	def plot_all_by_fit(self):
	
		#Get some fits and embryos if not already existent
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		cnaughts=[]
		ynaughts=[]
		ks=[]
		names=[]
		for fit	in self.curr_mol.sel_fits:
			
			ks.append(fit.k_opt)
			ynaughts.append(fit.ynaught_opt)
			cnaughts.append(fit.cnaught_opt)
			names.append(fit.embryo.name)
		
		N=shape(self.curr_mol.sel_fits)[0]
		ind= arange(N)
		width=0.3
		
		self.create_plot_tab("bar")
		
		self.ax.bar(ind-width, ks, width, color='m',label='k')
		self.ax2=self.ax.twinx()
		
		self.ax2.bar(ind, cnaughts, width, color='c',label='c0')
		self.ax2.bar(ind+width, ynaughts, width, color='y',label='y0')
		
		self.ax.legend(loc=2, borderaxespad=0.)	
		self.ax2.legend(loc=1, borderaxespad=0.)	
		
		self.ax.set_xticks(ind+0.5*width)
		self.ax.set_xticklabels(names,rotation=90)
		
		self.adjust_canvas()		
	
	def hist_k(self):
		
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		ks=[]
		for fit	in self.curr_mol.sel_fits:
			
			ks.append(fit.k_opt)
			
		self.hist_parm(ks,'k')
		
	def hist_taumin(self):
		
		if self.curr_mol.sel_fits==[]:
			self.sumup_molecule()
		
		hs=[]
		for fit	in self.curr_mol.sel_fits:
			
			hs.append(fit.halflife_min)
		
		minbin=min(hs)-mod(min(hs),self.bin_width_halfilfe_min)
		maxbin=max(hs)+mod(max(hs),self.bin_width_halfilfe_min)
		
		print minbin,maxbin
		
		bins=arange(minbin-self.bin_width_halfilfe_min,maxbin+self.bin_width_halfilfe_min,self.bin_width_halfilfe_min)
		
		self.hist_parm(hs,'halflife_min',lbl_x="halflife (min)",bin_vec=bins)	
		
	def hist_parm(self,parmvec,name,names=[],lbl_x="",lbl_y="",bin_vec=None):
		
		self.create_plot_tab("bar")	
		
		if bin_vec==None:
			n, bins, patches = self.ax.hist(parmvec, 10, histtype='bar',label=name)
		else:
			n, bins, patches = self.ax.hist(parmvec, bin_vec, histtype='bar',label=name)
		self.ax.set_xlabel(lbl_x)
		self.ax.set_ylabel(lbl_y)
		self.adjust_canvas()
	
	def bar_parm(self,parmvec,names=[],errs=[],lbl_y="",ax=None):
		
		if ax==None:
			self.create_plot_tab("bar")	
			ax=self.ax
		
		ind=arange(len(parmvec))
		width=0.2
		
		errs=asarray(errs)
		
		ax.bar(ind-width,parmvec,width,yerr=errs,error_kw=dict(ecolor='k'))
		
		ax.set_xlim([min(ind)-0.5,max(ind)+0.5])
		
		ax.set_xticks(ind-0.5*width)
		ax.set_xticklabels(names,rotation=90)
		ax.set_ylabel(lbl_y)
		self.adjust_canvas()
				
	def norm_data_fit_plot(self,ax=None,mol=None):
		
		if mol==None:
			mol=self.curr_mol
		
		#Get some fits and embryos if not already existent
		if mol.sel_fits==[]:
			self.sumup_molecule(mol=mol)
		
		#Bin data, pin data and fit
		mol=pyfdap_fit.fit_binned_mol(mol,True,plot=False)
		
		#Draw plot
		self.error_bar_mol(mol=mol,ax=ax)
	
	def unpin_data_fit_plot(self,ax=None,mol=None):
		
		if mol==None:
			mol=self.curr_mol
		
		#Get some fits and embryos if not already existent
		if mol.sel_fits==[]:
			self.sumup_molecule(mol=mol)
		
		#Bin data, pin data and fit
		self.curr_mol=pyfdap_fit.fit_binned_mol(mol,False,plot=False)
		
		#Draw plot
		self.error_bar_mol(mol=mol,ax=ax)
	
	def error_bar_mol(self,ax=None,mol=None):
		
		#Make plot tab
		if ax==None:
			self.create_plot_tab("err")
			ax=self.ax
		
		if mol==None:
			mol=self.curr_mol
			
		last_fit=mol.sel_fits[0]
		
		#Finally plot
		if last_fit.fit_ext==1:
			ax.errorbar(mol.tvec_avg,mol.data_av,xerr=mol.tvec_errors,yerr=mol.data_errors,fmt='ro',label='data_av_ext')
			ax.plot(mol.tvec_avg,mol.fit_av,'k-',label='fit_av_ext')	
		if last_fit.fit_int==1:
			ax.errorbar(mol.tvec_avg,mol.data_av,xerr=mol.tvec_errors,yerr=mol.data_errors,fmt='bo',label='data_av_int')
			ax.plot(mol.tvec_avg,mol.fit_av,'k-',label='fit_av_int')
		if last_fit.fit_slice==1:
			ax.errorbar(mol.tvec_avg,mol.data_av,xerr=mol.tvec_errors,yerr=mol.data_errors,fmt='go',label='data_av_slice')
			ax.plot(mol.tvec_avg,mol.fit_av,'k-',label='fit_av_slice')
			
		#ax.autoscale(enable=True, axis='x', tight=True)
		ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
		#ax.set_ylim([0,1.1])
		ax.set_xlim(0,mol.tvec_avg[-1]+mol.embryos[0].framerate)
		self.adjust_canvas()
	
	def sel_molecules(self):
		
		names=[]
		for mol in self.molecules:
			names.append(mol.name)
		
		sel_dialog=pyfdap_subwin.listSelectorDialog(self,names,leftTitle="Available",rightTitle="Selected",itemsRight=[])
		
		selectedMols=[]
		
		if sel_dialog.exec_():
			selected = sel_dialog.getSelection()
			
			for mol in self.molecules:
				if mol.name in selected:
					selectedMols.append(mol)
			
		return selectedMols
	
	def perform_ttest(self):
		
		selected=self.sel_molecules()
		
		kopt1=selected[0].get_kopts()
		kopt2=selected[1].get_kopts()
		
		stat,pval=pyfdap_stats.ttest_standard(kopt1,kopt2)
		
	def perform_welch(self):
		
		selected=self.sel_molecules()
		
		kopt1=selected[0].get_kopts()
		kopt2=selected[1].get_kopts()
		
		stat,pval=pyfdap_stats.ttest_welch(kopt1,kopt2)
			
	def perform_mann_whitney(self):
		
		selected=self.sel_molecules()
		
		kopt1=selected[0].get_kopts()
		kopt2=selected[1].get_kopts()
		
		stat,pval=pyfdap_stats.mann_whitney_test(kopt1,kopt2)
		
	def perform_wilcoxon(self):
		
		selected=self.sel_molecules()
		
		kopt1=selected[0].get_kopts()
		kopt2=selected[1].get_kopts()
		
		stat,pval=pyfdap_stats.wilcoxon_test(kopt1,kopt2)
		
	def perform_sharipo(self):
		
		pyfdap_stats.sharipo_test(self.curr_mol.get_kopts())
	
	def compare_molecule_unpin(self):
		
		selected=self.sel_molecules()
		
		self.create_plot_tab("err",tabname="Molecule Comparison")
		
		for mol in selected:
			self.unpin_data_fit_plot(ax=self.ax,mol=mol)	
			
	def compare_molecule_norm(self):
		
		selected=self.sel_molecules()
		
		self.create_plot_tab("err",tabname="Molecule Comparison")
		
		for mol in selected:
			self.norm_data_fit_plot(ax=self.ax,mol=mol)	
				
	def compare_molecule_ks(self):
		
		selected=self.sel_molecules()
		
		self.create_plot_tab("bar",tabname="ks Comparison")
		
		ks=[]
		errs=[]
		names=[]
		for mol in selected:
			ks.append(mol.k_av)
			errs.append(mol.k_std)
			names.append(mol.name)
			
		self.bar_parm(ks,names=names,errs=errs,lbl_y="Decay rate (1/s)",ax=self.ax)	
	
#-------------------------------------------------------------------------------------------------------------------------------------
#Main
			
def main():
	
	try:
		ignWarnings=bool(int(sys.argv[1]))
	except:
		ignWarnings=True
		
	try:
		redirect=bool(int(sys.argv[2]))
	except:
		redirect=True
	
	print "Launching PyFDAP GUI with options:"
	print "ignWarnings:" ,ignWarnings
	print "redirect:", redirect
		
	#Creating application	
	app = QtGui.QApplication(sys.argv)
	font=app.font()
	font.setPointSize(12)
	app.setFont(font)
	main_win = pyfdp(ignWarnings=ignWarnings,redirect=redirect)
	main_win.show()

	sys.exit(app.exec_())
		
if __name__ == '__main__':
	main()
