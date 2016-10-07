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

#pyfdap_conf is a simple class to save PyFDAP GUI configurations.

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

from numpy import *
import pickle
from pyfdap_misc_module import *
from pyfdap_img_module import *
from pyfdap_fit_module import *
from pyfdap_stats_module import *
import matplotlib.pyplot as plt
import sys
import platform

#-------------------------------------------------------------------------------------------------------------------------------------
#Config object

class pyfdap_conf:
	
	#Creates new molecule object
	def __init__(self):
		
		#Recently open files
		self.recent=[]
		
		#Last view
		self.plot_hidden=False
		self.term_hidden=False
		self.prop_hidden=False
		self.backup_to_file=False
		self.backup_to_mem=False
		
	
			
	def save_conf(self,fn_save):
		
		#Open and write file using pickle
		with open(fn_save, 'wb') as output:
			pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
	
	def load_conf(self,fn_load):
		if platform.system() in ["Darwin","Linux"]:
			#Open file
			filehandler=open(fn_load, 'r')
			
			#Load file using pickle
			loaded_conf=pickle.load(filehandler)
			
		elif platform.system() in ["Windows"]:
			#Open file
			filehandler=open(fn_load, 'rb')
			
			#Load file using pickle
			loaded_conf=pickle.load(filehandler)
		
		return loaded_conf	
