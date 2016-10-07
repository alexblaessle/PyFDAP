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
#Main embryo class saving all important data, parameters and results for particular molecule.
#Includes subclass pre_dataset, noise_dataset and fit.

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================
from numpy import *
import pickle
import pyfdap_misc_module as pmisc
from pyfdap_img_module import *
from pyfdap_fit_module import *
from pyfdap_stats_module import *
import matplotlib.pyplot as plt
import sys
import platform

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Embryo object

class embryo:
	
	#Creates new embryo object
	def __init__(self,name,mode):
		if mode=="fdap":
			
			#Dataset
			self.name = name
			self.mode = mode
			self.dataset = ""
			#self.emb = ""
			self.fn_datafolder=""    
			self.fn_maskfolder=""    
			self.fn_resultfolder=""
			self.data_enc='uint16'
			self.data_ft='tif'
			self.data_res_px=512.
			self.framerate=10*60
			self.file_list=[]
			self.nframes=2
			self.tstart=0
			self.tend=self.framerate*(self.nframes-1)
			self.post_delay=0
			self.tvec_data=linspace(self.tstart,self.tend,self.nframes)
			self.tvec_data[1:]=self.tvec_data[1:]+self.post_delay
			self.tend=self.tend+self.post_delay
			self.steps_data=self.nframes
			self.dt_data=self.framerate
			self.ignored=[]
			self.tvec_pooled=self.tvec_data
			self.thresh_meth="Otsu"
			self.threshs=[]
			self.fill_mask=0
			self.thresh_masked=False
			
			#Embryo dimensions
			self.centers_embr_px=[]
			self.radiuses_embr_px=[]
			
			#Predataset
			self.pre=pre_dataset(self)
			
			#Analysis results
			self.mask_embryo=None
			self.masks_ext=None
			self.masks_int=None
			self.slice_av_data_d=[]
			self.int_av_data_d=[]
			self.ext_av_data_d=[]
			self.vals_slice=[]
			self.noise=noise_dataset(self)
			
			#Fitting parameters
			self.fit_number=0
			self.fits=[]
			self.debug_fit=0
			
			self.F_ext=None
			self.F_int=None
			self.F_slice=None
			
			#Debugging flags
			self.debug_analysis=0
						
	#Adds new fit				
	def add_fit(self,fit_number,name,mode):
		
		#Increase fit_number
		self.fit_number=self.fit_number+1
		
		#Append new fit
		self.fits.append(fit(self,fit_number,name,mode))
	
	#Removes fit
	def delete_fit(self,this_fit):
		
		#Decrease fit_number
		self.fit_number=self.fit_number-1
		
		#Append new fit
		self.fits.pop(this_fit)

	#Saves embryo
	def save_embryo(self,fn_save):
		#Composing filename
		
		if fn_save==None:
			fn_save=self.fn_resultfolder+"/"+self.name+".pk"
		
		#Open and write file using pickle
		with open(fn_save, 'wb') as output:
			pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
	
	#Load embryo
	def load_embryo(self,fn_load):
		
		if platform.system() in ["Darwin","Linux"]:
			#Open file
			filehandler=open(fn_load, 'r')
			
			#Load file using pickle
			loaded_emb=pickle.load(filehandler)
			
		elif platform.system() in ["Windows"]:
			#Open file
			filehandler=open(fn_load, 'rb')
			
			#Load file using pickle
			loaded_emb=pickle.load(filehandler)
			
			
		return loaded_emb	
	
	#Copies embryo
	def copy_embryo(self,other_embryo):
		
		#Going through all attributes of source embryo
		for item in vars(other_embryo):
			#Don't copy name
			if item is not "name":
				vars(self)[str(item)] = vars(other_embryo)[str(item)]
		
	#Print out embryo
	def print_embryo(self):
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Don't print out large arrays
			if isinstance(vars(self)[str(item)],(int,float,str)):
				print item, " = ", vars(self)[str(item)]
			elif isinstance(vars(self)[str(item)],list) and shape(vars(self)[str(item)])<5:
				print item, " = ", vars(self)[str(item)]
	
	#Print out embryo memory usage
	def print_mem_usage(self):
		
		memusage_name=[]
		memusage_size=[]
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			memusage_name.append(item)
			memusage_size.append(sys.getsizeof(vars(self)[str(item)]))
			
		zipped=zip(memusage_name, memusage_size)
		
		print shape(zipped)
		
		sorted_by_size=sorted(zipped, key=lambda zipped:zipped[1])
		
		for i in range(shape(memusage_name)[0]):
			print sorted_by_size[i][0], sorted_by_size[i][1]
	
	def update_version(self):
		
		#Create temporarly a blank molecule file
		tempemb=embryo("temp","fdap")
		
		#Update embryo file
		pmisc.update_obj(tempemb,self)
		
	def load_csv_timeseries(self,typ,arr):
		
		if typ=="ext":
			self.ext_av_data_d=arr
		elif typ=="int":
			self.int_av_data_d=arr
		elif typ=="slice":
			self.slice_av_data_d=arr
		else:
			print "Warning: Series type " + typ + " unknown." 
			
			
		
		
		
#-------------------------------------------------------------------------------------------------------------------------------------
#Fit object as subobject of embryo
		
class fit(embryo):
	
	#Create new fit object
	def __init__(self,embryo,fit_number,name,mode):
		
		if mode=="empty":
			self.fit_number=fit_number
			
		elif mode=="default":
			
			if embryo.mode=="fdap":
				
				#Settings
				self.embryo=embryo
				self.name=name
				self.fit_number=fit_number
				self.opt_meth="Constrained Nelder-Mead"
				self.fit_ext=1
				self.fit_slice=0
				self.fit_int=0
				self.fit_cnaught=1
				self.fit_ynaught=1
				self.x0=[0,200,50]
				self.LB_cnaught=0.
				self.UB_cnaught=None
				self.LB_ynaught=0.
				self.UB_ynaught=400.
				self.LB_k=0.
				self.UB_k=None
				self.debug_fit=0
				self.maxfun=1000
				self.opt_tol=1e-15
				self.save_track=0
				self.Fynaught=None
				self.model="exp"
				self.npower=2
				self.tvec_pooled=[]
				
				#Results
				self.Rsq=None
				self.ssd=None
				self.k_opt=None
				self.ynaught_opt=None
				self.cnaught_opt=None
				self.success=None
				self.track_parms=[]
				self.track_fit=[]
				self.fit_av_d=[]
				self.iterations=0
				self.fcalls=0
				self.halflife_min=0
				
				
	#Print out fit
	def print_fit(self):
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Don't print out timeseries
			if "_d" not in item:
				print item, " = ", vars(self)[str(item)]
				
	#Print results
	def print_results(self):
		print "=============================="
		print "Results of fit: ", self.name
		res_parms=["k_opt","ynaught_opt","cnaught_opt","halflife_min","success","fcalls","iterations","Rsq","ssd"]
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Don't print out timeseries
			if item in res_parms:
				print item, " = ", vars(self)[str(item)]
		print
		
	def update_version(self):
		
		#Create temporarly a fit and embryo
		tempemb=embryo("temp","fdap")
		tempfit=fit(tempemb,100000,"temp","fdap")
		
		#Update fit file
		pmisc.update_obj(tempfit,self)
		
#-------------------------------------------------------------------------------------------------------------------------------------
#pre object as subobject of embryo
		
class pre_dataset(embryo):
	
	#Create new fit object
	def __init__(self,embryo):
			
		#Settings
		self.embryo=embryo
		self.fn_datafolder=""
		self.fn_maskfolder=""
		self.data_ft="tif"
		self.data_enc="uint16"
		self.res_px=512.
		self.center_embr_px=[256,256]
		self.radius_embr_px=253
		self.pre_slice=None
		self.pre_ext=None
		self.pre_int=None
		self.mask_embryo=[]
		self.mask_ext=None
		self.mask_int=None
		
									
	#Print out fit
	def print_pre(self):
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Don't print out timeseries
			if "_d" not in item:
				print item, " = ", vars(self)[str(item)]
	
	def update_version(self):
		
		#Create temporarly a pre and embryo
		tempemb=embryo("temp","fdap")
		temppre=pre_dataset(tempemb)
		
		#Update fit file
		pmisc.update_obj(temppre,self)
	
#-------------------------------------------------------------------------------------------------------------------------------------
#noise object as subobject of embryo
		
class noise_dataset(embryo):
	
	#Create new fit object
	def __init__(self,embryo):
			
		#Settings
		self.embryo=embryo
		self.mode="out"
		self.fn_datafolder=""
		self.data_ft="tif"
		self.data_enc="uint16"
		self.res_px=512.
		self.noise=None
		
	#Print out fit
	def print_noise(self):
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Don't print out timeseries
			if "_d" not in item:
				print item, " = ", vars(self)[str(item)]							
	
	def update_version(self):
		
		#Create temporarly a pre and embryo
		tempemb=embryo("temp","fdap")
		tempnoise=noise_dataset(tempemb)
		
		#Update fit file
		pmisc.update_obj(tempnoise,self)
	