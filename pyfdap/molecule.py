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

#Main molecule class saving all important data, parameters and results for particular molecule.
#Includes subclass bkgd_dataset and pre_dataset
	
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
import gc

#-------------------------------------------------------------------------------------------------------------------------------------
#Molecule object

class molecule:
	
	#Create new molecule object
	def __init__(self,name):
		
		#Data
		self.name=name
		self.embryos=[]
		self.version="1.1a"
		
		#Bkgds
		self.bkgds=[]
		self.bkgd_number=0
		self.bkgd_slice_av=None
		self.bkgd_ext_av=None
		self.bkgd_int_av=None
		
		self.F_ext=None
		self.F_int=None
		self.F_slice=None
		
		self.bkgd_pre_slice=None
		self.bkgd_pre_ext=None
		self.bkgd_pre_int=None
		
		#Results
		self.sel_fits=[]
		self.k_av=None
		self.ynaught_av=None
		self.cnaught_av=None
		self.halflife_s_av=None
		self.halflife_min_av=None
		self.Rsq_av=None
		self.ssd_av=None
		
		self.k_std=None
		self.ynaught_std=None
		self.cnaught_std=None
		self.Rsq_std=None
		self.ssd_std=None
		
		
		self.fitting_parms=[]
			
	def add_embryo(self,embryo):
		
		self.embryos.append(embryo)
		
	def remove_embryo(self,embryo):
		
		self.embryos.remove(embryo)
				
	def add_bkgd(self,bkgd_number,name,mode):
		
		self.bkgd_number=self.bkgd_number+1	
		self.bkgds.append(bkgd_dataset(self,bkgd_number,name,mode))

	def delete_bkgd(self,this_bkgd):
		
		self.bkgd_number=self.bkgd_number-1	
		self.bkgds.pop(this_bkgd)
	
	def save_molecule(self,fn_save):
		
		gc.collect()
		
		#Open and write file using pickle
		with open(fn_save, 'wb') as output:
			pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)
	
	def load_molecule(self,fn_load,update=True):
		if platform.system() in ["Darwin","Linux"]:
			#Open file
			filehandler=open(fn_load, 'r')
			
			#Load file using pickle
			loaded_mol=pickle.load(filehandler)
			
		elif platform.system() in ["Windows"]:
			#Open file
			filehandler=open(fn_load, 'rb')
			
			#Load file using pickle
			loaded_mol=pickle.load(filehandler)
		
		#Update to latest version if selected
		if update:
			loaded_mol.update_version()
		
		return loaded_mol	
	
	def update_version(self):
		
		#Create temporarly a blank molecule file
		moltemp=molecule("temp")
		
		#Update molecule file
		pmisc.update_obj(moltemp,self)
		
		#Update all embryos, pres and fits
		for emb in self.embryos:
			
			emb.update_version()
			
			for fit in emb.fits:
				fit.update_version()
			
			emb.pre.update_version()	
			
		#Update all backgrounds
		for b in self.bkgds:
			b.update_version()
	
	def sumup_results(self):
		
		self.fitting_parms=[]
		
		last_fit=self.sel_fits[0]					
		
		#Check that all fits have roughly the same parameters
		for fit in self.sel_fits:
			
			problem_parms=[]
			
			#Make sure fit object has model
			if hasattr(fit,"model"):
				pass
			else:
				fit.model="exp"
				fit.npower=1
				
			self.fitting_parms.append(fit.opt_meth)
			if fit.model!=last_fit.model:
				problem_parms.append("model")
			if fit.fit_slice!=last_fit.fit_slice:
				problem_parms.append("fit_slice")
			if fit.fit_ext!=last_fit.fit_ext:
				problem_parms.append("fit_ext")
			if fit.fit_int!=last_fit.fit_int:
				problem_parms.append("fit_int")
			if fit.fit_cnaught!=last_fit.fit_cnaught:
				problem_parms.append("fit_c0")
			if fit.fit_ynaught!=last_fit.fit_ynaught:
				problem_parms.append("fit_y0")
			if fit.LB_k==None and last_fit.LB_k==None:
				pass
			elif fit.LB_k!=None and last_fit.LB_k!=None:
				pass
			else:
				problem_parms.append("LB_k")
			if fit.UB_k==None and last_fit.UB_k==None:
				pass
			elif fit.UB_k!=None and last_fit.UB_k!=None:
				pass
			else:
				problem_parms.append("UB_k")
			if fit.LB_cnaught==None and last_fit.LB_cnaught==None:
				pass
			elif fit.LB_cnaught!=None and last_fit.LB_cnaught!=None:
				pass
			else:
				problem_parms.append("LB_cnaught")
			if fit.UB_cnaught==None and last_fit.UB_cnaught==None:
				pass
			elif fit.UB_cnaught!=None and last_fit.UB_cnaught!=None:
				pass
			else:
				problem_parms.append("UB_cnaught")
			if fit.LB_ynaught==None and last_fit.LB_ynaught==None:
				pass
			elif fit.LB_ynaught!=None and last_fit.LB_ynaught!=None:
				pass
			else:
				problem_parms.append("LB_ynaught")
			if fit.UB_ynaught==None and last_fit.UB_ynaught==None:
				pass
			elif fit.UB_ynaught!=None and last_fit.UB_ynaught!=None:
				pass
			else:
				problem_parms.append("UB_ynaught")
				
			if shape(problem_parms)[0]>0:
			
				print "WARNING: Fit", fit.name, "does not have the same parameters than", last_fit.name
				print "Problematic parameters:", problem_parms
				self.fitting_parms=[]
				
						
			if fit.fit_slice==1:
				self.fitting_parms.append("fit_slice")	
			elif fit.fit_ext==1:
				self.fitting_parms.append("fit_ext")
			elif fit.fit_int==1:
				self.fitting_parms.append("fit_int")
			if fit.fit_cnaught==1:
				self.fitting_parms.append("fit_c0")
			if fit.fit_ynaught==1:
				self.fitting_parms.append("fit_y0")
			if fit.LB_k==None:
				self.fitting_parms.append("LB_k=No")
			else:
				self.fitting_parms.append("LB_k=Yes")
			if fit.LB_cnaught==None:
				self.fitting_parms.append("LB_cnaught=No")
			else:
				self.fitting_parms.append("LB_cnaught=Yes")
			if fit.LB_ynaught==None:
				self.fitting_parms.append("LB_ynaught=No")
			else:
				self.fitting_parms.append("LB_ynaught=Yes")
			
			last_fit=fit
		
		#Finally averaging results
		ssd_av_temp=[]
		Rsq_av_temp=[]
		k_av_temp=[]
		ynaught_av_temp=[]
		cnaught_av_temp=[]
		halflife_min_temp=[]
		halflife_s_temp=[]
		for fit in self.sel_fits:
			
			k_av_temp.append(fit.k_opt)
			ynaught_av_temp.append(fit.ynaught_opt)
			cnaught_av_temp.append(fit.cnaught_opt)
			Rsq_av_temp.append(fit.Rsq)
			ssd_av_temp.append(fit.ssd)
		
		
		self.k_av=mean(k_av_temp)
		self.cnaught_av=mean(cnaught_av_temp)
		self.ynaught_av=mean(ynaught_av_temp)
		self.Rsq_av=mean(Rsq_av_temp)
		self.ssd_av=mean(ssd_av_temp)
		
		self.halflife_min_av=log(2)/(self.k_av*60)
		self.halflife_s_av=log(2)/(self.k_av)
		
		self.k_std=std(k_av_temp)
		self.cnaught_std=std(cnaught_av_temp)
		self.ynaught_std=std(ynaught_av_temp)
		self.Rsq_std=std(Rsq_av_temp)
		self.ssd_std=std(ssd_av_temp)
		
	def get_kopts(self):
		
		kopts=[]
		for fit in self.sel_fits:
			kopts.append(fit.k_opt)
	
		return kopts
				
#-------------------------------------------------------------------------------------------------------------------------------------
#bkgd object as subobject of molecule
		
class bkgd_dataset(molecule):
	
	#Create new fit object
	def __init__(self,molecule,bkgd_number,name,mode):
		
		if mode=="empty":
			self.bkgd_number=bkgd_number
			
		elif mode=="default":
			
			#Settings
			self.molecule=molecule
			self.bkgd_number=bkgd_number
			self.name=name
			self.fn_bkgdfolder=""
			self.fn_maskfolder=""
			self.data_ft="tif"
			self.data_enc="uint16"
			self.res_px=512.
			self.centers_embr_px=[]
			self.radiuses_embr_px=[]
			self.bkgd_folder=None
			self.bkgd_slice_av=None
			self.bkgd_ext_av=None
			self.bkgd_int_av=None
			self.bkgd_slice_vec=None
			self.bkgd_ext_vec=None
			self.bkgd_int_vec=None
			self.bkgd_vals_slice=[]
			self.masks_embryo=[]
			self.masks_ext=None
			self.masks_int=None
			self.bkgd_slice=None
			self.bkgd_ext=None
			self.bkgd_int=None
			self.ignored=[]
			self.thresh_masked=False
			
			self.thresh_meth="Otsu"
			self.threshs=[]
			
			#Pre dataset
			self.pre=pre_dataset(self)
			
	#Print out fit
	def print_bkgd(self):
		
		#Going through all attributes of source embryo
		for item in vars(self):
			
			#Don't print out timeseries
			if "_d" not in item:
				print item, " = ", vars(self)[str(item)]
	
	def update_version(self):
		
		#Create temporarly a blank molecule file
		moltemp=molecule("temp")
		bkgdtemp=bkgd_dataset(moltemp,100000,"temp","default")
		
		#Update molecule file
		pmisc.update_obj(bkgdtemp,self)
		
		del bkgdtemp
		gc.collect()
		
		
	
#-------------------------------------------------------------------------------------------------------------------------------------
#pre object as subobject of bkgd
		
class pre_dataset(bkgd_dataset):
	
	#Create new fit object
	def __init__(self,bkgd_dataset):
			
		#Settings
		self.bkgd=bkgd_dataset
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
		
		#Create temporarly a blank molecule file
		moltemp=molecule("temp")
		bkgdtemp=bkgd_dataset(moltemp,100000,"temp","default")
		pretemp=pre_dataset(bkgdtemp)
		
		#Update molecule file
		pmisc.update_obj(pretemp,self)
	