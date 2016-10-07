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

#Parameter fitting module for PyFDAP, including following functions:
#(1) LB_ynaught_func: Returns the LB for y0 according to the F function in Mueller (2012)
#(2) comp_corr_F: Computes F value according to Mueller (2012)
#(3) comp_av_corr_F: Computes average F values over complete molecule according to Mueller (2012)
#(4) opt_parms_fit: Returns optimal fitting parameters according to Mueller (2012)
#(5) correct_ignored_vecs: Adjusts analyzed time series to time vector with ignored values
#(6) interp_fit: Interpolates/Extrapolates fit to have the same length as embryo timevecs without ignored time points
#(7) fdap_fitting: Selects the right optimization settings given by fit object, calls optimzation algorithm and puts results into right objects
#(7) calc_exp_ssd: Objective function of optimzation algorithm. Returns ssd of fit-data

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

import sys
from numpy import *
#from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.optimize as sopt
from pyfdap_stats_module import *
import matplotlib.pyplot as plt
import time

#=====================================================================================================================================
#Module Functions
#=====================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------
#LB scale function for ynaught

def LB_ynaught_func(molecule,embryo,this_fit):
	
	molecule=comp_corr_F(molecule,embryo)
	
	print embryo.noise.noise
	print embryo.F_ext
	
	if embryo.fits[this_fit].fit_ext==1:	
		print type(embryo.pre.pre_ext)
		embryo.fits[this_fit].LB_ynaught=embryo.F_ext*(embryo.pre.pre_ext-embryo.noise.noise)+embryo.noise.noise
			
	if embryo.fits[this_fit].fit_int==1:	
		embryo.fits[this_fit].LB_ynaught=embryo.F_int*(embryo.pre.pre_int-embryo.noise.noise)+embryo.noise.noise
		
	if embryo.fits[this_fit].fit_slice==1:	
		embryo.fits[this_fit].LB_ynaught=embryo.F_slice*(embryo.pre.pre_slice-embryo.noise.noise)+embryo.noise.noise
	
	return embryo

#-------------------------------------------------------------------------------------------------------------------------------------
#Compute correction function F

#def comp_corr_F(molecule,embryo):
	
	#min_ext=[]
	#min_int=[]
	#min_slice=[]
	
	#for bkgd in molecule.bkgds:
		
		#if hasattr(bkgd,'bkgd_ext_vec_ign') and shape(bkgd.bkgd_ext_vec_ign)[0]>0:
			#bkgd_ext_vec_pre=list(bkgd.bkgd_ext_vec_ign)
			#bkgd_ext_vec_pre.insert(0,bkgd.pre.pre_ext)
			#bkgd_slice_vec_pre=list(bkgd.bkgd_slice_vec_ign)
			#bkgd_slice_vec_pre.insert(0,bkgd.pre.pre_slice)
			#bkgd_int_vec_pre=list(bkgd.bkgd_int_vec_ign)
			#bkgd_int_vec_pre.insert(0,bkgd.pre.pre_int)
		#else:
			#bkgd_ext_vec_pre=list(bkgd.bkgd_ext_vec)
			#bkgd_ext_vec_pre.insert(0,bkgd.pre.pre_ext)
			#bkgd_slice_vec_pre=list(bkgd.bkgd_slice_vec)
			#bkgd_slice_vec_pre.insert(0,bkgd.pre.pre_slice)
			#bkgd_int_vec_pre=list(bkgd.bkgd_int_vec)
			#bkgd_int_vec_pre.insert(0,bkgd.pre.pre_int)
		
		#Fnew_ext=(asarray(bkgd_ext_vec_pre)-embryo.noise.noise)/(bkgd.pre.pre_ext-embryo.noise.noise)
		#Fnew_int=(asarray(bkgd_int_vec_pre)-embryo.noise.noise)/(bkgd.pre.pre_int-embryo.noise.noise)
		#Fnew_slice=(asarray(bkgd_slice_vec_pre)-embryo.noise.noise)/(bkgd.pre.pre_slice-embryo.noise.noise)
		
		#min_ext.append(min(Fnew_ext))
		#min_int.append(min(Fnew_int))
		#min_slice.append(min(Fnew_slice))
			
	#embryo.F_ext=mean(min_ext)
	#embryo.F_int=mean(min_int)
	#embryo.F_slice=mean(min_slice)
	
	#return molecule	

def comp_corr_F(molecule,embryo):
	
	regions=["ext","int","slice"]
	
	for r in regions:
	
		try:
			setattr(embryo,'F_'+r,comp_corr_F_region(molecule,embryo,r))
		except TypeError:
			print "Could not compute F value for region ", r, "for embryo", embryo.name
	
	
	
def comp_corr_F_region(molecule,embryo,region):
	
	min_F=[]
	
	#Loop throuh bkgds
	for bkgd in molecule.bkgds:
		
		#Grab right bkgd vector
		if hasattr(bkgd,'bkgd_'+region+'_vec_ign') and shape(getattr(bkgd,'bkgd_'+region+'_vec_ign'))[0]>0:
			bkgd_vec=list(getattr(bkgd,'bkgd_'+region+'_vec_ign'))	
		else:
			bkgd_vec=list(getattr(bkgd,'bkgd_'+region+'_vec'))
		
		#Insert preconversion value at the start
		bkgd_pre=getattr(bkgd.pre,'pre_'+region)	
		bkgd_vec.insert(0,bkgd_pre)
		
		#Grab preconversion value
		pre=getattr(embryo.pre,'pre_'+region)
		
		F=corr_F(bkgd_vec,bkgd_pre,pre,embryo.noise.noise)
		
		min_F.append(min(F))
		
	return mean(min_F)	

def corr_F(bkgd_vec,bkgd_pre,pre,noise):
	return (asarray(bkgd_vec)-noise)/(bkgd_pre-noise)
	

#-------------------------------------------------------------------------------------------------------------------------------------
#Compute average correction function F over molecule object

def comp_av_corr_F(molecule):
	
	F_exts=[]
	F_ints=[]
	F_slices=[]
	
	for embryo in molecule.embryos:
	
		molecule=comp_corr_F(molecule,embryo)
		
		F_exts.append(embryo.F_ext)
		F_ints.append(embryo.F_int)
		F_slices.append(embryo.F_slice)
	
	molecule.F_ext=mean(F_exts)
	molecule.F_int=mean(F_ints)
	molecule.F_slice=mean(F_slices)
	
	return molecule

#-------------------------------------------------------------------------------------------------------------------------------------
#Return optimal parameters according to Mueller (2012)

def opt_parms_fit(fit,molecule):
	
	#Get F value as LB for ynaught
	fit.embryo=LB_ynaught_func(molecule,fit.embryo,fit.fit_number)
	
	#Optimal initial guesses
	if fit.fit_ext==1:
		fit.x0[1]=fit.embryo.ext_av_data_d[0]-fit.embryo.pre.pre_ext
		fit.x0[2]=fit.embryo.ext_av_data_d[0]
	elif fit.fit_int==1:
		fit.x0[1]=fit.embryo.int_av_data_d[0]-fit.embryo.pre.pre_int
		fit.x0[2]=fit.embryo.int_av_data_d[0]
	elif fit.fit_slice==1:
		fit.x0[1]=fit.embryo.slice_av_data_d[0]-fit.embryo.pre.pre_slice
		fit.x0[2]=fit.embryo.slice_av_data_d[0]
	
	#Optimal UB for y0
	try:
		fit.UB_ynaught=max(fit.embryo.int_av_data_d)
	except ValueError: 
		print "Cannot set optimal upper bound for y0. Intracellular concentrations not available."
		
	if fit.UB_ynaught<=fit.x0[2]:
		fit.x0[2]=fit.UB_ynaught
	
	return fit

#-------------------------------------------------------------------------------------------------------------------------------------
#Functions to handle data series with ignored time points

def correct_ignored_vecs(embryo):
	
	if hasattr(embryo,'ignored'):
		pass
	else:
		embryo.ignored=[]
	
	embryo.tvec_ignored=delete(embryo.tvec_data,embryo.ignored)	
	
	if shape(embryo.ext_av_data_d)[0]>0:
		embryo.ext_av_data_ign=delete(embryo.ext_av_data_d,embryo.ignored)	
		embryo.int_av_data_ign=delete(embryo.int_av_data_d,embryo.ignored)	
		embryo.slice_av_data_ign=delete(embryo.slice_av_data_d,embryo.ignored)	
	
	return embryo

def interp_fit(fit):
	print fit.embryo.name
	print len(fit.embryo.tvec_data),len(fit.embryo.tvec_ignored),len(fit.fit_av_d)
	
	fit.fit_av_int=interp(fit.embryo.tvec_data,fit.embryo.tvec_ignored, fit.fit_av_d)
	
	return fit

def extrap_fit_to_end(fit,tstart,tend):
	
	#tstart=fit.embryo.tvec_data[-1]+fit.embryo.framerate
	
	tvec_extrap=arange(tstart,tend+fit.embryo.framerate,fit.embryo.framerate)
	#tvec_extrap=concatenate((fit.embryo.tvec_data,tvec_extrap))
	
	if fit.k_opt!=None:
		f_extrap=fit.cnaught_opt*exp(-fit.k_opt*tvec_extrap)+fit.ynaught_opt
	else:
		print "Warning, you need to perform the fit before you can use this plot."
		return False
	
	return tvec_extrap,f_extrap

	
#-------------------------------------------------------------------------------------------------------------------------------------
#Fits exponential function to data

def fdap_fitting(embryo,this_fit,gui=None):
	
	globals()['embryo']=embryo
	globals()['this_fit']=this_fit
	globals()['gui']=gui
	
	#For good measure, check if ignored vectors are correct
	embryo=correct_ignored_vecs(embryo)
	
	#Counter for function calls
	global iterations
	iterations=0
	
	#Check if constrained and if we need xtransform
	#embryo.fits[this_fit],x0=check_constrained(embryo.fits[this_fit])
	
	#-------------------------------------------------------------------------------------------------------------------------------------
	#Calling optimization algorithms
	#-------------------------------------------------------------------------------------------------------------------------------------
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#fit_cnaught==1 and fit_ynaught==1
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	if embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==1:
		
		#Copying x0 into local variable to pass to solver
		#if embryo.fits[this_fit].transform==0:
		x0=list(embryo.fits[this_fit].x0)
		
		#Getting bounds
		bnds = ((embryo.fits[this_fit].LB_k, embryo.fits[this_fit].UB_k), (embryo.fits[this_fit].LB_cnaught, embryo.fits[this_fit].UB_cnaught),(embryo.fits[this_fit].LB_ynaught,embryo.fits[this_fit].UB_ynaught))
		
		#Calling optimizers
		if embryo.fits[this_fit].opt_meth=='brute':
			kstep=(embryo.fits[this_fit].UB_k-embryo.fits[this_fit].LB_k)/50
			cstep=(embryo.fits[this_fit].UB_cnaught-embryo.fits[this_fit].LB_cnaught)/50
			ystep=(embryo.fits[this_fit].UB_ynaught-embryo.fits[this_fit].LB_ynaught)/50
			
			ranges=(slice(embryo.fits[this_fit].LB_k,embryo.fits[this_fit].UB_k,kstep),slice(embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].UB_cnaught,cstep),slice(embryo.fits[this_fit].LB_ynaught,embryo.fits[this_fit].UB_ynaught,ystep))
			
			res=sopt.brute(calc_exp_ssd, ranges, full_output=bool(embryo.debug_fit),finish=sopt.fmin)
		
		elif embryo.fits[this_fit].opt_meth=='Constrained Nelder-Mead':
			x0=transform_x0(embryo.fits[this_fit].x0,[embryo.fits[this_fit].LB_k,embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].LB_ynaught],[embryo.fits[this_fit].UB_k,embryo.fits[this_fit].UB_cnaught,embryo.fits[this_fit].UB_ynaught])
			res=sopt.fmin(constr_calc_exp_ssd,x0,ftol=embryo.fits[this_fit].opt_tol,maxiter=embryo.fits[this_fit].maxfun,disp=bool(embryo.debug_fit),full_output=True)	
			
		else:
			if embryo.fits[this_fit].opt_meth=='Anneal':
				random.seed(555)
				res=sopt.minimize(calc_exp_ssd, x0, method='Anneal')
			else:
				
				res=sopt.minimize(calc_exp_ssd,x0,method=embryo.fits[this_fit].opt_meth,tol=embryo.fits[this_fit].opt_tol,bounds=bnds,options={'maxiter': embryo.fits[this_fit].maxfun, 'disp': bool(embryo.debug_fit)})
				
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#fit_cnaught==1 and fit_ynaught==0
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				
	elif embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==0:
		
		#Copying x0 into local variable to pass to solver
		#if embryo.fits[this_fit].transform==0:
		x0=list(embryo.fits[this_fit].x0)
		x0.pop(2)
		
		#Getting bounds
		bnds = ((embryo.fits[this_fit].LB_k, embryo.fits[this_fit].UB_k), (embryo.fits[this_fit].LB_cnaught, embryo.fits[this_fit].UB_cnaught))
		
		#Calling optimizers
		if embryo.fits[this_fit].opt_meth=='brute':
			ranges=(slice(embryo.fits[this_fit].LB_k,embryo.fits[this_fit].UB_k,1),slice(embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].UB_cnaught,10))	
			res=sopt.brute(calc_exp_ssd, ranges, full_output=True,finish=sopt.fmin)
		
		elif embryo.fits[this_fit].opt_meth=='Constrained Nelder-Mead':
			x0=transform_x0(embryo.fits[this_fit].x0,[embryo.fits[this_fit].LB_k,embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].LB_ynaught],[embryo.fits[this_fit].UB_k,embryo.fits[this_fit].UB_cnaught,embryo.fits[this_fit].UB_ynaught])
			res=sopt.fmin(constr_calc_exp_ssd,x0,ftol=embryo.fits[this_fit].opt_tol,maxiter=embryo.fits[this_fit].maxfun,disp=bool(embryo.debug_fit),full_output=True)	
		
		else:
			if embryo.fits[this_fit].opt_meth=='Anneal':
				random.seed(555)
				res=sopt.minimize(calc_exp_ssd, x0, method='Anneal')
			else:
				res=sopt.minimize(calc_exp_ssd,x0,method=embryo.fits[this_fit].opt_meth,tol=embryo.fits[this_fit].opt_tol,bounds=bnds,options={'maxiter': embryo.fits[this_fit].maxfun, 'disp': bool(embryo.debug_fit)})
			
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#fit_cnaught==0 and fit_ynaught==1
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
				
	elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==1:
		
		#Copying x0 into local variable to pass to solver
		#if embryo.fits[this_fit].transform==0:
		x0=list(embryo.fits[this_fit].x0)
		x0.pop(1)
		
		#Getting bounds
		bnds = ((embryo.fits[this_fit].LB_k, embryo.fits[this_fit].UB_k), (embryo.fits[this_fit].LB_ynaught, embryo.fits[this_fit].UB_ynaught))
		
		#Calling optimizers
		if embryo.fits[this_fit].opt_meth=='brute':
			ranges=(slice(embryo.fits[this_fit].LB_k,embryo.fits[this_fit].UB_k,1),slice(embryo.fits[this_fit].LB_ynaught,embryo.fits[this_fit].UB_ynaught,1))
			
			res=sopt.brute(calc_exp_ssd, ranges, full_output=bool(embryo.debug_fit),finish=sopt.fmin)
		
		elif embryo.fits[this_fit].opt_meth=='Constrained Nelder-Mead':
			x0=transform_x0(embryo.fits[this_fit].x0,[embryo.fits[this_fit].LB_k,embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].LB_ynaught],[embryo.fits[this_fit].UB_k,embryo.fits[this_fit].UB_cnaught,embryo.fits[this_fit].UB_ynaught])
			res=sopt.fmin(constr_calc_exp_ssd,x0,ftol=embryo.fits[this_fit].opt_tol,maxiter=embryo.fits[this_fit].maxfun,disp=bool(embryo.debug_fit),full_output=True)	
			
		else:
			if embryo.fits[this_fit].opt_meth=='Anneal':
				random.seed(555)
				res=sopt.minimize(calc_exp_ssd, x0, method='Anneal')
			else:	
				res=sopt.minimize(calc_exp_ssd,x0,method=embryo.fits[this_fit].opt_meth,tol=embryo.fits[this_fit].opt_tol,bounds=bnds,options={'maxiter': embryo.fits[this_fit].maxfun, 'disp': bool(embryo.debug_fit)})
			
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#fit_cnaught==0 and fit_ynaught==0
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==0:
		
		#Copying x0 into local variable to pass to solver
		#if embryo.fits[this_fit].transform==0:
		x0=list(embryo.fits[this_fit].x0)
	
		#Getting bounds
		bnds = ((embryo.fits[this_fit].LB_k, embryo.fits[this_fit].UB_k),)
		
		#Calling optimizers
		if embryo.fits[this_fit].opt_meth=='brute':

			ranges=(1,400)
			res=sopt.brute(calc_exp_ssd, (ranges,), full_output=bool(embryo.debug_fit),finish=sopt.fmin)
			
		elif embryo.fits[this_fit].opt_meth=='Constrained Nelder-Mead':
			x0=transform_x0(embryo.fits[this_fit].x0,[embryo.fits[this_fit].LB_k,embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].LB_ynaught],[embryo.fits[this_fit].UB_k,embryo.fits[this_fit].UB_cnaught,embryo.fits[this_fit].UB_ynaught])
			
			print "bla"
			raw_input()
			res=sopt.fmin(constr_calc_exp_ssd,x0,ftol=embryo.fits[this_fit].opt_tol,maxiter=embryo.fits[this_fit].maxfun,disp=bool(embryo.debug_fit),full_output=True)	
			
		
		else:
			if embryo.fits[this_fit].opt_meth=='Anneal':
				random.seed(555)
				res=sopt.minimize(calc_exp_ssd, x0, method='Anneal')
			else:
				
				res=sopt.minimize(calc_exp_ssd,x0,method=embryo.fits[this_fit].opt_meth,tol=embryo.fits[this_fit].opt_tol,bounds=bnds,options={'maxiter': embryo.fits[this_fit].maxfun, 'disp': bool(embryo.debug_fit)})
	
	#-------------------------------------------------------------------------------------------------------------------------------------
	#Saving results in embryo object
	#-------------------------------------------------------------------------------------------------------------------------------------	
	
	if embryo.fits[this_fit].opt_meth=='brute':
		
		if embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==1:
			
			embryo.fits[this_fit].k_opt=res[0]
			embryo.fits[this_fit].cnaught_opt=res[1]
			embryo.fits[this_fit].ynaught_opt=res[2]
			
		elif embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==0:
			embryo.fits[this_fit].k_opt=res[0]
			embryo.fits[this_fit].cnaught_opt=res[1]
			embryo.fits[this_fit].ynaught_opt=embryo.fits[this_fit].x0[2]
			
		elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==1:
			embryo.fits[this_fit].k_opt=res[0]
			embryo.fits[this_fit].cnaught_opt=embryo.fits[this_fit].x0[1]
			embryo.fits[this_fit].ynaught_opt=res[1]
			
		elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==0:
			embryo.fits[this_fit].k_opt=res[0]
			embryo.fits[this_fit].cnaught_opt=embryo.fits[this_fit].x0[1]
			embryo.fits[this_fit].ynaught_opt=embryo.fits[this_fit].x0[2]
			
		embryo.fits[this_fit].ssd=res[1]
		embryo.fits[this_fit].success=True
		embryo.fits[this_fit].halflife_s=log(2)/embryo.fits[this_fit].k_opt
		embryo.fits[this_fit].halflife_min=embryo.fits[this_fit].halflife_s/60
		embryo.fits[this_fit].iterations=iterations
		#In bruteforce iterations = fcalls???
		embryo.fits[this_fit].fcalls=iterations
		
		if embryo.fits[this_fit].fit_ext==1:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.ext_av_data_d,embryo.fits[this_fit].ssd)
		elif embryo.fits[this_fit].fit_int==1:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.int_av_data_d,embryo.fits[this_fit].ssd)
		else:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.slice_av_data_d,embryo.fits[this_fit].ssd)
		
	elif embryo.fits[this_fit].opt_meth=='Constrained Nelder-Mead':
		
		res_new=xtransform(res[0],[embryo.fits[this_fit].LB_k,embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].LB_ynaught],[embryo.fits[this_fit].UB_k,embryo.fits[this_fit].UB_cnaught,embryo.fits[this_fit].UB_ynaught])
		
		if embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==1:
			embryo.fits[this_fit].k_opt=res_new[0]
			embryo.fits[this_fit].cnaught_opt=res_new[1]
			embryo.fits[this_fit].ynaught_opt=res_new[2]
			
		elif embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==0:
			embryo.fits[this_fit].k_opt=res_new[0]
			embryo.fits[this_fit].cnaught_opt=res_new[1]
			embryo.fits[this_fit].ynaught_opt=embryo.fits[this_fit].x0[2]
			
		elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==1:
			embryo.fits[this_fit].k_opt=res_new[0]
			embryo.fits[this_fit].cnaught_opt=embryo.fits[this_fit].x0[1]
			embryo.fits[this_fit].ynaught_opt=res_new[1]
			
		elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==0:
			embryo.fits[this_fit].k_opt=res_new[0]
			embryo.fits[this_fit].cnaught_opt=embryo.fits[this_fit].x0[1]
			embryo.fits[this_fit].ynaught_opt=embryo.fits[this_fit].x0[2]	
		
		embryo.fits[this_fit].ssd=res[1]
		embryo.fits[this_fit].success=not bool(res[4])
		embryo.fits[this_fit].iterations=res[2]
		embryo.fits[this_fit].fcalls=res[3]
		
		if embryo.fits[this_fit].model=="exp":
			embryo.fits[this_fit].halflife_s=log(2)/embryo.fits[this_fit].k_opt
		elif embryo.fits[this_fit].model=="power":
			embryo.fits[this_fit].halflife_s=((2**(embryo.fits[this_fit].npower-1)-1)*embryo.fits[this_fit].cnaught_opt**(1-embryo.fits[this_fit].npower))/(embryo.fits[this_fit].k_opt*(embryo.fits[this_fit].npower-1))
		
		embryo.fits[this_fit].halflife_min=embryo.fits[this_fit].halflife_s/60
		
		if embryo.fits[this_fit].fit_ext==1:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.ext_av_data_d,embryo.fits[this_fit].ssd)
		elif embryo.fits[this_fit].fit_int==1:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.int_av_data_d,embryo.fits[this_fit].ssd)
		else:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.slice_av_data_d,embryo.fits[this_fit].ssd)
		
	else:	
		if embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==1:
			embryo.fits[this_fit].k_opt=res.x[0]
			embryo.fits[this_fit].cnaught_opt=res.x[1]
			embryo.fits[this_fit].ynaught_opt=res.x[2]
			
		elif embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==0:
			embryo.fits[this_fit].k_opt=res.x[0]
			embryo.fits[this_fit].cnaught_opt=res.x[1]
			embryo.fits[this_fit].ynaught_opt=embryo.fits[this_fit].x0[2]
			
		elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==1:
			embryo.fits[this_fit].k_opt=res.x[0]
			embryo.fits[this_fit].cnaught_opt=embryo.fits[this_fit].x0[1]
			embryo.fits[this_fit].ynaught_opt=res.x[1]
			
		elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==0:
			embryo.fits[this_fit].k_opt=res.x[0]
			embryo.fits[this_fit].cnaught_opt=embryo.fits[this_fit].x0[1]
			embryo.fits[this_fit].ynaught_opt=embryo.fits[this_fit].x0[2]
			
		embryo.fits[this_fit].ssd=res.fun
		embryo.fits[this_fit].success=res.success
		embryo.fits[this_fit].iterations=res.nit
		embryo.fits[this_fit].fcalls=res.nfev
		
		if embryo.fits[this_fit].model=="exp":
			embryo.fits[this_fit].halflife_s=log(2)/embryo.fits[this_fit].k_opt
		elif embryo.fits[this_fit].model=="power":
			embryo.fits[this_fit].halflife_s=((2**(embryo.fits[this_fit].npower-1)-1)*embryo.fits[this_fit].cnaught_opt**(1-embryo.fits[this_fit].npower))/(embryo.fits[this_fit].k_opt*(embryo.fits[this_fit].npower-1))
		
		embryo.fits[this_fit].halflife_min=embryo.fits[this_fit].halflife_s/60
		
		if embryo.fits[this_fit].fit_ext==1:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.ext_av_data_d,embryo.fits[this_fit].ssd)
		elif embryo.fits[this_fit].fit_int==1:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.int_av_data_d,embryo.fits[this_fit].ssd)
		else:
			embryo.fits[this_fit].Rsq=fit_Rsq(embryo.slice_av_data_d,embryo.fits[this_fit].ssd)
		
	return embryo

#-------------------------------------------------------------------------------------------------------------------------------------
#Objective function for fdap fitting

def calc_exp_ssd(x):	
	
	#Counting function calls
	global iterations
	iterations=iterations+1
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Checking input
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Check if any variable is negative
	if iterations==0:
		for entry in x:
			if entry<0:
				#If initial guess is out of bounds, just do nothing
				print "Check your initial guess, some entry is negative"
				raw_input("Press Enter to continue")
			
	if embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==1:
		
		knew=x[0]
		cnaught=x[1]
		ynaught=x[2]
		
	elif embryo.fits[this_fit].fit_cnaught==1 and embryo.fits[this_fit].fit_ynaught==0:
	
		knew=x[0]
		cnaught=x[1]
		ynaught=embryo.fits[this_fit].x0[2]
			
	elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==1:
		
		knew=x[0]
		ynaught=x[1]
		cnaught=embryo.fits[this_fit].x0[1]
		
	elif embryo.fits[this_fit].fit_cnaught==0 and embryo.fits[this_fit].fit_ynaught==0:
	
		knew=x[0]
		ynaught=embryo.fits[this_fit].x0[2]
		cnaught=embryo.fits[this_fit].x0[1]
	
	if embryo.debug_fit==1:
		
		print "------------------------------------------"
		print "knew=",knew, "ynaught=", ynaught, "cnaught=", cnaught
	
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Exponential/Power decay
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#If no model selected, take exponential decay
		if hasattr(embryo.fits[this_fit],"model"):
			pass
		else:
			embryo.fits[this_fit].model="exp"
		
	if shape(embryo.ignored)[0]>0:
			
		if embryo.fits[this_fit].model=="exp":
			embryo.fits[this_fit].fit_av_d=cnaught*exp(-knew*embryo.tvec_ignored)+ynaught
		elif embryo.fits[this_fit].model=="power":
			#embryo.fits[this_fit].fit_av_d=((knew*embryo.tvec_ignored+cnaught**((1-embryo.fits[this_fit].npower)/embryo.fits[this_fit].npower))*(embryo.fits[this_fit].npower-1))**(1/(1-embryo.fits[this_fit].npower))+ynaught
			embryo.fits[this_fit].fit_av_d=(cnaught**(1-embryo.fits[this_fit].npower)-knew*embryo.tvec_ignored*(1-embryo.fits[this_fit].npower))**(1/(1-embryo.fits[this_fit].npower))+ynaught
	else:
		if embryo.fits[this_fit].model=="exp":
			embryo.fits[this_fit].fit_av_d=cnaught*exp(-knew*embryo.tvec_data)+ynaught
		elif embryo.fits[this_fit].model=="power":
			#embryo.fits[this_fit].fit_av_d=((knew*embryo.tvec_data+cnaught**((1-embryo.fits[this_fit].npower)/embryo.fits[this_fit].npower))*(embryo.fits[this_fit].npower-1))**(1/(1-embryo.fits[this_fit].npower))+ynaught
			embryo.fits[this_fit].fit_av_d=(cnaught**(1-embryo.fits[this_fit].npower)-knew*embryo.tvec_data*(1-embryo.fits[this_fit].npower))**(1/(1-embryo.fits[this_fit].npower))+ynaught
			
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	#Compute ssd
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	#Residuals
	if shape(embryo.ignored)[0]==0:
		if embryo.fits[this_fit].fit_ext==1 and embryo.fits[this_fit].fit_slice==0 and embryo.fits[this_fit].fit_int==0:
			res=embryo.ext_av_data_d-embryo.fits[this_fit].fit_av_d
		elif embryo.fits[this_fit].fit_ext==0 and embryo.fits[this_fit].fit_slice==1 and embryo.fits[this_fit].fit_int==0:
			res=embryo.slice_av_data_d-embryo.fits[this_fit].fit_av_d
		elif embryo.fits[this_fit].fit_ext==0 and embryo.fits[this_fit].fit_slice==0 and embryo.fits[this_fit].fit_int==1:
			res=embryo.int_av_data_d-embryo.fits[this_fit].fit_av_d
		else:	
			print "You have selected to fit to slice and ext. This won't work"
			raw_input("Wait")
	else:
		if embryo.fits[this_fit].fit_ext==1 and embryo.fits[this_fit].fit_slice==0 and embryo.fits[this_fit].fit_int==0:
			res=embryo.ext_av_data_ign-embryo.fits[this_fit].fit_av_d
		elif embryo.fits[this_fit].fit_ext==0 and embryo.fits[this_fit].fit_slice==1 and embryo.fits[this_fit].fit_int==0:
			res=embryo.slice_av_data_ign-embryo.fits[this_fit].fit_av_d
		elif embryo.fits[this_fit].fit_ext==0 and embryo.fits[this_fit].fit_slice==0 and embryo.fits[this_fit].fit_int==1:
			res=embryo.int_av_data_ign-embryo.fits[this_fit].fit_av_d
		else:	
			print "You have selected to fit to slice and ext. This won't work"
			raw_input("Wait")
	#SSD
	ssd=sum(res**2)
	
	if embryo.debug_fit==1:
		
		if iterations==1:
			if gui==None:
			
				global fig_mon
				fig_mon=plt.figure()
				fig_mon.show()
				global ax_mon
				
				ax_mon=fig_mon.add_subplot(111)
			else:
				pass
		else:
			if gui==None:
				ax_mon.cla()
		
		if gui==None:
			
			if embryo.fits[this_fit].fit_ext==0 and embryo.fits[this_fit].fit_slice==1 and embryo.fits[this_fit].fit_int==0:
				
				ax_mon.plot(embryo.tvec_data,embryo.slice_av_data_d,'g-')
				ax_mon.plot(embryo.tvec_data,embryo.fits[this_fit].fit_av_d,'g--')
				
			elif embryo.fits[this_fit].fit_ext==1 and embryo.fits[this_fit].fit_slice==0 and embryo.fits[this_fit].fit_int==0:
				
				ax_mon.plot(embryo.tvec_data,embryo.ext_av_data_d,'r-')
				ax_mon.plot(embryo.tvec_data,embryo.fits[this_fit].fit_av_d,'r--')
			
			elif embryo.fits[this_fit].fit_ext==0 and embryo.fits[this_fit].fit_slice==0 and embryo.fits[this_fit].fit_int==1:
			
				ax_mon.plot(embryo.tvec_data,embryo.int_av_data_d,'b-')
				ax_mon.plot(embryo.tvec_data,embryo.fits[this_fit].fit_av_d,'b--')
				
			plt.draw()
			plt.pause(0.0001)
			
			#---------------
			#NOTE: this command fixes the automatic live plotting problems with matplotlib v1.3, however generates the error: 
			#/usr/local/lib/python2.7/dist-packages/matplotlib/backend_bases.py:2407: MatplotlibDeprecationWarning: Using default event loop until function specific to this GUI is implemented
			#warnings.warn(str, mplDeprecation)
			#This is not a problem, however annoying
			#---------------
			
		else:
			#---------------
			#Automatic update in plotting doesn't work with gui yet, since there is no plt.pause() for canvas object, need to write something in matplotlib/stackoverflow
			#---------------
			pass
			#gui.canvas.draw()
		
	if embryo.fits[this_fit].save_track==1:
		embryo.fits[this_fit].track_parms.append([knew,cnaught,ynaught])
		embryo.fits[this_fit].track_fit.append(embryo.fits[this_fit].fit_av_d)
		
	return ssd

def constr_calc_exp_ssd(x):
	
	x=xtransform(x,[embryo.fits[this_fit].LB_k,embryo.fits[this_fit].LB_cnaught,embryo.fits[this_fit].LB_ynaught],[embryo.fits[this_fit].UB_k,embryo.fits[this_fit].UB_cnaught,embryo.fits[this_fit].UB_ynaught])
	
	ssd=calc_exp_ssd(x)
	
	return ssd
	
def xtransform(x,LB,UB):
	
	#Determine number of parameters to be fitted
	nparams=len(x)

	
	#Make empty vector
	xtrans = zeros(shape(x))
	
	# k allows some variables to be fixed, thus dropped from the
	# optimization.
	k=0

	for i in range(nparams):

		#Upper bound only
		if UB[i]!=None and LB[i]==None:
			xtrans[i]=UB[i]-x[k]**2
			k=k+1
			
		#Lower bound only	
		elif UB[i]==None and LB[i]!=None:
			xtrans[i]=LB[i]+x[k]**2
			k=k+1
		
		#Both bounds
		elif UB[i]!=None and LB[i]!=None:

			xtrans[i] = (sin(x[k])+1.)/2.*(UB[i] - LB[i]) + LB[i]

			xtrans[i] = max([LB[i],min([UB[i],xtrans[i]])])

			k=k+1
		
		#No bounds
		elif UB[i]==None and LB[i]==None:
			xtrans[i] = x[k]
			k=k+1
			
		#Note: The original file has here another case for fixed variable, but since we made the decision earlier which when we call fdap_fitting, we don't need this here.
	
	return xtrans	
		
def transform_x0(x0,LB,UB):

	# transform starting values into their unconstrained surrogates. Check for infeasible starting guesses.
	x0u = list(x0)
	
	nparams=len(x0)
	
	k=0
	for i in range(nparams):
		
		#Upper bound only
		if UB[i]!=None and LB[i]==None:
			if UB[i]<=x0[i]:
				x0u[k]=0
			else:
				x0u[k]=sqrt(UB[i]-x0[i])	
			k=k+1
			
		#Lower bound only
		elif UB[i]==None and LB[i]!=None:
			if LB[i]>=x0[i]:
				x0u[k]=0
			else:
				x0u[k]=sqrt(x0[i]-LB[i])	
			k=k+1
		
		
		#Both bounds
		elif UB[i]!=None and LB[i]!=None:
			if UB[i]<=x0[i]:
				x0u[k]=pi/2
			elif LB[i]>=x0[i]:
				x0u[k]=-pi/2
			else:
				x0u[k] = 2*(x0[i] - LB[i])/(UB[i]-LB[i]) - 1;
				#shift by 2*pi to avoid problems at zero in fminsearch otherwise, the initial simplex is vanishingly small
				x0u[k] = 2*pi+arcsin(max([-1,min(1,x0u[k])]));
			k=k+1
		
		#No bounds
		elif UB[i]==None and LB[i]==None:
			x0u[k] = x[i]
			k=k+1
	
	return x0u
	
def check_constrained(fit):
	
	constr=False
	
	if fit.fit_cnaught==1 and (fit.LB_ynaught!=None or fit.UB_ynaught!=None):
		constr=True	
	if fit.fit_cnaught==1 and (fit.LB_cnaught!=None or fit.UB_cnaught!=None):
		constr=True	
	
	if constr==True and fit.opt_meth not in ["TNC","L-BFGS-B","SLSQP","brute"]:
		x0=transform_x0(fit.x0,[fit.LB_k,fit.LB_cnaught,fit.LB_ynaught],[fit.UB_k,fit.UB_cnaught,fit.UB_ynaught])		
		fit.transform=1
	else:
		fit.transform=0
		x0=fit.x0
		
	return fit,x0


def fit_binned_mol(mol,pinned,plot=False,fit_cnaught=True):
	
	#Make molecule global so we can do stuff with it
	globals()['mol']=mol
	
	#Grab first fit as reference
	fit=mol.sel_fits[0]
	
	if fit.fit_ext==1 and fit.fit_slice==0 and fit.fit_int==0:
		region="ext"
	elif fit.fit_ext==0 and fit.fit_slice==1 and fit.fit_int==0:
		region="slice"
	elif fit.fit_ext==0 and fit.fit_slice==0 and fit.fit_int==1:
		region="int"
	
	#Get binned vectors and put it into molecule
	tvec_bin,r_bin=bin_tvec_data(mol,pinned,region,plot=plot)
	
	#Pasting into molecule object and compute errors
	mol.tvec_avg=mean_bin(tvec_bin)
	mol.tvec_errors=std_bin(tvec_bin)
	
	setattr(mol,region+"_avg",mean_bin(r_bin))
	setattr(mol,region+"_err",std_bin(r_bin))
	
	#Grab data and error for further use
	mol.data_av=getattr(mol,region+"_avg")
	mol.data_errors=getattr(mol,region+"_err")
	
	#Plot if selected
	if plot:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
	
		ax.errorbar(mol.tvec_avg,mol.data_av,yerr=mol.data_errors,xerr=mol.tvec_errors,fmt='ro')
	
		plt.draw()
		raw_input()
	
	#Define x0 and other paramters
	mol.x0=[fit.k_opt,mol.data_av[0]-mol.data_av[-1],mol.data_av[-1]]
	mol.LB_k=0
	mol.UB_k=None
	mol.LB_ynaught=0
	mol.UB_ynaught=None
	mol.LB_cnaught=0
	mol.UB_cnaught=None
	
	#Transform x0
	mol.x0=transform_x0(mol.x0,[mol.LB_k,mol.LB_cnaught,mol.LB_ynaught],[mol.UB_k,mol.UB_cnaught,mol.UB_ynaught])
	
	#Throw out x0 for cnaught
	x0=list(mol.x0)
	if not fit_cnaught:
		x0.pop(1)
		
	#Pass to optimization algorithm 
	res=sopt.fmin(fit_simple_obj,x0,full_output=True)	
	
	#Put results into molecule to be passed back
	res_new=xtransform(res[0],[mol.LB_k,mol.LB_cnaught,mol.LB_ynaught],[mol.UB_k,mol.UB_cnaught,mol.UB_ynaught])
	
	if fit_cnaught:
		mol.k_opt_refit=res_new[0]
		mol.cnaught_opt_refit=res_new[1]
		mol.ynaught_opt_refit=res_new[2]
	else:	
		mol.k_opt_refit=res_new[0]
		mol.ynaught_opt_refit=res_new[1]
				
	mol.ssd_refit=res[1]
	mol.success_refit=not bool(res[4])
	mol.iterations_refit=res[2]
	mol.fcalls_refit=res[3]
	
	if plot:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
	
		ax.errorbar(mol.tvec_avg,mol.data_av,yerr=mol.data_errors,xerr=mol.tvec_errors,fmt='ro')
		ax.plot(mol.tvec_avg,mol.fit_av,'b-')
		
		plt.draw()
		raw_input()
	
	return mol
	
def fit_simple_obj(x):
	
	#Grab first embryo and fit as reference
	emb=mol.sel_fits[0].embryo
	fit=mol.sel_fits[0]
	
	#Transform for constrained NM
	x=xtransform(x,[mol.LB_k,mol.LB_cnaught,mol.LB_ynaught],[mol.UB_k,mol.UB_cnaught,mol.UB_ynaught])
	
	#Fitting all three parameters here, could add options though	
	if len(x)==2:
		knew=x[0]
		cnaught=mol.x0[1]
		ynaught=x[1]
	elif len(x)==3:
		knew=x[0]
		cnaught=x[1]
		ynaught=x[2]
	
	#Generate model series
	if fit.model=="exp":
		mol.fit_av=cnaught*exp(-knew*mol.tvec_avg)+ynaught
	elif fit.model=="power":
		mol.fit_av=(cnaught**(1-fit.npower)-knew*mol.tvec_avg*(1-fit.npower))**(1/(1-fit.npower))+ynaught
	
	#Compute residuals
	res=mol.data_av-mol.fit_av
	res=res[isfinite(res)]
	
	#SSD
	ssd=sum(res**2)
	
	return ssd

def get_common_tvec(mol):
	
	#Empty lists
	maxs=[]
	mins=[]
	lens=[]
	tvecs=[]
	
	#Find out range of all tvecs and save everything in lists
	for fit in mol.sel_fits:

		emb=fit.embryo
		maxs.append(max(emb.tvec_data))
		mins.append(min(emb.tvec_data))
		lens.append(len(emb.tvec_data))
		dt=emb.tvec_data[-1]-emb.tvec_data[-2]
		
		if len(emb.ignored)>0:
			tvecs.append(emb.tvec_ignored)
		else:
			tvecs.append(emb.tvec_data)

	tmin=min(mins)
	tmax=max(maxs)
	tlen=max(lens)

	#Generating bin vector (adding small percentage of dt to make sure that bounds get into bins too)
	tvec_bin_edges=arange(tmin-0.00001*dt,tmax+1.00001*dt,dt)
	
	return tvec_bin_edges, tvecs

def bin_tvec_data(mol,pinned,region,plot=False):
	
	#Get common tvec
	tvec_bin_edges,tvecs=get_common_tvec(mol)
	
	#Get hist and mapping
	h,mappings=simple_hist(tvec_bin_edges,tvecs,plot=plot)
	
	#Make empty vectors filled with lists
	tvec_bin=list(zeros([len(tvec_bin_edges)-1]))
	r_bin=list(zeros(shape(tvec_bin)))
	
	#Fill in empty lists
	for i in range(len(r_bin)):
		r_bin[i]=[]
		tvec_bin[i]=[]
	
	#Go through each selected fit, grab data vectors, and assign data to bin vector according to mapping from simple_hist
	for i,fit in enumerate(mol.sel_fits):

		m=mappings[i]				
		emb=fit.embryo
		if len(emb.ignored)>0:
			for j in range(len(getattr(emb,region+'_av_data_ign'))):
			
				if pinned:
					r_bin[m[j]].append(pin_dataseries(getattr(emb,region+'_av_data_ign'),fit.ynaught_opt,fit.cnaught_opt)[j])
				else:
					r_bin[m[j]].append(getattr(emb,region+'_av_data_ign')[j])
				
				tvec_bin[m[j]].append(emb.tvec_ignored[j])
				
		else:
			for j in range(len(getattr(emb,region+'_av_data_d'))):
				
				if pinned:
					
					r_bin[m[j]].append(pin_dataseries(getattr(emb,region+'_av_data_d'),fit.ynaught_opt,fit.cnaught_opt)[j])
					
				else:
					r_bin[m[j]].append(getattr(emb,region+'_av_data_d')[j])
				
				tvec_bin[m[j]].append(emb.tvec_data[j])
										
	return tvec_bin,r_bin

def pin_dataseries(datavec,ynaught,cnaught):
	pinned=(asarray(datavec)-ynaught)/cnaught
	return pinned
				
def simple_hist(bins_edges,datavecs,plot=False):
	
	#Empty vector for histogram
	h=zeros([len(bins_edges)-1])
	
	mappings=[]
	
	#Loop through all datavectors
	for data in datavecs:
		
		m=[]
		#Loop through datavector
		for d in data:
			found=False
			#Loop through bin edges
			for i in range(len(bins_edges)-1):
				
				#Check if data point is in this bin
				if bins_edges[i]<=d and d<bins_edges[i+1]:
				
					m.append(i)
					h[i]=h[i]+1
					found=True
					break
				
			if not found:	
				print "not found", d	, bins_edges,m,i
				raw_input()
					
		#Put in mappings list
		mappings.append(m)
	
	if plot:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
		bins=edge_to_bin_vec(bins_edges)
		
		ax.bar(bins,h)
		plt.draw()
		raw_input()
	return h, mappings

def edge_to_bin_vec(edge_vec):
	
	bin_vec=[]
	
	for i in range(len(edge_vec)-1):
		bin_vec.append((edge_vec[i]+edge_vec[i+1])/2)
		
	return bin_vec	

def mean_bin(vec):

	mvec=[]
	for v in vec:
		
		mvec.append(nanmean(v))
	
	mvec=asarray(mvec)
	
	return mvec
	

def std_bin(vec):
	
	mvec=[]
	for v in vec:
		mvec.append(nanstd(v))
	
	mvec=asarray(mvec)
	
	return mvec
							
		
		
		
		
			

