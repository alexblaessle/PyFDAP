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

#Image analysis module for PyFDAP toolbox, including following functions:
#1) analyze_fdap_data: Analyzes FDAP data given by embryo object
#2) analyze_fdap_bkgds: Analyzes all FDAP bkgd data sets and puts them back into molecule object
#3) get_ext_mask: Uses the Otsu algorithm to compute extracellular and intracellular masks
#4) otsu_imagej: Reimplementation of ImageJ's Otsu algorithm. Adjusted for 16bit images.
#5) get_bkgd: Computes background values for all bkgd objects of a molecule
#6) get_pre: Analyzes pre image data sets
#7) get_embryo_mask: Computes mask of embryo defined by its radius and center
#8) get_noise: Computes noise values via different methods defined in noise subobject
#9) oval_to_circle: Computes center and radius of circle from oval data returned by imageJ

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

#numpy
from numpy import *

#Plotting
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib
import pylab as plab
from pyfdap_misc_module import *

#Misc
import time
import csv
import os, os.path
import sys

#Image processing
import matplotlib.image as mpimg


#=====================================================================================================================================
#Module Functions
#=====================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------
#Load and analyze FDAP data set

def analyze_fdap_data(embryo):
	
	#Get sorted file list
	fn_data_files=get_sorted_folder_list(embryo.fn_datafolder,embryo.data_ft)
	
	#If folder string ends without /, add /
	if embryo.fn_datafolder[-1]=="/":
		pass
	else:
		embryo.fn_datafolder=embryo.fn_datafolder+"/"
	
	#Creating embryo mask
	embryo.masks_embryo=[]
	for i in range(shape(fn_data_files)[0]):
		mask_embryo=get_embryo_mask(embryo.radiuses_embr_px[i],embryo.centers_embr_px[i],embryo.data_res_px,0,fill=embryo.fill_mask)
		embryo.masks_embryo.append(mask_embryo)
		
	#Creating masks for exterior and interior of cells
	masks_ext,masks_int=get_ext_mask(embryo.fn_maskfolder,embryo.data_ft,embryo.thresh_meth,embryo.masks_embryo,embryo.thresh_masked,embryo.threshs,0)
	
	print len(embryo.masks_embryo), len(masks_ext), len(fn_data_files)
	
	print "analyzed masks"
	
	#Some empty result vectors
	ext_av=[]
	int_av=[]
	slice_av=[]
	embryo.vals_slice=[]
	
	#Looping trough all data images
	for i in range(shape(fn_data_files)[0]):
	
		#Load current data img
		fn_load=embryo.fn_datafolder+fn_data_files[i]
		data_img = mpimg.imread(fn_load).astype(embryo.data_enc)
		data_vals=data_img.real
		data_vals=data_vals.astype('float')
		
		#Multiplying with embryo mask and substracting bkgd
		data_vals_slice=data_vals*embryo.masks_embryo[i]
		embryo.vals_slice.append(data_vals_slice)	
		
		#Multiplying with embryo mask and exterior/interior and subtracting bkgd
		data_vals_ext=data_vals*mask_embryo*masks_ext[i]
		data_vals_int=data_vals*mask_embryo*masks_int[i]
		
		#Computing average concentrations for all 3 regions
		slice_av.append(float(sum(data_vals_slice))/float(sum(mask_embryo)))
		ext_av.append(float(sum(data_vals_ext))/float(sum(mask_embryo*masks_ext[i])))
		int_av.append(float(sum(data_vals_int))/float(sum(mask_embryo*masks_int[i])))
		
		print "analyzed", fn_data_files[i]
	
	#Getting pre image
	embryo.pre=get_pre(embryo.pre,0)
	
	print "analyzed pre"
	
	#Getting noise level
	embryo.noise=get_noise(embryo.noise,0)
	
	print "analyzed noise"
	
	#Mapping results back to embryo object
	embryo.masks_ext=masks_ext
	embryo.masks_int=masks_int

	embryo.slice_av_data_d=slice_av
	embryo.int_av_data_d=int_av
	embryo.ext_av_data_d=ext_av
	
	print "Done with analysis of embryo:", embryo.name
	
	return embryo

#-------------------------------------------------------------------------------------------------------------------------------------
#Load and analyze FDAP background data set

def analyze_fdap_bkgds(molecule):

	#Getting background concentrations from seperate bkgd experiments
	molecule.bkgd_slice,molecule.bkgd_ext,molecule.bkgd_int,molecule.bkgd_pre_slice,molecule.bkgd_pre_ext,molecule.bkgd_pre_int,molecule.bkgds=get_bkgd(molecule.bkgds,0)
	
	return molecule

def compute_av_bkgd(molecule):
	
	bkgd_slice=[]
	bkgd_ext=[]
	bkgd_int=[]
	bkgd_pre_slice=[]
	bkgd_pre_ext=[]
	bkgd_pre_int=[]
	
	for bkgd in molecule.bkgds:
		bkgd_slice.append(bkgd.bkgd_slice_av)
		bkgd_ext.append(bkgd.bkgd_ext_av)
		bkgd_int.append(bkgd.bkgd_int_av)
		bkgd_pre_slice.append(bkgd.pre.pre_slice)
		bkgd_pre_ext.append(bkgd.pre.pre_ext)
		bkgd_pre_int.append(bkgd.pre.pre_int)
	

	molecule.bkgd_slice=mean_ignore_none(bkgd_slice)
	molecule.bkgd_ext=mean_ignore_none(bkgd_ext)
	molecule.bkgd_int=mean_ignore_none(bkgd_int)
	molecule.bkgd_pre_slice=mean_ignore_none(bkgd_pre_slice)
	molecule.bkgd_pre_ext=mean_ignore_none(bkgd_pre_ext)
	molecule.bkgd_pre_int=mean_ignore_none(bkgd_pre_int)
	
	return molecule

def mean_ignore_none(vec):
	
	try:
		return mean(vec)
	except TypeError:
		return None
		
#-------------------------------------------------------------------------------------------------------------------------------------
#Take folder where masks are stored and use Otsu algorithm to find contours
		
def get_ext_mask(fn_maskfolder,data_ft,thresh_meth,masks_embryo,thresh_masked,threshs,debug_opt):
	print "in ext mask"
	#Get sorted file list
	fn_mask_files=get_sorted_folder_list(fn_maskfolder,data_ft)
	
	#If folder string ends without /, add /
	if fn_maskfolder[-1]=="/":
		pass
	else:
		fn_maskfolder=fn_maskfolder+"/"
	
	#Empty list for masks
	masks_ext=[]
	masks_int=[]
	
	#Loop through all mask files and get masks
	for i in range(shape(fn_mask_files)[0]):
		
		#Compose filename
		fn_load=fn_maskfolder+fn_mask_files[i]
	
		#Some debugging plots
		if debug_opt==1:
			
			data_img = mpimg.imread(fn_load).astype('float')
			mask_img=data_img.real
			mask_img=mask_img.astype('float')
			
			#Plot image
			fig=plt.figure()
			fig.show()
			ax=fig.add_subplot(131)
			con=ax.contourf(mask_img)
			plt.colorbar(con)
			plt.title("Original image")
			plt.draw()
		
		#Load current image file as grayscale image
		mask_img = mpimg.imread(fn_load).astype("uint16")
		
		#If thresh_masked is selected, apply masking before threshholding
		if thresh_masked:
			mask_img=mask_img*masks_embryo[i]
		
		#Threshholding
		if thresh_meth=="Otsu":
			opt_thresh2,mask2=otsu_imagej(mask_img,1,0,0)
		elif thresh_meth=="Adaptive":
			#NOTE: still developmental
			mask2=adaptive_thresh(mask_img,5)
		elif thresh_meth=="Manual":
			#NOTE: still developmental
			mask2=fixed_thresh(mask_img,threshs[i])
		
		if debug_opt==1:
			print "curr_img=", fn_load
			print "opt_thresh cv2", opt_thresh
			print "opt_thresh imagej", opt_thresh2
			print "opt_thresh own", opt_thresh3
			
			fig=plt.figure()
			fig.show()
			ax=fig.add_subplot(131)
			ax.contourf(mask)
			ax=fig.add_subplot(132)
			ax.contourf(mask2)
			ax=fig.add_subplot(133)
			ax.contourf(mask-mask2)
			
			plt.draw()
		
			raw_input()
		
		mask=mask2
		
		#Some debugging plots
		if debug_opt==1:
			#Plot found contours
			ax=fig.add_subplot(132)
			con=ax.contourf(mask)
			plt.colorbar(con)
			plt.title("Mask")
			plt.draw()
		
		#Inverting contours
		mask_int=mask
		mask_ext=(1-mask)
			
		#Some debugging plots
		if debug_opt==1:
			ax=fig.add_subplot(133)
			con=ax.contourf(mask)
			plt.colorbar(con)
			plt.title("Inverted Mask")
			plt.draw()
			raw_input()
		
		#Appending new mask
		masks_ext.append(mask_ext)	
		masks_int.append(mask_int)	
	
	print "done ext mask"
	
	return masks_ext, masks_int	

#-------------------------------------------------------------------------------------------------------------------------------------
#Scikit adaptive threshholding

def adaptive_thresh(img,block_size,offset=10):
	from skimage.filters import threshold_adaptive
	mask = threshold_adaptive(img, block_size, offset=10)
	return mask

#-------------------------------------------------------------------------------------------------------------------------------------
#Fixed threshholding

def fixed_thresh(img,thresh):
	mask=zeros(shape(img))
	mask=where(img>thresh,mask,1)
	return mask

#-------------------------------------------------------------------------------------------------------------------------------------
#Otsu algorithm of imagej

def otsu_imagej(img,maxval,minval,debug_opt):
	
	#Initialize values
	#L = img.max()
	L = 256
	S = 0 
	N = 0
	
	#Make bin vector (histogram autobinning does not work if there are nans in img)
	bins=linspace(nanmin(img),nanmax(img),L+1)
	
	#Compute histogram
	data,bin_edges=histogram(img,bins)
	bin_width=diff(bin_edges)[0]
	
	#Debugging plot for histogram
	if debug_opt==1:
		bin_vec=arange(L)
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(121)
		ax.bar(bin_vec,data)
		ax.plot(bin_vec,other_data,'g')
		#ax.plot(bin_vec,data,'r')
		
		plt.draw()
		
	for k in range(L):
		#Total histogram intensity
		S = S+ k * data[k]
		#Total number of data points
		N = N + data[k]		
	

	#Temporary variables
	Sk = 0
	BCV = 0
	BCVmax=0
	kStar = 0
	
	#The entry for zero intensity
	N1 = data[0] 
	
	#Look at each possible threshold value,
	#calculate the between-class variance, and decide if it's a max
	for k in range (1,L-1): 
		#No need to check endpoints k = 0 or k = L-1
		Sk = Sk + k * data[k]
		N1 = N1 + data[k]

		#The float casting here is to avoid compiler warning about loss of precision and
		#will prevent overflow in the case of large saturated images
		denom = float(float(N1) * (float(N) - float(N1))) 

		if denom != 0:
			#Float here is to avoid loss of precision when dividing
			num = float(( float(N1) / float(N) ) * S) - Sk 
			BCV = float((num * num)) / float(denom)
		
		else:
			BCV = 0

		if BCV >= BCVmax: 
			#Assign the best threshold found so far
			BCVmax = BCV
			kStar = k
	
	kStar=bin_edges[0]+kStar*bin_width
	
	#Now manipulate the image
	bin_img=zeros(shape(img))
	for i in range(shape(img)[0]):
		for j in range(shape(img)[1]):
			if isnan(img[i,j]):
				bin_img[i,j]=minval
			else:
				if img[i,j]<=kStar:
					bin_img[i,j]=minval
				else:				
					bin_img[i,j]=maxval
	
	if debug_opt==1:
		print "Optimal threshold = ", kStar
		print "#Pixels above threshold = ", sum(bin_img)/float(maxval)
		print "#Pixels below threshold = ", shape(img)[0]**2-sum(bin_img)/float(maxval)
		raw_input()	
		ax2=fig.add_subplot(122)
		ax2.contourf(bin_img)
		plt.draw()
		raw_input()
			
	return kStar,bin_img
	
#-------------------------------------------------------------------------------------------------------------------------------------
#Takes folder where background images are stored and computes background concentration
		
def get_bkgd(bkgds,debug_opt):
	
	#Note: Multiple bkgd objects can be given in bkgds_val
	bkgds_val_slice_temp=[]
	bkgds_val_ext_temp=[]
	bkgds_val_int_temp=[]
	
	bkgds_pre_ext=[]
	bkgds_pre_int=[]
	bkgds_pre_slice=[]
	
	
		
	
	for j in range(shape(bkgds)[0]):
		
		#Get current bkgd folder
		curr_bkgd=bkgds[j]
		print "Analyzing", curr_bkgd.name
		#Get sorted file list
		fn_bkgd_files=get_sorted_folder_list(curr_bkgd.fn_bkgdfolder,curr_bkgd.data_ft)
		
		#If folder string ends without /, add /
		if curr_bkgd.fn_bkgdfolder[-1]=="/":
			pass
		else:
			curr_bkgd.fn_bkgdfolder=curr_bkgd.fn_bkgdfolder+"/"
		
		#Some empty vectors to keep track of bkgd intensity in each time step
		bkgds_val_slice=[]
		bkgds_val_ext=[]
		bkgds_val_int=[]
		
		#Get embryo mask
		curr_bkgd.masks_embryo=[]
		for i in range(shape(fn_bkgd_files)[0]):
			mask_embryo=get_embryo_mask(curr_bkgd.radiuses_embr_px[i],curr_bkgd.centers_embr_px[i],curr_bkgd.res_px,0)
			curr_bkgd.masks_embryo.append(mask_embryo)
		
		#Create Masks using Otsu algorithm
		masks_ext,masks_int=get_ext_mask(curr_bkgd.fn_maskfolder,curr_bkgd.data_ft,curr_bkgd.thresh_meth,curr_bkgd.masks_embryo,curr_bkgd.thresh_masked,curr_bkgd.threshs,0)
		
		#Save masks to bkgd object
		curr_bkgd.masks_ext=masks_ext
		curr_bkgd.masks_int=masks_int
		
		#Create empty vector to save image values to current bkgd
		curr_bkgd.bkgd_vals_slice=[]
		
		#Create empty vector to save embryo masks to current bkgd
		curr_bkgd.masks_embryo=[]
		
		#Loop through all bkgd images and compute bkgd concentration
		for i in range(shape(fn_bkgd_files)[0]):
			
			#Load current image file as color image
			fn_load=curr_bkgd.fn_bkgdfolder+fn_bkgd_files[i]
			
			bkgd_img = mpimg.imread(fn_load).astype(curr_bkgd.data_enc)
			bkgd_vals=bkgd_img.real
			bkgd_vals=bkgd_vals.astype('float')
			
			#Multiplying with embryo mask
			bkgd_vals_slice=mask_embryo*bkgd_vals
			bkgd_vals_ext=mask_embryo*masks_ext[i]*bkgd_vals
			bkgd_vals_int=mask_embryo*masks_int[i]*bkgd_vals
			
			if debug_opt==1:
			
				fig=plt.figure()
				fig.show()
				ax=fig.add_subplot(131)
				ax.contourf(bkgd_vals)
				ax=fig.add_subplot(132)
				ax.contourf(mask_embryo)
				ax=fig.add_subplot(133)
				ax.contourf(bkgd_vals_slice)
				plt.draw()
				raw_input()
				
			curr_bkgd.bkgd_vals_slice.append(bkgd_vals_slice)
			
			#Computing background in slice
			bkgd_slice=sum(bkgd_vals_slice)/sum(mask_embryo)
			bkgd_ext=sum(bkgd_vals_ext)/sum(mask_embryo*masks_ext[i])
			bkgd_int=sum(bkgd_vals_int)/sum(mask_embryo*masks_int[i])
			
			#Appending new bkgd	
			bkgds_val_slice.append(bkgd_slice)
			bkgds_val_ext.append(bkgd_ext)
			bkgds_val_int.append(bkgd_int)
			
			print "Analyzed ", fn_bkgd_files[i]
			
		#Averaging over all bkgds_val
		bkgds_val_slice_temp.append(float(sum(bkgds_val_slice))/float(shape(bkgds_val_slice)[0]))
		bkgds_val_ext_temp.append(float(sum(bkgds_val_ext))/float(shape(bkgds_val_ext)[0]))
		bkgds_val_int_temp.append(float(sum(bkgds_val_int))/float(shape(bkgds_val_int)[0]))
		
		curr_bkgd.bkgd_slice_av=bkgds_val_slice_temp[j]
		curr_bkgd.bkgd_ext_av=bkgds_val_ext_temp[j]
		curr_bkgd.bkgd_int_av=bkgds_val_int_temp[j]
		
		curr_bkgd.bkgd_slice_vec=bkgds_val_slice
		curr_bkgd.bkgd_ext_vec=bkgds_val_ext
		curr_bkgd.bkgd_int_vec=bkgds_val_int
		
		#Getting pre image
		curr_bkgd.pre=get_pre(curr_bkgd.pre,0)
		
		bkgds_pre_slice.append(curr_bkgd.pre.pre_slice)
		bkgds_pre_ext.append(curr_bkgd.pre.pre_ext)
		bkgds_pre_int.append(curr_bkgd.pre.pre_int)
		
		print "Analyzed background dataset:", curr_bkgd.name
	
	#Averaging over all bkgds_val_temp
	bkgd_slice=float(sum(bkgds_val_slice_temp))/float(shape(bkgds_val_slice_temp)[0])
	bkgd_ext=float(sum(bkgds_val_ext_temp))/float(shape(bkgds_val_ext_temp)[0])
	bkgd_int=float(sum(bkgds_val_int_temp))/float(shape(bkgds_val_int_temp)[0])
	
	#Averaging over all bkgds_val_temp
	bkgd_pre_slice=float(sum(bkgds_pre_slice))/float(shape(bkgds_pre_slice)[0])
	bkgd_pre_ext=float(sum(bkgds_pre_ext))/float(shape(bkgds_pre_ext)[0])
	bkgd_pre_int=float(sum(bkgds_pre_int))/float(shape(bkgds_pre_int)[0])
	
	return bkgd_slice,bkgd_ext,bkgd_int, bkgd_pre_slice,bkgd_pre_ext,bkgd_pre_int, bkgds

#-------------------------------------------------------------------------------------------------------------------------------------
#Takes folder where pre images are stored and computes background concentration of pre images
		
def get_pre(pre,debug_opt):
	
	print "in pre"
	
	#Get sorted file list
	fn_pre_files=get_sorted_folder_list(pre.fn_datafolder,pre.data_ft)
	
	#If folder string ends without /, add /
	if pre.fn_datafolder[-1]=="/":
		pass
	else:
		pre.fn_datafolder=pre.fn_datafolder+"/"
	
	pres_slice=[]
	pres_ext=[]
	pres_int=[]
	
	#Creating pre mask
	pre.masks_embryo=[]
	mask_embryo=get_embryo_mask(pre.radius_embr_px,pre.center_embr_px,pre.res_px,0)
	pre.masks_embryo.append(mask_embryo)
	
	#Get extracellular/intracellular mask
	if hasattr(pre,'bkgd'):
		masks_ext,masks_int=get_ext_mask(pre.fn_maskfolder,pre.data_ft,pre.bkgd.thresh_meth,pre.masks_embryo,pre.bkgd.thresh_masked,pre.bkgd.threshs,0)
	else:
		masks_ext,masks_int=get_ext_mask(pre.fn_maskfolder,pre.data_ft,pre.embryo.thresh_meth,pre.masks_embryo,pre.embryo.thresh_masked,pre.embryo.threshs,0)
	
	pre.masks_ext=masks_ext
	pre.masks_int=masks_int
	
	
	#Check if folder is empty
	if len(fn_pre_files)<1:
		pre.pre_ext=0.
		pre.pre_slice=0.
		pre.pre_int=0.
		print "WARNING: Can't find any pre-conversion images in folder " + pre.fn_datafolder + " with filetype " + pre.data_ft + "!"
		return pre
	
	#Loop through all pre images and compute pre concentration
	for i in range(shape(fn_pre_files)[0]):
		
		#Load current image file as color image
		fn_load=pre.fn_datafolder+fn_pre_files[i]
		
	
		pre_img = mpimg.imread(fn_load).astype(pre.data_enc)
		pre_vals=pre_img.real
		pre_vals=pre_vals
		
		#Multiplying with pre mask
		pre_vals_slice=mask_embryo*pre_vals
		pre_vals_ext=mask_embryo*masks_ext[i]*pre_vals
		pre_vals_int=mask_embryo*masks_int[i]*pre_vals
		
		#Computing background in slice
		pres_slice.append(float(sum(pre_vals_slice))/float(sum(mask_embryo)))
		pres_ext.append(float(sum(pre_vals_ext))/float(sum(mask_embryo*masks_ext[i])))
		pres_int.append(float(sum(pre_vals_int))/float(sum(mask_embryo*masks_int[i])))
	
	#Averaging over all pre
	pre_slice=float(sum(pres_slice))/float(shape(pres_slice)[0])
	pre_ext=float(sum(pres_ext))/float(shape(pres_ext)[0])
	pre_int=float(sum(pres_int))/float(shape(pres_int)[0])
	
	pre.pre_ext=pre_ext
	pre.pre_slice=pre_slice
	pre.pre_int=pre_int
	
	return pre
		
#-------------------------------------------------------------------------------------------------------------------------------------
#Takes radius and resolution and returns mask of embryo
		
def get_embryo_mask(radius,center,res,debug_opt,fill=0):
	
	
	
	#Converting res to int if not already
	res=int(res)
	
	mask_embryo=zeros((res,res))
	
	#Going through res*res and see what's within embryo boundaries
	for i in range(res):
		for j in range(res):
			if sqrt((i-center[0])**2 + (j-center[1])**2) < radius:
				mask_embryo[i,j]=1
			else:
				mask_embryo[i,j]=fill
				
	mask_embryo=mask_embryo.T
	
	#Debugging plot of embryo mask
	if debug_opt==1:
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
		con=ax.contourf(mask_embryo)
		plt.colorbar(con)
		plt.title("Embryo mask")
		plt.draw()
		raw_input()	
	
	return mask_embryo

#-------------------------------------------------------------------------------------------------------------------------------------
#Generate noise img dataset

def get_noise(noise,debug_opt):
	
	if noise.mode=="seperate":
	
		fn_noise_files=get_sorted_folder_list(noise.fn_datafolder,noise.data_ft)
		
		mean_sum=0
		
		for i in range(shape(fn_noise_files)[0]):
			
			fn_load=noise.fn_datafolder+fn_noise_files[i]
			noise_vals = mpimg.imread(fn_load).astype(noise.data_enc)
			
			curr_mean=mean(noise_vals)
			mean_sum=mean_sum+curr_mean
		
		mean_final=float(mean_sum)/float(shape(fn_noise_files)[0])
	
		noise.noise=mean_final
	
	elif noise.mode=="predefined":
		pass
	
	elif noise.mode=="outside":
		
		fn_data_files=get_sorted_folder_list(noise.embryo.fn_datafolder,noise.embryo.data_ft)
		
		out_av=[]
		
		#Looping trough all data images
		for i in range(shape(fn_data_files)[0]):
		
			#Load current data img
			fn_load=noise.embryo.fn_datafolder+fn_data_files[i]
			data_img = mpimg.imread(fn_load).astype(noise.embryo.data_enc)
			data_vals=data_img.real
			data_vals=data_vals.astype('float')
			
			mask_out=(1-noise.embryo.masks_embryo[i])
			
			#Multiplying with embryo mask and substracting bkgd
			data_vals_out=data_vals*mask_out
			
			#Computing average concentrations for all 3 regions
			out_av.append(float(sum(data_vals_out))/float(sum(mask_out)))
			
		noise.noise=sum(out_av)/float(shape(out_av)[0])	
		
	return noise

#-------------------------------------------------------------------------------------------------------------------------------------
#Approximate circle properties from Mueller (2012) oval coordinates

def oval_to_circle(left,right,hor,ver):
	
	center=[left+float(hor)/2.,right+float(ver)/2.]
	radius=(float(hor)/2.+float(ver)/2.)/2.
	
	return center,radius

#-------------------------------------------------------------------------------------------------------------------------------------
#Own histogram function

def hist_own(arr,bins=None,bin_min=None,bin_max=None):
	
	#If array is not two dimensional, return nothing
	#if shape(arr)!=2:
		#raise Warning("Array has the wrong size")
		#return
	
	#Check input
	if bins==None:
		bins=10
	if bin_min==None:
		bin_min=arr.min()
	if bin_max==None:	
		bin_max=arr.max()
		
	#Make bin vector
	binvec=linspace(bin_min,bin_max,bins+1)
	binned_data=zeros(shape(binvec)[0]-1,)
	
	#Looping through image and assign to binned_data vector
	for i in range(shape(arr)[0]):
		for j in range(shape(arr)[1]):
			for k in range(bins):
				if arr[i,j]<=binvec[k+1] and binvec[k]<=arr[i,j]:
					binned_data[k]=binned_data[k]+1
					break
				
	return binned_data, binvec

#-------------------------------------------------------------------------------------------------------------------------------------
#Get pixels with maximum intensity

def get_max_int_pxs(img,maxint):
	
	max_pxs=where(img==maxint)
	return max_pxs

#-------------------------------------------------------------------------------------------------------------------------------------
#Return image with pxs labeled

def label_pxs(img,pxs,fill=nan):
	
	lab=fill*ones(shape(img))
	
	lab[pxs]=1
	
	return lab
	
	