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

#Misc module for PyFDAP, including following functions:

#1) get_sorted_folder_list: Returns list of files in alphabetical order in selected folder
#2) get_folders_in_folder: Returns list of folders in folder
#3) create_mol_folder_struct: Create meaningful folder structure for molecule object 
#4) create_emb_folder_struct: Create meaningful folder structure for embryo object 
#5) create_nonexist_dir: Create directory if not already existent
#6) clear_data_structure: Clear molecule data structure
#7) clear_folder: Clear folder without deleting subfolders
#8) remove_existent: Remove file if existent
#9) check_folders_full: Check if folder is full
#10) check_folder_empty: Check if folder is empty
#11) write_csv_object: Write arbitrary object properties into csv file
#12) write_csv_embryo: Write embryo object properties into csv file
#13) write_csv_molecule: Write molecule object properties into csv file
#14) write_csv_timeseries: Write a list of time series into csv file
#15) adjust_folderpaths: Adjust folder paths in object if root folder is different


#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

from numpy import *
import os
import csv
import sys
import platform
from embryo import *
import pyfdap_img_module

#=====================================================================================================================================
#Module Functions
#=====================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------
#Get sorted file list from folder for files of type f_type

def get_sorted_folder_list(fn_folder,f_type):
	
	if fn_folder=='':
		return []
	
	if not os.path.exists(fn_folder):
		return []
	
	#Getting files in data folder
	file_list=os.listdir(fn_folder)
		
	file_list_new=[]
	
	#Going through all file names 
	for i in range(shape(file_list)[0]):
		#Check if it's the right file type
		if f_type in file_list[i]:
			file_list_new.append(file_list[i])
	
	#Sorting
	file_list_new.sort()
	file_list=file_list_new
	
	return file_list

#-------------------------------------------------------------------------------------------------------------------------------------
#Get folder list in folder

def get_folders_in_folder(fn_folder):
	
	if fn_folder=='':
		
		return []
	
	if not os.path.exists(fn_folder):
	
		return []
	
	#Getting files in data folder
	file_list=os.listdir(fn_folder)
		
	file_list_new=[]
	
	#Going through all file names 
	for i in range(shape(file_list)[0]):
		#Check if it's the right file type
		if os.path.isdir(fn_folder+"/"+file_list[i]):
			file_list_new.append(file_list[i])
	
	#Sorting
	file_list_new.sort()
	file_list=file_list_new
	
	return file_list

#-------------------------------------------------------------------------------------------------------------------------------------
#Create folder structure suitable for FDAP data structure

def create_mol_folder_struct(fn_folder,n_embryos,n_bkgds,name):
	
	fn_root=fn_folder+"/"+name
	
	create_nonexist_dir(fn_root)
	
	fn_embryos=fn_root+"/embryos"
	fn_bkgds=fn_root+"/bkgds"
	
	create_nonexist_dir(fn_embryos)
	create_nonexist_dir(fn_bkgds)
	
	for i in range(n_embryos):	
		curr_fn=fn_embryos+"/"+"emb"+str(i)
		create_nonexist_dir(curr_fn)
		create_nonexist_dir(curr_fn+"/"+"green")
		create_nonexist_dir(curr_fn+"/"+"red")
		create_nonexist_dir(curr_fn+"/"+"pre")
		create_nonexist_dir(curr_fn+"/"+"bright")
		create_nonexist_dir(curr_fn+"/"+"pre/red")
		create_nonexist_dir(curr_fn+"/"+"pre/green")
		create_nonexist_dir(curr_fn+"/"+"pre/bright")
		
	for i in range(n_bkgds):	
		curr_fn=fn_bkgds+"/"+"bkgd"+str(i)
		create_nonexist_dir(curr_fn)
		create_nonexist_dir(curr_fn+"/"+"green")
		create_nonexist_dir(curr_fn+"/"+"red")
		create_nonexist_dir(curr_fn+"/"+"pre")
		create_nonexist_dir(curr_fn+"/"+"bright")
		create_nonexist_dir(curr_fn+"/"+"pre/red")
		create_nonexist_dir(curr_fn+"/"+"pre/green")
		create_nonexist_dir(curr_fn+"/"+"pre/bright")
	
	create_nonexist_dir(fn_root+"/results")

def create_emb_folder_struct(fn_folder,name,green_str,red_str,bright_str,green_pre_str,red_pre_str,bright_pre_str,green_post_str,red_post_str,bright_post_str):
	
	curr_fn=fn_folder+"/"+name
	create_nonexist_dir(curr_fn+"/"+"green")
	create_nonexist_dir(curr_fn+"/"+"red")
	create_nonexist_dir(curr_fn+"/"+"pre")
	create_nonexist_dir(curr_fn+"/"+"bright")
	create_nonexist_dir(curr_fn+"/"+"pre/red")
	create_nonexist_dir(curr_fn+"/"+"pre/green")
	create_nonexist_dir(curr_fn+"/"+"pre/bright")
	create_nonexist_dir(curr_fn+"/"+"lsm")
		
	os.chdir(curr_fn)
	os.system("mv "+ green_str + " green/" )
	os.system("mv "+ red_str + " red/" )
	os.system("mv "+ bright_str + " bright/" )
	
	os.system("mv "+ green_post_str + " green/" )
	os.system("mv "+ red_post_str + " red/" )
	os.system("mv "+ bright_post_str + " bright/" )
	
	os.system("mv "+ green_pre_str + " pre/green/" )
	os.system("mv "+ red_pre_str + " pre/red/" )
	os.system("mv "+ bright_pre_str + " pre/bright/" )
	
	os.system("mv *.lsm lsm/" )
		
def create_nonexist_dir(directory):
	if not os.path.exists(directory):
		os.mkdir(directory)

def clear_data_structure(fn_folder,name):
	curr_fn=fn_folder+"/"+name
	
	clear_folder(curr_fn+"/"+"green")
	clear_folder(curr_fn+"/"+"red")
	clear_folder(curr_fn+"/"+"bright")
	clear_folder(curr_fn+"/"+"pre/red")
	clear_folder(curr_fn+"/"+"pre/green")
	clear_folder(curr_fn+"/"+"pre/bright")
		
def clear_folder(fn_folder):
	
	file_list=get_sorted_folder_list(fn_folder,'.tif')
	for fn in file_list:
		curr_fn=fn_folder+"/"+fn
		remove_existent(curr_fn)
	
def remove_existent(fn):
	
	if os.path.isfile(fn):
		os.remove(fn)
		print "Successfully removed ", fn
	else:
		print "Couldn't remove ", fn

def check_folders_full(fn_folder,name):
	
	curr_fn=fn_folder+"/"+name
	
	check_folder_empty(curr_fn+"/"+"green")
	check_folder_empty(curr_fn+"/"+"red")
	check_folder_empty(curr_fn+"/"+"pre")
	check_folder_empty(curr_fn+"/"+"bright")
	check_folder_empty(curr_fn+"/"+"pre/red")
	check_folder_empty(curr_fn+"/"+"pre/green")
	check_folder_empty(curr_fn+"/"+"pre/bright")

def check_folder_empty(fn):

	files=get_sorted_folder_list(fn,".tif")
	if files==[]:
		print "Folder:", fn, " empty"
		return True
	else:
		return False
	
#-------------------------------------------------------------------------------------------------------------------------------------
#Write objects to csv

def write_csv_sel_props(objs,fn,names=[]):
	
	#Grab selected properties of first object
	props=objs[0].selected_props
	props_print=list(props)
	props_print.insert(0,' ')
	
	wfile=csv.writer(open(fn,'wb'), delimiter=';')
	wfile.writerow(props_print)
	
	for i,obj in enumerate(objs):
		
		vals=[]
		if len(names)==len(objs):
			vals.append(names[i])
		else:
			if hasattr(obj,'name'):
				vals.append(obj.name)
			else:
				vals.append('')
				
		#Put values of properties in list 
		for item in props:
			vals.append(vars(obj)[str(item)])

		wfile.writerow(vals)
	
	return True	
			

def write_csv_object(wfile,obj):
	
	names=[]
	vals=[]
	
	#Going through all attributes of obj
	for item in vars(obj):
		
		#Don't print out large arrays
		if isinstance(vars(obj)[str(item)],(int,float,str)):
			names.append(item)
			vals.append(vars(obj)[str(item)])
				
		elif isinstance(vars(obj)[str(item)],list) and shape(vars(obj)[str(item)])<=2:
			names.append(item)
			vals.append(vars(obj)[str(item)])
	wfile.writerow(names)
	wfile.writerow(vals)	
	wfile.writerow([""])	
	
def write_csv_embryo(fn_save,embryo):
	
	wfile=csv.writer(open(fn_save,'wb'), delimiter=';')
	
	wfile.writerow([embryo.name])
	write_csv_object(wfile,embryo)
	
	wfile.writerow(["Pre"])
	write_csv_object(wfile,embryo.pre)
	
	wfile.writerow(["Noise"])
	write_csv_object(wfile,embryo.noise)
	
	wfile.writerow(["Fits"])
	for fit in embryo.fits:
		wfile.writerow([fit.name])
		write_csv_object(wfile,fit)
		
	print "Done writing ", embryo.name, "into", fn_save	
	
def write_csv_molecule(fn_save,molecule):
	
	wfile=csv.writer(open(fn_save,'wb'), delimiter=';')
	
	wfile.writerow([molecule.name])
	write_csv_object(wfile,molecule)
	
	wfile.writerow(["Embryos"])
	for embryo in molecule.embryos:
		wfile.writerow([embryo.name])
		write_csv_object(wfile,embryo)
		
		wfile.writerow(["Pre"])
		write_csv_object(wfile,embryo.pre)
		
		wfile.writerow(["Noise"])
		write_csv_object(wfile,embryo.noise)
		
		wfile.writerow(["Fits"])
		for fit in embryo.fits:
			wfile.writerow([fit.name])
			write_csv_object(wfile,fit)
			
	wfile.writerow(["Bkgds"])
	for bkgd in molecule.bkgds:
		wfile.writerow([bkgd.name])
		write_csv_object(wfile,bkgd)
		
		wfile.writerow(["Pre"])
		write_csv_object(wfile,bkgd.pre)
			
	print "Done writing ", molecule.name, "into", fn_save	
	
def write_csv_timeseries(fn_save,ts_list,ts_names):
	
	wfile=csv.writer(open(fn_save,'wb'), delimiter=';')
	
	for i in range(shape(ts_list)[0]):
		wfile.writerow([ts_names[i]])
		wfile.writerow(ts_list[i])
		
	print "Done writing ", ts_names, "into", fn_save		
		
def print_mem_usage(obj):
		
	memusage_name=[]
	memusage_size=[]
	
	#Going through all attributes of object
	for item in vars(obj):
		
		memusage_name.append(item)
		memusage_size.append(get_size_of(vars(obj)[str(item)]))
		
	zipped=zip(memusage_name, memusage_size)
	
	sorted_by_size=sorted(zipped, key=lambda zipped:zipped[1])
	
	for i in range(shape(memusage_name)[0]):
		print sorted_by_size[i][0], sorted_by_size[i][1]
		
	return sorted_by_size	

def get_size_of(obj):
	if isinstance(obj,list):
		s=0
		for l in obj:
			s=s+get_size_of(l)
	elif 'numpy' in str(type(obj)):
		s=obj.nbytes
	else:
		s=sys.getsizeof(obj)
	
	return s	
	

def adjust_folderpaths(obj,keystr,debug):
	cwd=os.getcwd()
	for item in vars(obj):
		if "datafolder" in item or "maskfolder" in item or "bkgdfolder" in item:
			
			if platform.system() in ["Windows"] and platform.architecture()[0]=="64bit":
				tempstr=vars(obj)[str(item)].split("/"+keystr)
				datafolder=cwd+"\\"+keystr+"\\"+tempstr[-1][1:]		
			else:
				tempstr=vars(obj)[str(item)].split("/"+keystr)
				datafolder=cwd+"/"+keystr+"/"+tempstr[-1][1:]	
		
			vars(obj)[str(item)] = datafolder 
		
			if debug==1:
				print "Adjusted file path of ", item, " to ", datafolder
				if not os.path.exists(datafolder):
					raise Warning(datafolder + "does not exist, please check if you have downloaded the test dataset properly")
				if len(get_sorted_folder_list(datafolder,".tif"))<1:
					raise Warning(datafolder + "seems not to contain any .tif images, please check if you have downloaded the test dataset properly")
	
	return obj

def equ_extend_vec(vec,n):
	
	d=vec[-1]-vec[-2]
	
	for i in range(1,n+1):
		if isinstance(vec,list):
			vec=vec.append(vec[-1]+i*d)
		else:
			vec=append(vec,i)
			
	return vec

def update_obj(obj_blank,obj,debug=False):
	
	if debug:
		print "======update_obj: updated properties======"
	
	#Going through all attributes blank object
	for item in vars(obj_blank):
		
		if not hasattr(obj,str(item)):
			setattr(obj, str(item), vars(obj_blank)[str(item)])
			
			if debug:
				print item, " = ", vars(self)[str(item)]
	
	return obj

def gen_embryos_from_csv(fn,ident=["time,ext","int","slice","pre_ext","pre_int","pre_slice"],delim=',',debug=True,mol=None,inSameCol=True,idxPre=[0,1],noiseVal=1):
	
	"""Generates a list of embryo objects from a csv sheet.
	
	csv sheet should look like
	
	+------+---------+---------+-----------+------+---------+---------+-----------+
	| time | ext     | int     | slice     | time | ext     | int     | slice     |  
	+======+=========+=========+===========+======+=========+=========+===========+
	| t1   | ext(t1) | int(t1) | slice(t1) | t1   | ext(t1) | int(t1) | slice(t1) |
	+------+---------+---------+-----------+------+---------+---------+-----------+
	| t2   | ext(t2) | int(t2) | slice(t2) | t2   | ext(t2) | int(t2) | slice(t2) |
	+------+---------+---------+-----------+------+---------+---------+-----------+
	| ...  | ...     | ...     | ...       | ...  | ...     | ...     | ...       |
	+------+---------+---------+-----------+------+---------+---------+-----------+
	
	An example csv sheet can be found on our website.
	
	.. note:: Identifiers should have same order as header.
	
	If molecule object is given, will add newly created embryos to molecule objects
	embryos list.
	
	Args:
		fn (str): Filepath to csv sheet.
		
	Keyword Args:
		ident (list): List of identifiers. 
		delim (str): csv delimiter.
		debug (bool): Print debugging messages.
		mol (pyfdap.molecule): Molecule object. 
		
	Returns:
		list: List of created embryos.
	
	"""
	
	header,data=read_csv(fn,delim=delim)
	
	emb=None
	
	embryos=[]
	
	j=0
	
	#Loop through all columns
	for i in range(data.shape[1]):
		
		#If time column, save tvec for further use
		if ident[0] in header[i]:
			if inSameCol:
				tvec=data[idxPre[1]+1:,i]
			else:
				tvec=data[:,i]
				
			emb=None
			continue
		
		#If header is second identifier, start new embryo object
		if ident[1] in header[i]:
			
			j=j+1		
			emb=embryo("genFromCSV"+str(j),"fdap")
			
			#Set noise value
			emb.noise.noise=noiseVal
			
			if mol!=None:
				mol.add_embryo(emb)
		
		
		if inSameCol:
			
			"""Prevalues and timeseries are in same column:
			
			Grab mean of valus between idxPre[0] and idxPre[1] and
			set it as prevalue. Then take all remaining values in column
			and set it as timeseries.
			"""
			
			if len(data[idxPre[0]:idxPre[1]+1,i])>1:
				preval=mean(data[idxPre[0]:idxPre[1]+1,i])
			else:
				preval=data[idxPre[0]:idxPre[1]+1,i]
				
			setattr(emb.pre,"pre_"+header[i],preval)	
				
			setattr(emb,header[i]+'_av_data_d',data[idxPre[1]+1:,i])
			emb.tvec_data=tvec
				
		else:
			
			"""Prevalues and timeseries are not in same column:
			
			If pre is in header, take mean of that column and apply as
			prevalue. Else set column as timeseries.
			"""
			
			if "pre" in header[i]:
			
				try:
					#If more than 1 pre value is given, take the mean
					if len(data[:,i])>1:
						preval=mean(data[:,i])
					else:
						preval=data[0,i]
						
					setattr(emb.pre,header[i],preval)		
				
				except:
					print "Warning, cannot set pre values."
					
			else:	
				try:
					setattr(emb,header[i]+'_av_data_d',data[:,i])
					emb.tvec_data=tvec
				except:
					print "Warning, cannot set timeseries with header" + header[i] + " . Embryo does not have timeseries."
		
		# Check for ignored timepoints
		vecs=[emb.ext_av_data_d,emb.int_av_data_d,emb.slice_av_data_d]
		lens=[len(emb.ext_av_data_d),len(emb.int_av_data_d),len(emb.slice_av_data_d)]
		
		emb.ignored=where(isnan(vecs[lens.index(max(lens))]))[0]
		emb.ext_av_data_ign=delete(emb.ext_av_data_d,emb.ignored)	
		emb.int_av_data_ign=delete(emb.int_av_data_d,emb.ignored)	
		emb.slice_av_data_ign=delete(emb.slice_av_data_d,emb.ignored)	
		emb.tvec_ignored=delete(emb.tvec_data,emb.ignored)	
		
		# Add some dummy values in centers and radiuses list 
		# (some functions in GUI use it to generate some stuff, so there better be something in)
		emb.centers_embr_px=len(emb.tvec_data)*[[256,256]]
		emb.radiuses_embr_px=len(emb.tvec_data)*[300]
		
		# Append
		embryos.append(emb)
				
	return embryos

def gen_bkgds_from_csv(fn,mol,ident=["time,ext","int","slice","pre_ext","pre_int","pre_slice"],delim=',',debug=True,inSameCol=True,idxPre=[0,1]):
	
	"""Generates a list of background objects from a csv sheet.
	
	csv sheet should look like
	
	+------+---------+---------+-----------+------+---------+---------+-----------+
	| time | ext     | int     | slice     | time | ext     | int     | slice     |  
	+======+=========+=========+===========+======+=========+=========+===========+
	| t1   | ext(t1) | int(t1) | slice(t1) | t1   | ext(t1) | int(t1) | slice(t1) |
	+------+---------+---------+-----------+------+---------+---------+-----------+
	| t2   | ext(t2) | int(t2) | slice(t2) | t2   | ext(t2) | int(t2) | slice(t2) |
	+------+---------+---------+-----------+------+---------+---------+-----------+
	| ...  | ...     | ...     | ...       | ...  | ...     | ...     | ...       |
	+------+---------+---------+-----------+------+---------+---------+-----------+
	
	An example csv sheet can be found on our website.
	
	.. note:: Identifiers should have same order as header.
	
	Args:
		fn (str): Filepath to csv sheet.
		
	Keyword Args:
		ident (list): List of identifiers. 
		delim (str): csv delimiter.
		debug (bool): Print debugging messages.
		mol (pyfdap.molecule): Molecule object. 
		
	Returns:
		list: List of created embryos.
	
	"""
	
	header,data=read_csv(fn,delim=delim)
	
	bkgds=[]
	
	j=0
	
	#Loop through all columns
	for i in range(data.shape[1]):
		
		#If time column, save tvec for further use
		if ident[0] in header[i]:
			if inSameCol:
				tvec=data[idxPre[1]+1:,i]
			else:
				tvec=data[:,i]
				
			
			continue
		
		#If header is second identifier, start new bkgd object
		if ident[1] in header[i]:
			
			j=j+1
			mol.add_bkgd(mol.bkgd_number,"genFromCSV"+str(j),"default")
			bkgd=mol.bkgds[-1]
				
		if inSameCol:
			
			"""Prevalues and timeseries are in same column:
			
			Grab mean of valus between idxPre[0] and idxPre[1] and
			set it as prevalue. Then take all remaining values in column
			and set it as timeseries.
			"""
			
			if len(data[idxPre[0]:idxPre[1]+1,i])>1:
				preval=mean(data[idxPre[0]:idxPre[1]+1,i])
			else:
				preval=data[idxPre[0]:idxPre[1]+1,i]
				
			setattr(bkgd.pre,"pre_"+header[i],preval)	
				
			setattr(bkgd,"bkgd_"+header[i]+'_vec',data[idxPre[1]+1:,i])
			setattr(bkgd,"bkgd_"+header[i]+"_av",mean(data[idxPre[1]+1:,i]))
				
		else:
			
			"""Prevalues and timeseries are not in same column:
			
			If pre is in header, take mean of that column and apply as
			prevalue. Else set column as timeseries.
			"""
			
			if "pre" in header[i]:
			
				
				
				try:
					#If more than 1 pre value is given, take the mean
					if len(data[:,i])>1:
						preval=mean(data[:,i])
					else:
						preval=data[0,i]
						
					setattr(bkgd.pre,"pre_"+header[i],preval)	
							
				except:
					print "Warning, cannot set pre values."
					
			else:	
				try:
						
					setattr(bkgd,"bkgd_"+header[i]+'_vec',data[:,i])
					setattr(bkgd,"bkgd_"+header[i]+"_av",mean(data[:,i]))
					
				except:
					print "Warning, cannot set timeseries with header" + header[i] + " . Embryo does not have timeseries."
		
		# Check for ignored timepoints
		vecs=[bkgd.ext_vec,bkgd.int_vec,bkgd.slice_vec]
		lens=[len(bkgd.ext_vec),len(bkgd.int_vec),len(bkgd.slice_vec)]
		
		bkgd.ignored=where(isnan(bkgd.ext_vec))[0]
		bkgd.ext_av_data_ign=delete(bkgd.ext_vec,bkgd.ignored)	
		bkgd.int_av_data_ign=delete(bkgd.int_vec,bkgd.ignored)	
		bkgd.slice_av_data_ign=delete(bkgd.slice_vec,bkgd.ignored)	
		
		# Add some dummy values in centers and radiuses list 
		# (some functions in GUI use it to generate some stuff, so there better be something in)
		bkgd.centers_embr_px=len(bkgd.ext_vec)*[[256,256]]
		bkgd.radiuses_embr_px=len(bkgd.ext_vec)*[300]
		
		bkgds.append(bkgd)
		
		
	#Update average bkgd values
	mol=pyfdap_img_module.compute_av_bkgd(mol)
	
	return bkgds
	
	
def read_csv(fn,delim=',',hasheader=True,toarray=True,dtype=float):		
	
	"""Reads in csv sheet and returns data and header.
	
	Args:
		fn (str): Filepath to csv sheet.
		
	Keyword Args:
		delim (str): csv delimiter.
		hasreader (bool): Read out header.
		toarray (bool): Convert data into numpy.ndarray.
		
	Returns:
		tuple: Tuple containing:
	
			* header (list): Header.
			* data (numpy.ndarray): Data.
			
	"""
	
	with open(fn,'r') as f:
		
		fcsv=csv.reader(f, delimiter=delim)
		
		data=[]
		for i,row in enumerate(fcsv):
			
			if i==0:
				if hasheader:
					header=row
				else:
					header=[]
					data.append(replace_list(row,'',nan))
			else:
				data.append(replace_list(row,'',nan))
		
		if asarray:
			data=asarray(data).astype(dtype)
			
	return header,data

def replace_list(l,vOld,vNew):
	
	"""Replaces all occurences of vOld with vNew in list.
	"""
	
	for i,x in enumerate(l):
		if x==vOld:
			l[i]=vNew
		
	return l	
	

	