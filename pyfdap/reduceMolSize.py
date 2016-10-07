###Script to delete embryo.vals_slice property to make molecule files substantial smaller

#Importing modules
import pyfdap_misc_module 
from molecule import *
import gc

#Define folder with molecule files
fnFolder="../../FDAP/chd_tld/results/10min/"

#Get list of molecule files in folder
files=pyfdap_misc_module.get_sorted_folder_list(fnFolder,'.pk')



#Loop through files
for f in files:
	
	#Create dummy molecule
	mol=molecule('bla')
	
	#Load molecule file
	try:
		mol=mol.load_molecule(fnFolder+f)
	except EOFError:
		print "File ", fnFolder+f , "seems to be corrupted."
	
	print "Reducing filesize of molecule:  " , mol.name, " in file ", f
	
	#Loop through embryos and write None into vals_slice
	for emb in mol.embryos:
		emb.vals_slice=None
		emb.pre.masks_embryo=None
		emb.pre.masks_ext=None
		emb.pre.masks_int=None
	
	#Loop through backgrounds and write None into masks
	for bkgd in mol.bkgds:
		bkgd.pre.masks_embryo=None
		bkgd.pre.masks_ext=None
		bkgd.pre.masks_int=None
		
	
	#Save molecule again
	mol.save_molecule(fnFolder+f)
	
	#Remove mol to free up RAM
	del mol
	gc.collect()
	
	
print "Done"	
	
	