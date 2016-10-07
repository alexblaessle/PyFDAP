import os

#Settings
version="1.0"
target_dir="/home/alex_loc/Documents/Research/PyFDAP/pyfdap_v"+version

#Making necessary folders 
try:
	os.mkdir(target_dir+"/Doc")
except OSError:
	pass
try:
	os.mkdir(target_dir+"/TestDataset")
except OSError:
	pass
try:
	os.mkdir(target_dir+"/TestDataset/results")
except OSError:
	pass
	
file_list=["pyfdap_app.py","pyfdap_img_module.py","pyfdap_misc_module.py","pyfdap_fit_module.py","pyfdap_misc_module.py",
	   "pyfdap_stats_module.py","embryo.py","molecule.py","pyfdap_subwin.py","pyfdap_term.py","pyfdap_conf.py","Doc/manual.pdf","TestDataset/results/TestDataset_20min.pk"]

folder_list=["TestDataset/squint-background_20min-interval","TestDataset/squint-dendra2_20min-interval"]

for fn in file_list:
	fn_to=""
	if "/" in fn:
		fn_to=fn.split("/")
		fn_to="/".join(fn_to[:-1])
		
	os.system("cp -v "+ " " + fn + " " + target_dir +"/"+fn_to )

for fn in folder_list:
	if "/" in fn:
		fn_to=fn.split("/")
		fn_to="/".join(fn_to[:-1])
	
	os.system("cp -vr "+ " " +fn + " " + target_dir  +"/"+fn_to )








