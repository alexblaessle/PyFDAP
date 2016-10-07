#Since htlatex lacks some output options, this Python script fixes everything

import os
import sys

#Get output folder
fnout=sys.argv[1]

cwd=os.getcwd()
files_before=os.listdir(cwd)

os.system("htlatex manual.tex")



try:
	os.mkdir(fnout)
except OSError: 
	pass

files_after=os.listdir(cwd)

for fn in files_after:
	if fn not in files_before:
		os.system("mv " + fn +" "+fnout+"/")
		
print "Done."		