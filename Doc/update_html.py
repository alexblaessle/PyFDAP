#Since htlatex lacks some output options, this Python script fixes everything

import os
cwd=os.getcwd()
files_before=os.listdir(cwd)

os.system("htlatex manual.tex")

try:
	os.mkdir("html2")
except OSError: 
	pass

files_after=os.listdir(cwd)

for fn in files_after:
	if fn not in files_before:
		os.system("mv " + fn +" html2/")
		
print "Done."		