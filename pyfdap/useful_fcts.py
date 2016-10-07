def all_mol_hist(taus,splits,bin_width,ax,canvas,names=[],cum=False):
	
	import numpy as np
	
	#Clear axis
	ax.clear()
	canvas.draw()
	
	#Make bin vector
	minbin=min(taus)-np.mod(min(taus),bin_width)
	maxbin=max(taus)+np.mod(max(taus),bin_width)
	bins=np.arange(minbin,maxbin+bin_width,bin_width)
	print "bins will range from ", minbin, " to ", maxbin
	
	#Make splitted tau vector
	taus_splitted=[]
	for i in range(len(splits)-1):
		taus_splitted.append(taus[splits[i]:splits[i+1]])
		print len(taus_splitted[-1])
		
	#make hist
	
	
	if len(names)==0:
		
		ax.hist(taus_splitted,bins,cumulative=cum)
		
			
	else:
		ax.hist(taus_splitted,bins,label=names,cumulative=cum)
		ax.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)	
		
	canvas.draw()
	
	return True
	
def ignore_extra_frames(emb,frames):
	
	import numpy as np
	
	#Append new frames
	for f in frames:
		if f not in emb.ignored:
			emb.ignored.append(f)
	
	#Delete entries
	emb.tvec_ignored=np.delete(emb.tvec_data,emb.ignored)		
	emb.ext_av_data_ign=np.delete(emb.ext_av_data_d,emb.ignored)	
	emb.int_av_data_ign=np.delete(emb.int_av_data_d,emb.ignored)	
	emb.slice_av_data_ign=np.delete(emb.slice_av_data_d,emb.ignored)		
	
	return emb

def all_mol_taus_to_list(l,ind):
	taus=[]
	for mol in l:
		for emb in mol.embryos:
			if len(emb.fits)>ind:
				taus.append(emb.fits[ind].halflife_min)
				
	return taus	

def grab_ith_fits(mol,i):
	
	mol.sel_fits=[]
	for emb in mol.embryos:
		mol.sel_fits.append(emb.fits[i])
	
	return mol

	
	
	
	