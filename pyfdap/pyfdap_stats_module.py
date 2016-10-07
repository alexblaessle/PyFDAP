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

#Statistics module for PyFRAP toolbox, including following functions:
#(1) fit_Rsq: Calulates Rsq value of fit

#=====================================================================================================================================
#Importing necessary modules
#=====================================================================================================================================

from numpy import *
import scipy.stats

#=====================================================================================================================================
#Module Functions
#=====================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------
#Calulates Rsq value of fit
	
def fit_Rsq(squ_av_data_d,ssd_opt):
	
	#Get mean of squ_av_data_d
	mean_data=mean(squ_av_data_d)
	
	#Get derivation of data points from mean
	SStot=sum((squ_av_data_d-mean_data)**2)
	
	#Calculate RSq value
	Rsq=1-ssd_opt/SStot
	
	return Rsq

def wilcoxon_test(x,y,zero_method='wilcox', correction=False,printOut=True):
	stat,pval=scipy.stats.wilcoxon(x, y=y, zero_method=zero_method, correction=correction)
	
	if printOut:
		print "Results of Wilcoxon-Test:"
		print "p-Value: ", pval
		print "Wilcoxon-Statistics:", stat
	
	return stat,pval

def mann_whitney_test(x,y,printOut=True):
	stat,pval=scipy.stats.mannwhitneyu(x, y)
	
	if printOut:
		print "Results of Mann-Whitney-U-Test:"
		print "p-Value: ", pval
		print "U-Statistics:", stat
	
	return stat,pval

def ttest_standard(x,y):
	stat,pval=scipy.stats.ttest_ind(x, y, equal_var=True)
	
	if printOut:
		print "Results of Standard t-Test:"
		print "p-Value: ", pval
		print "Statistics:", stat
		
	return stat,pval
	
def ttest_welch(x,y):
	stat,pval=scipy.stats.ttest_ind(x, y, equal_var=False)
	
	if printOut:
		print "Results of Welch's t-Test:"
		print "p-Value: ", pval
		print "Statistics:", stat	
	
	return stat,pval
	
def sharipo_test(x):

	stat,pval=scipy.stats.sharipo(x)
	
	if printOut:
		print "Results of Sharipo-Test:"
		print "p-Value: ", pval
		print "Statistics:", stat
	
	return stat,pval


	
	
	
	