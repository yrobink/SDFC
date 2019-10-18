# -*- coding: utf-8 -*-

#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

###############
## Libraries ##
##{{{

import sys,os
import pickle as pk
import multiprocessing as mp


## Scientific libraries
##=====================

import numpy as np
import scipy.stats as sc
import pandas as pd
import xarray as xr


## Plot libraries ##
##==================

import matplotlib as mpl
try:
	import matplotlib.pyplot as plt
except:
	mpl.use("Qt5Agg")
	import matplotlib.pyplot as plt
## from mpl_toolkits.mplot3d import Axes3D
## import cartopy.crs as ccrs
## import cartopy.feature as cfeature

#mpl.rcParams['font.size'] = 30
#plt.rc('text',usetex=True)
#plt.rcParams['text.latex.unicode'] = True

##}}}
###############

import SDFC as sd
import SDFC.NonParametric as sdn
import SDFC.tools as sdt
import texttable as tt


from SDFC.NonParametric.__lmoments import _lmoments2
from SDFC.NonParametric.__lmoments import _lmoments3

###############
## Fonctions ##
###############

def test_normal():##{{{
	Law = sd.Normal
	
	## Dataset
	size  = 2000
	t,X_loc,X_scale,_ = sdt.Dataset.covariates(size)
	loc   = 1. + 0.8 * X_loc# - 0.5 * X_loc**2
	scale = 0.08 * X_scale
	
	Y = np.random.normal( loc = loc , scale = scale )
	
	## Fit for Normal law
	law = Law( method = "MLE" , n_bootstrap = 10 )
	law.fit( Y , c_loc = X_loc , c_scale = X_scale )
	print(law)
##}}}

def test_exponential():##{{{
	
	## Dataset
	size  = 2000
	t,X_loc,X_scale,_ = sdt.Dataset.covariates(size)
	loc   = 1. + 0.8 * X_loc# - 0.5 * X_loc**2
	scale = 0.08 * X_scale
	
	Y = np.random.exponential( scale = scale )
	
	Law = sd.Exponential
	law = Law( method = "MLE" , n_bootstrap = 10 )
	law.fit( Y , c_loc = X_loc , c_scale = X_scale )
	print(law)
##}}}

def test_gamma():##{{{
	## Dataset
	size  = 2000
	t,X_loc,X_scale,X_shape = sdt.Dataset.covariates(size)
	loc   = 1. + 0.8 * X_loc# - 0.5 * X_loc**2
	scale = 0.08 * X_scale
	shape = 1. + 0.3 * X_shape
	
	Y = np.random.gamma( scale = scale , shape = shape )
	
	Law = sd.Gamma
	law = Law( method = "MLE" , n_bootstrap = 10 )
	law.fit( Y , f_scale = scale , c_shape = X_shape )
	print(law)
##}}}


##########
## main ##
##########

if __name__ == "__main__":
	
#	test_normal()
#	test_exponential()
#	test_gamma()
	
	X = sc.genextreme.rvs( size = 10000 , loc = 1 , scale = 0.5 , c = 0.3 )
	
	print(sdn.lmoments(X))
	print(np.mean(X))
	print(_lmoments2(X))
	print(_lmoments3(X))
	
	print("Done")
