# -*- coding: utf-8 -*-

#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

###############
## Libraries ##
###############

import sys,os
import pickle as pk
import multiprocessing as mp

import numpy as np
import sklearn.datasets as skd
import pandas as pd
import scipy.stats as sc
import scipy.optimize as sco

import matplotlib as mpl
try:
	import matplotlib.pyplot as plt
except:
	mpl.use("Qt5Agg")
	import matplotlib.pyplot as plt

#import SDFC                as sdo
import SDFC               as sd
import SDFC.tools         as sdt
import SDFC.NonParametric as sdnp


####################
## Paramètres mpl ##
####################

#mpl.rcParams['font.size'] = 30
#plt.rc('text',usetex=True)
#plt.rcParams['text.latex.unicode'] = True


###############
## Fonctions ##
###############




## Tests for SDFC parametric laws
##===============================

def test_law( law , generator , size ):##{{{
	
	name = str(type(law)).split(".")[-1][:-2]
	
	print("Test {} law         ".format(name))
	
	## Covariates
	_,X_loc,X_scale,X_shape = sdt.Dataset.covariates(size)
	loc   = 1.  + 0.8  * X_loc
	scale = 0.2 + 0.08 * X_scale
	shape = 1.  + 0.3  * X_shape
	
	## Dataset
	Y = generator( loc , scale , shape )
	
	## Full fit
	print("==> Full fit..." , end  = "\r" )
	try:
		law.fit( Y , c_loc = X_loc , c_scale = X_scale , c_shape = X_shape )
		print("==> Full fit. (OK)" , end  = "\n" )
	except:
		print("==> Full fit. (FAIL)" , end  = "\n" )
	
	## Stationary fit
	print("==> Stationary fit..." , end  = "\r" )
	try:
		law.fit( Y )
		print("==> Stationary fit. (OK)" , end  = "\n" )
	except:
		print("==> Stationary fit. (FAIL)" , end  = "\n" )
	
	lp = law.kinds_params
	## Fit by fixing one parameter
	print("==> Fit with one parameter fixed..." , end  = "\r" )
	try:
		if "scale" in lp or "shape" in lp:
			law.fit( Y , f_loc = loc   , c_scale = X_scale , c_shape = X_shape )
		if "loc" in lp or "shape" in lp:
			law.fit( Y , c_loc = X_loc , f_scale = scale   , c_shape = X_shape )
		if "loc" in lp or "scale" in lp:
			law.fit( Y , c_loc = X_loc , c_scale = X_scale , f_shape = shape   )
		print("==> Fit with one parameter fixed. (OK)" , end  = "\n" )
	except:
		print("==> Fit with one parameter fixed. (FAIL)" , end  = "\n" )
	
	## Fit by fixing two parameters
	print("==> Fit with two parameters fixed..." , end  = "\r" )
	try:
		if "loc" in lp:
			law.fit( Y , c_loc = X_loc , f_scale = scale   , f_shape = shape   )
		if "scale" in lp:
			law.fit( Y , f_loc = loc   , c_scale = X_scale , f_shape = shape   )
		if "shape" in lp:
			law.fit( Y , f_loc = loc   , f_scale = scale   , c_shape = X_shape )
		print("==> Fit with two parameters fixed. (OK)" , end  = "\n" )
	except:
		print("==> Fit with two parameters fixed. (FAIL)" , end  = "\n" )
##}}}

def test_gpd( size ):##{{{
	
	print("Test GPD law")
	
	## Covariates
	t,X_loc,X_scale,X_shape = sdt.Dataset.covariates(size)
	loc   = 1.  + 0.8  * X_loc
	scale = 0.2 + 0.08 * X_scale
	shape = 1.  + 0.3  * X_shape
	
	## Dataset
	Y = sc.genpareto.rvs( loc = loc , scale = scale , c = shape )
	law = sd.GPD( method = "mle" )
	
	## Full fit
	print("==> Full fit..." , end  = "\r" )
	try:
		law.fit( Y , f_loc = loc , c_scale = X_scale , c_shape = X_shape )
		print("==> Full fit. (OK)" , end  = "\n" )
	except:
		print("==> Full fit. (FAIL)" , end  = "\n" )
	
	## Stationary fit
	print("==> Stationary fit..." , end  = "\r" )
	try:
		law.fit( Y , f_loc = loc )
		print("==> Stationary fit. (OK)" , end  = "\n" )
	except:
		print("==> Stationary fit. (FAIL)" , end  = "\n" )
	
	## Fit by fixing one parameter
	print("==> Fit with one parameter fixed..." , end  = "\r" )
	try:
		law.fit( Y , f_loc = loc , f_scale = scale , c_shape = X_shape )
		law.fit( Y , f_loc = loc , c_scale = X_scale , f_shape = shape )
		print("==> Fit with one parameter fixed. (OK)" , end  = "\n" )
	except:
		print("==> Fit with one parameter fixed. (FAIL)" , end  = "\n" )
##}}}



## Tests for non-parametric tools
##===============================

def test_quantile_regression( size = 2500 ):##{{{
	
	print( "Test of QuantileRegression" )
	
	## Law
	t,X_loc,X_scale,_ = sdt.Dataset.covariates(size)
	loc   = 0.8 * X_loc
	scale = 0.08 * X_scale
	Y = np.random.normal( loc = loc , scale = scale )
	X = np.stack( (X_loc,X_scale) ).T
	
	## Fit
	ltau = np.arange( 0.05 , 0.96 , 0.01 )
	try:
		q = sdnp.quantile( Y , ltau , X )
		print( "......OK   (Fit)" )
	except:
		print( "......FAIL (Fit)" )
##}}}


## Run all tests in one function 
##==============================

def run_all_tests( size = 2500 ):##{{{

	## Test laws
	test_law( sd.Normal()      , lambda loc,scale,shape : np.random.normal( loc , scale )  , size )
	test_law( sd.Exponential() , lambda loc,scale,shape : np.random.exponential( scale )   , size )
	test_law( sd.Gamma()       , lambda loc,scale,shape : np.random.gamma( scale , shape ) , size )
	test_law( sd.GEV()         , lambda loc,scale,shape : sc.genextreme.rvs( loc = loc , scale = scale , c = -shape ) , size )
	
	## Test non parametric
	test_quantile_regression( size = size )
##}}}



#############
## Classes ##
#############



##########
## main ##
##########

if __name__ == "__main__":
	
	print(sd.__version__)
	run_all_tests()
	print("Done")


