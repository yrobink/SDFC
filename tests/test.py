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

def test_law( sdlaw , gen , size = 2500 , has_loc = True , has_scale = True , has_shape = True  , plot = True ): ##{{{
	
	law = sdlaw( n_bootstrap = 100 )
	name_law = str(law).split("\n")[0]
	print("Test of {}".format(str(sdlaw).split(".")[-1][:-2]) )
	
	## Generic co-variates
	t,X_loc,X_scale,X_shape = sdt.Dataset.covariates(size)
	if isinstance(law,sd.GammaLaw) : X_shape = 1 / ( 1 + np.exp( - 8 * (t-0.5) ) )
	
	## Build data
	loc   = 1. + 0.8 * X_loc
	scale = 0.08 * X_scale
	shape = 0.3 * X_shape
	Y = gen( loc = loc , scale = scale , shape = shape )
	
	## Fit all law
	try:
		if isinstance(law,sd.NormalLaw) : law.fit( Y , loc_cov = X_loc , scale_cov = X_scale )
		if isinstance(law,sd.ExpLaw)    : law.fit( Y ,                   scale_cov = X_scale )
		if isinstance(law,sd.GammaLaw)  : law.fit( Y ,                   scale_cov = X_scale , shape_cov = X_shape )
		if isinstance(law,sd.GPDLaw)    : law.fit( Y , loc = loc       , scale_cov = X_scale , shape_cov = X_shape )
		if isinstance(law,sd.GEVLaw)    : law.fit( Y , loc_cov = X_loc , scale_cov = X_scale , shape_cov = X_shape )
		print( "......OK   (Fit)" )
	except:
		print( "......FAIL (Fit)" )
	
	## Test bootstrap law
	try:
		law_bs = law.bootstrap_law(0)
		print( "......OK   (Bootstrap law)" )
	except:
		print( "......FAIL (Bootstrap law)" )
	
	## Test fit with fix values
	try:
		if isinstance(law,sd.NormalLaw):
			lawf = sdlaw( n_bootstrap = 0 )
			lawf.fit( Y , floc = loc , scale_cov = X_scale )
			lawf.fit( Y , loc_cov = X_loc , fscale = scale )
		if isinstance(law,sd.GammaLaw):
			lawf = sdlaw( n_bootstrap = 0 )
			lawf.fit( Y , fscale = scale , shape_cov = X_shape )
			lawf.fit( Y , scale_cov = X_scale , fshape = shape )
		if isinstance(law,sd.GPDLaw):
			lawf = sdlaw( n_bootstrap = 0 )
			lawf.fit( Y , loc = loc , fscale = scale , shape_cov = X_shape )
			lawf.fit( Y , loc = loc , scale_cov = X_scale , fshape = shape )
		if isinstance(law,sd.GEVLaw):
			lawf = sdlaw( n_bootstrap = 0 )
			lawf.fit( Y , floc = loc , scale_cov = X_scale , shape_cov = X_shape )
			lawf.fit( Y , loc_cov = X_loc , fscale = scale , shape_cov = X_shape )
			lawf.fit( Y , loc_cov = X_loc , scale_cov = X_scale , fshape = shape )
			lawf.fit( Y , floc = loc , fscale = scale , shape_cov = X_shape )
			lawf.fit( Y , loc_cov = X_loc , fscale = scale , fshape = shape )
			lawf.fit( Y , floc = loc , scale_cov = X_scale , fshape = shape )
		print( "......OK   (Fit with fix values)" )
	except:
		print( "......FAIL (Fit with fix values)" )
	
	## Test the predict
	##=================
	try:
		if has_loc   : law.predict_loc(   loc_cov   = 1.  +   2 * X_loc   )
		if has_scale : law.predict_scale( scale_cov = 0.1 + 0.4 * X_scale )
		if has_shape : law.predict_shape( shape_cov = 0.1 - 0.4 * X_shape )
		print( "......OK   (Predict functions)" )
	except:
		print( "......FAIL (Predict functions)" )
	
	
	## Plot
	if plot:
		
		nrow,ncol,fs = 2,3,5
		fig = plt.figure( figsize = (fs*ncol,fs*nrow) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.plot( t , Y , color = "blue" , linestyle = "" , marker = "." , alpha = 0.5 )
		ax.set_xlabel( r"$t$" )
		ax.set_ylabel( r"$Y$" )
		
		ax = fig.add_subplot( nrow , ncol , 4 )
		colors = [ c for c in plt.cm.Reds( np.linspace( 0.3 , 0.7 , 3 ) ) ]
		if has_loc   : ax.plot( t , X_loc   , color = colors[0] , linestyle = "-" , label = "loc cov"   )
		if has_scale : ax.plot( t , X_scale , color = colors[1] , linestyle = "-" , label = "scale cov" )
		if has_shape : ax.plot( t , X_shape , color = colors[2] , linestyle = "-" , label = "shape cov" )
		ax.set_xlabel( r"$t$" )
		ax.set_ylabel( "Co-variates" )
		ax.legend( loc = "lower right" )
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		sdt.plot_confidences_intervals( law , ax )
		labels = []
		if has_loc   : labels += [ r"$\mu_0$"    , r"$\mu_1$" ]
		if has_scale : labels += [ r"$\sigma_0$" , r"$\sigma_1$" ]
		if has_shape : labels += [ r"$\xi_0$" , r"$\xi_1$" ]
		ax.set_xticklabels( labels )
		
		if has_loc:
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax.plot( loc , law.loc , color = "red" , linestyle = "" , marker = "." , alpha = 0.5 )
			ax.plot( [min(loc.min(),law.loc.min()),max(loc.max(),law.loc.max())] , [min(loc.min(),law.loc.min()),max(loc.max(),law.loc.max())] , color = "black" )
			ax.set_xlabel( r"$\mu$" )
			ax.set_ylabel( r"$\hat{\mu}$" )
		
		if has_scale:
			ax = fig.add_subplot( nrow , ncol , 5 )
			ax.plot( scale , law.scale , color = "red" , linestyle = "" , marker = "." , alpha = 0.5 )
			ax.plot( [min(scale.min(),law.scale.min()),max(scale.max(),law.scale.max())] , [min(scale.min(),law.scale.min()),max(scale.max(),law.scale.max())] , color = "black" )
			ax.set_xlabel( r"$\sigma$" )
			ax.set_ylabel( r"$\hat{\sigma}$" )
		
		if has_shape:
			ax = fig.add_subplot( nrow , ncol , 6 )
			ax.plot( shape , law.shape , color = "red" , linestyle = "" , marker = "." , alpha = 0.5 )
			ax.plot( [min(shape.min(),law.shape.min()),max(shape.max(),law.shape.max())] , [min(shape.min(),law.shape.min()),max(shape.max(),law.shape.max())] , color = "black" )
			ax.set_xlabel( r"$\xi$" )
			ax.set_ylabel( r"$\hat{\xi}$" )
		
		fig.set_tight_layout(True)
		plt.show()
##}}}


## Tests for non-parametric tools
##===============================

def test_quantile_regression( size = 2500 , plot = True ):##{{{
	
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
	
	if plot:
		## Plot
		nrow,ncol = 1,2
		fig = plt.figure( figsize = (15,8) )
		fig.suptitle( "Quantile Regression" )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
		ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." )
		ylim = ax.get_ylim()
		
		
		ax = fig.add_subplot( nrow , ncol , 2 )
		for i in range(ltau.size):
			ax.plot( t , q[:,i] , color = "grey" , linestyle = "-" , marker = "" )
		ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." )
		ax.set_ylim()
		
		plt.tight_layout()
		plt.show()
##}}}


## Run all tests in one function 
##==============================

def run_all_tests( size = 2500 , plot = False ):##{{{

	## Test laws
	test_law( sd.NormalLaw , size = size , gen = lambda loc,scale,shape : sc.norm.rvs(  loc = loc      , scale = scale                   ) , has_loc = True  , has_scale = True , has_shape = False , plot = plot )
	test_law( sd.ExpLaw    , size = size , gen = lambda loc,scale,shape : sc.expon.rvs(                  scale = scale                   ) , has_loc = False , has_scale = True , has_shape = False , plot = plot )
	test_law( sd.GammaLaw  , size = size , gen = lambda loc,scale,shape : np.random.gamma(               scale = scale , shape = shape   ) , has_loc = False , has_scale = True , has_shape = True  , plot = plot )
	test_law( sd.GPDLaw    , size = size , gen = lambda loc,scale,shape : sc.genpareto.rvs( loc = loc  , scale = scale , c     = shape   ) , has_loc = False , has_scale = True , has_shape = True  , plot = plot )
	test_law( sd.GEVLaw    , size = size , gen = lambda loc,scale,shape : sc.genextreme.rvs( loc = loc , scale = scale , c     = - shape ) , has_loc = True  , has_scale = True , has_shape = True  , plot = plot )
	
	## Test non parametric
	test_quantile_regression( size = size , plot = plot )
##}}}



#############
## Classes ##
#############



##########
## main ##
##########

if __name__ == "__main__":
	
	print(sd.__version__)
	run_all_tests( plot = False )
	print("Done")


