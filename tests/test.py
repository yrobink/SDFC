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

def test_np( plot = True ):##{{{
	print( "Test non-parametric..." , end = "\r" )
	try:
		## Data
		size = 2000
		t,X = generic_data(size)
		
		loc   = X
		scale = 0.1 * X + 0.1
		Y = np.random.normal( loc = X , scale = scale )
		
		## Stats
		m   = sdnp.mean( Y , X )
		s   = sdnp.std( Y , X , m = m , linkFct = sdt.ExpLinkFct() )
		med = sdnp.median( Y , X )
		
		if plot:
			## Plot
			nrow,ncol = 2,2
			fig = plt.figure()
			fig.suptitle( "Non parametric" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." )
			ylim = ax.get_ylim()
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( t , m     , color = "red"    , linestyle = "-"  , marker = "" , label = "mean" )
			ax.plot( t , m - s , color = "red"    , linestyle = "--" , marker = "" , label = "std" )
			ax.plot( t , m + s , color = "red"    , linestyle = "--" , marker = "" , label = "std" )
			ax.plot( t , med   , color = "green"  , linestyle = "-"  , marker = "" , label = "median" )
			ax.set_ylim(ylim)
			ax.legend( loc = "upper left" )
			
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax.plot( loc , m , color = "blue" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 4 )
			ax.plot( scale , s , color = "blue" , linestyle = "" , marker = "." )
			
			plt.tight_layout()
			plt.show()
		print( "Test non-parametric (Done)" )
	except:
		print( "Test non-parametric (Fail)" )
##}}}


## Multivariate tests
##===================

def test_MultivariateNormalLaw( size = 2500 , plot = False ):##{{{
	
	## Just function for plot
	def kernel( dX , X , Y ):##{{{
		kde = sc.gaussian_kde( dX.T )
		XY = np.vstack([X.ravel(), Y.ravel()])
		Z = kde(XY).reshape(X.shape)
		return Z
	##}}}
	
	print( "Test of MultivariateNormalLaw" )
	
	## Stationary data
	mean_s = np.array( [5,5] )
	cov_s  = skd.make_spd_matrix(2)
	Y_s    = np.random.multivariate_normal( mean = mean_s , cov = cov_s , size = size )
	
	
	## Fit stationary
	try:
		method = "moments"
		mnorm_s = sd.MultivariateNormalLaw( method = method , n_bootstrap = 10 )
		mnorm_s.fit(Y_s)
		print( "......OK   (Stationary fit)" )
		Ye_s = np.random.multivariate_normal( mean = mnorm_s.mean[0,:] , cov = mnorm_s.cov[0,:,:] , size = size )
	except:
		print( "......Fail (Stationary fit)" )
	
	## Test stationary fit with fix parameters
	try:
		fmnorm_s = sd.MultivariateNormalLaw( method = method )
		fmnorm_s.fit( Y_s , fmean = mean_s )
		fmnorm_s.fit( Y_s , fcov  = cov_s  )
		print( "......OK   (Stationary fit with fix values)" )
	except:
		print( "......Fail (Stationary fit with fix values)" )
	
	## Test stationary predict
	try:
		mnorm_s.predict_mean()
		mnorm_s.predict_cov()
		print( "......OK   (Stationary predict)" )
	except:
		pass
		print( "......FAIL (Stationary predict)" )
	
	## Test stationary bootstrap law
	try:
		mnorm_s.bootstrap_law(0)
		print( "......OK   (Stationary bootstrap law)" )
	except:
		print( "......FAIL (Stationary bootstrap law)" )
	
	
	## Non stationary data
	t       = np.linspace( 0 , 1 , size )
	X       = t**2 + 0.1
	mean_ns = np.array( [X,-X] ).T
	cov0    = skd.make_spd_matrix(2)
	cov_ns  = np.array( [ cov0 * x for x in X ] )
	Y_ns    = np.array( [ np.random.multivariate_normal( mean = mean_ns[i,:] , cov = cov_ns[i,:,] , size = 1 ) for i in range(size) ] ).squeeze()
	
	## Fit non-stationary
	try:
		method = "moments"
		mnorm_ns = sd.MultivariateNormalLaw( method = method , n_bootstrap = 10 )
		mnorm_ns.fit( Y_ns , mean_cov = X , cov_cov = X )
		print( "......OK   (Non-stationary fit)" )
		ndata = 4
		ltime = np.array( (size-1) * np.linspace( 0 , 1 , ndata ) , dtype = np.int )
		Yn_ns = np.array( [ np.random.multivariate_normal( mean = mean_ns[i,:]       , cov = cov_ns[i,:,]       , size = size ) for i in ltime ] ).squeeze()
		Ye_ns = np.array( [ np.random.multivariate_normal( mean = mnorm_ns.mean[i,:] , cov = mnorm_ns.cov[i,:,] , size = size ) for i in ltime ] ).squeeze()
	except:
		print( "......Fail (Non-stationary fit)" )
	
	## Test non-stationary fit with fix parameters
	try:
		fmnorm_ns = sd.MultivariateNormalLaw( method = method )
		fmnorm_ns.fit( Y_ns , fmean = mean_ns , cov_cov = X)
		fmnorm_ns.fit( Y_ns , mean_cov = X , fcov  = cov_ns  )
		print( "......OK   (Non-stationary fit with fix values)" )
	except:
		print( "......Fail (Non-stationary fit with fix values)" )
	
	## Test non-stationary predict
	try:
		mnorm_ns.predict_mean( X**2 )
		mnorm_ns.predict_cov( X**2 )
		print( "......OK   (Non-stationary predict)" )
	except:
		pass
		print( "......FAIL (Non-stationary predict)" )
	
	## Test non-stationary bootstrap law
	try:
		mnorm_ns.bootstrap_law(0)
		print( "......OK   (Non-stationary bootstrap law)" )
	except:
		print( "......FAIL (Non-stationary bootstrap law)" )
	
	
	if plot:
		## Plot stationary
		xymin   = min(Y_s.min(),Ye_s.min())
		xymax   = max(Y_s.max(),Ye_s.max())
		XX,YY   = np.mgrid[xymin:xymax:100j,xymin:xymax:100j]
		ZY_s    = kernel( Y_s  , XX , YY )
		ZYe_s   = kernel( Ye_s , XX , YY )
		
		
		## Plot
		cmap0 = plt.cm.Blues
		cmap1 = plt.cm.Reds
		
		fig = plt.figure( figsize = (8,8) )
		
		ax = fig.add_subplot( 1 , 1, 1 )
		ax.contourf( XX , YY , ZY_s  , cmap = cmap0 )
		ax.contour(  XX , YY , ZYe_s , cmap = cmap1 )
		ax.plot( Y_s[:,0] , Y_s[:,1] , linestyle = "" , marker = "." , color = "black" )
		
		plt.tight_layout()
		plt.show()
		
		## Plot non-stationary
		Zn,Ze = [],[]
		XX = [None for _ in range(ndata)]
		YY = [None for _ in range(ndata)]
		
		
		for i in range(ndata):
			xymin = min( Yn_ns[i,:,:].min() , Ye_ns[i,:,:].min() )
			xymax = max( Yn_ns[i,:,:].max() , Ye_ns[i,:,:].max() )
			XX[i],YY[i] = np.mgrid[xymin:xymax:100j,xymin:xymax:100j]
			Zn.append( kernel( Yn_ns[i,:,:] , XX[i] , YY[i] ) )
			Ze.append( kernel( Ye_ns[i,:,:] , XX[i] , YY[i] ) )
		
		
		## Plot
		cmap0 = plt.cm.Blues
		cmap1 = plt.cm.inferno
		
		nrow,ncol = 2,2
		
		fig = plt.figure()
		
		for i in range(ndata):
			ax = fig.add_subplot( nrow , ncol , i + 1 , aspect = "equal" )
			ax.contourf( XX[i] , YY[i] , Zn[i] , cmap = cmap0 , extent = [xymin,xymax,xymin,xymax] )
			ax.contour(  XX[i] , YY[i] , Ze[i] , cmap = cmap1 )
			ax.plot( Yn_ns[i,:,0] , Yn_ns[i,:,1] , linestyle = "" , marker = "." , color = "black" , markersize = 1 )
			ax.plot( Ye_ns[i,:,0] , Ye_ns[i,:,1] , linestyle = "" , marker = "." , color = "red" , markersize = 1 )
		
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
#	test_np(plot)
	
	## Test multivariate
	test_MultivariateNormalLaw( size = size , plot = plot )
##}}}



#############
## Classes ##
#############



##########
## main ##
##########

if __name__ == "__main__":
	
	print(sd.__version__)
#	run_all_tests( plot = False )
	
	
	print("Done")


