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

def generic_data(size):##{{{
	t    = np.linspace( 0 , 1 , size )
	X    = t**2 + 0.2 * np.cos( 2 * np.pi * t )
	return t,X
##}}}

def kernel( dX , X , Y ):##{{{
	kde = sc.gaussian_kde( dX.T )
	XY = np.vstack([X.ravel(), Y.ravel()])
	Z = kde(XY).reshape(X.shape)
	return Z
##}}}


###########
## Tests ##
###########

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

def test_normal( plot = True ):##{{{
	
	print( "Test NormalLaw..." , end = "\r" )
	
	try:
		## Data
		size = 2000
		t,X,_ = sdt.Dataset.normal_cst_scale(size)
		
		## Law
		loc   = X
		scale = 0.1 + 0.1 * X
		Y = np.random.normal( loc = loc , scale = scale )
		
		## Fit
		law = sd.NormalLaw( n_bootstrap = 100 )
		law.fit( Y , loc_cov = X , scale_cov = X )
		
		
		## Fit with fix loc
		lawl = sd.NormalLaw()
		lawl.fit( Y , scale_cov = X , floc = loc )
		
		## Fit with fix scale
		laws = sd.NormalLaw()
		laws.fit( Y , loc_cov = X , fscale = scale )
		
		if plot:
			## Plot
			nrow,ncol = 2,2
			fig = plt.figure( figsize = ( 5 * ncol , 3 * nrow ) )
			fig.suptitle( "Normal Law" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." , alpha = 0.5 )
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( loc , law.loc , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( loc , laws.loc , color = "black" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax.plot( scale , law.scale , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( scale , lawl.scale , color = "black" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 4 )
			ax = sdt.plot_confidences_intervals( law , ax )
			ax.set_xticklabels( [r"$\mu_0$",r"$\mu_1$",r"$\sigma_0$",r"$\sigma_1$"] )
			
			plt.tight_layout()
			plt.show()
		print( "Test NormalLaw (Done)" )
	except:
		print( "Test NormalLaw (Fail)" )
##}}}

def test_exp( plot = True ):##{{{
	
	print( "Test ExpLaw..." , end = "\r" )
	
	try:
		## Data
		size = 2000
		t,X = generic_data(size)
		
		## Law
		scale = 0.1 + 0.3 * X
		Y = np.random.exponential( scale = scale )
		
		## Fit
		lf = sdt.InverseLinkFct()
		law = sd.ExpLaw( method = "MLE" , n_bootstrap = 100 )
		law.fit( Y , scale_cov = X )
		fscale = law.scale
		
		
		if plot:
			## Plot
			nrow,ncol = 2,2
			fig = plt.figure( figsize = ( 5 * ncol , 3 * nrow ) )
			fig.suptitle( "Exponential Law" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." , alpha = 0.5 )
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( scale , fscale , color = "blue" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax = sdt.plot_confidences_intervals( law , ax )
			ax.set_xticklabels( [r"$\sigma_0$",r"$\sigma_1$"] )
			
			
			plt.tight_layout()
			plt.show()
		print( "Test ExpLaw (Done)" )
	except:
		print( "Test ExpLaw (Fail)" )
##}}}

def test_gamma( plot = True ):##{{{
	
	print( "Test GammaLaw..." , end = "\r" )
	
	try:
		## Data
		size = 2000
		t,X = generic_data(size)
		
		## Law
		shape = 0.1 + X #np.repeat( 1.2 , size )
		scale = 0.1 + 0.1 * X
		Y = np.random.gamma( shape = shape , scale = scale )
		
		## Fit
		law = sd.GammaLaw( method = "MLE" , n_bootstrap = 100 )
		law.fit( Y , scale_cov = X , shape_cov = X )
		
		
		## Fit with fix scale
		lawsc = sd.GammaLaw()
		lawsc.fit( Y , shape_cov = X , fscale = scale )
		
		## Fit with fix shape
		lawsh = sd.GammaLaw()
		lawsh.fit( Y , scale_cov = X , fshape = shape )
		
		if plot:
			## Plot
			nrow,ncol = 2,2
			fig = plt.figure( figsize = ( 5 * ncol , 3 * nrow ) )
			fig.suptitle( "Gamma Law" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." , alpha = 0.5 )
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( scale , law.scale , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( scale , lawsh.scale , color = "black" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax.plot( shape , law.shape , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( shape , lawsc.shape , color = "black" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 4 )
			ax = sdt.plot_confidences_intervals( law , ax )
			ax.set_xticklabels( [r"$\sigma_0$",r"$\sigma_1$",r"$\xi_0$",r"$\xi_1$"] )
			
			plt.tight_layout()
			plt.show()
		print( "Test GammaLaw (Done)" )
	except:
		print( "Test GammaLaw (Fail)" )
##}}}

def test_gpd( plot = True ):##{{{
	print( "Test GPDLaw..." , end = "\r" )
	
	try:
		## Data
		size = 2000
		t,X = generic_data(size)
		
		## Law
		scale = 2 * X + 0.1
		Xs    = np.linspace( -1 , 1 , size )
		shape = -0.3 * Xs
		Y = sc.genpareto.rvs( loc = np.zeros( size ) , scale = scale , c = shape )
		
		## Fit
		gpd = sd.GPDLaw( "MLE" , n_bootstrap = 100 ) #, link_fct_shape = sdt.LogitLinkFct( -0.5 , 0.5 , 10 ) )
		gpd.fit( Y , loc = 0 , scale_cov = X , shape_cov = Xs )
		
		## Fit, with fix scale
		gpdsc = sd.GPDLaw( "MLE" ) #, link_fct_shape = sdt.LogitLinkFct( -0.5 , 0.5 , 10 ) )
		gpdsc.fit( Y , loc = 0 , fscale = scale , shape_cov = Xs )
		
		## Fit, with fix shape
		gpdsh = sd.GPDLaw( "MLE" ) #, link_fct_shape = sdt.LogitLinkFct( -0.5 , 0.5 , 10 ) )
		gpdsh.fit( Y , loc = 0 , scale_cov = X , fshape = shape )
		
		if plot:
			## Plot
			nrow,ncol = 2,2
			fig = plt.figure()
			fig.suptitle( "GPD Law" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." )
			ax.set_ylim( np.quantile( Y , [0.05,0.95] ) )
			
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( scale , gpd.scale , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( scale , gpdsh.scale , color = "blue" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax.plot( shape , gpd.shape , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( shape , gpdsc.shape , color = "blue" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 4 )
			ax = sdt.plot_confidences_intervals( gpd , ax )
			ax.set_xticklabels( [r"$\sigma_0$",r"$\sigma_1$",r"$\xi_0$",r"$\xi_1$"] )
			
			plt.tight_layout()
			plt.show()
		print( "Test GPDLaw (Done)" )
	except:
		print( "Test GPDLaw (Fail)" )
##}}}

def test_gev( plot = True ):##{{{
	print( "Test GEVLaw..." , end = "\r" )
	
	try:
		np.random.seed(0)
		
		## Data
		size = 250
		t,X = generic_data(size)
		
		## Law
		loc   = - X + 1
		scale = 0.1 + 0.1 * X
		Xs    = np.linspace( -0.5 , 0.5 , size )
		shape = Xs
		Y = sc.genextreme.rvs( loc = loc , scale = scale , c = -shape )
		
		## Fit
		law = sd.GEVLaw( method = "MLE" , n_bootstrap = 100 )# , link_fct_shape = sdt.LogitLinkFct( -0.5 , 0.5 ) )
		law.fit( Y , loc_cov = X , scale_cov = X , shape_cov = 3 * Xs + 2. )
		
		## Fit by fixing loc
		lawl = sd.GEVLaw( method = "MLE" )
		lawl.fit( Y , floc = loc , scale_cov = X , shape_cov = 3 * Xs + 2. )
		
		## Fit by fixing scale
		lawsc = sd.GEVLaw( method = "MLE" )
		lawsc.fit( Y , loc_cov = X , fscale = scale , shape_cov = 3 * Xs + 2. )
		
		## Fit by fixing shape
		lawsh = sd.GEVLaw( method = "MLE" )
		lawsh.fit( Y , loc_cov = X , scale_cov = X , fshape = shape )
		
		## Fit by fixing loc and scale
		lawlsc = sd.GEVLaw( method = "MLE" )
		lawlsc.fit( Y , floc = loc , fscale = scale , shape_cov = 3 * Xs + 2. )
		
		## Fit by fixing loc and shape
		lawlsh = sd.GEVLaw( method = "MLE" )
		lawlsh.fit( Y , floc = loc , scale_cov = X , fshape = shape )
		
		## Fit by fixing scale and shape
		lawscsh = sd.GEVLaw( method = "MLE" )
		lawscsh.fit( Y , loc_cov = X , fscale = scale , fshape = shape )
		
		if plot:
			## Plot
			nrow,ncol = 2,3
			fig = plt.figure( figsize = ( 4 * ncol , 3 * nrow ) )
			fig.suptitle( "GEV Law" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." , alpha = 0.5 )
			ax.set_ylim( np.quantile( Y , [0.05,0.95] ) )
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( loc , law.loc , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( loc , lawsc.loc , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( loc , lawsh.loc , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( loc , lawscsh.loc , color = "blue" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 3 )
			ax = sdt.plot_confidences_intervals( law , ax )
			ax.set_xticklabels( [r"$\mu_0$",r"$\mu_1$",r"$\sigma_0$",r"$\sigma_1$",r"$\xi_0$",r"$\xi_1$"] )
			
			
			ax = fig.add_subplot( nrow , ncol , 4 )
			ax.plot( scale , law.scale , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( scale , lawl.scale , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( scale , lawsh.scale , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( scale , lawlsh.scale , color = "blue" , linestyle = "" , marker = "." )
			
			ax = fig.add_subplot( nrow , ncol , 5 )
			ax.plot( shape , law.shape , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( shape , lawl.shape , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( shape , lawsc.shape , color = "blue" , linestyle = "" , marker = "." )
			ax.plot( shape , lawlsc.shape , color = "blue" , linestyle = "" , marker = "." )
			
			
			plt.tight_layout()
			plt.show()
		print( "Test GEVLaw (Done)" )
	except:
		print( "Test GEVLaw (Fail)" )
##}}}

def test_qr( plot = True ):##{{{
	
	print( "Test QuantileRegression..." , end = "\r" )
	
	try:
		## Data
		size = 2000
		t,X = generic_data(size)
		
		## Law
		loc = X
		scale = np.repeat( 0.1 , size )
		Y = np.random.normal( loc = loc , scale = scale )
		
		## Fit
		ltau = np.arange( 0.05 , 0.96 , 0.01 )
		q = sdnp.quantile( Y , ltau , X )
		
		if plot:
			## Plot
			nrow,ncol = 2,2
			fig = plt.figure()
			fig.suptitle( "Quantile Regression" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." )
			ylim = ax.get_ylim()
			
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			for i in range(ltau.size):
				ax.plot( t , q[:,i] , color = "black" , linestyle = "-" , marker = "" )
			ax.set_ylim()
			
			plt.tight_layout()
			plt.show()
		print( "Test QuantileRegression (Done)" )
	except:
		print( "Test QuantileRegression (Fail)" )
##}}}

def test_predict( plot = True ):##{{{
	print( "Test predict for NormalLaw..." , end = "\r" )
	
	try:
		size = 2000
		t,X,_ = sdt.Dataset.normal_cst_scale(size)
		
		## Law
		loc   = X + X**2
		scale = 0.1 + 0.1 * X
		Y = np.random.normal( loc = loc , scale = scale )
		
		## Fit
		law = sd.NormalLaw()
		Z = np.vstack( (X,X**2) ).T
		law.fit( Y , loc_cov = Z , scale_cov = X )
		
		new_loc   = law.predict_loc( loc_cov = np.vstack( (-X,X**2) ).T )
		new_scale = law.predict_scale( scale_cov = X**2 )
		Yn = np.random.normal( loc = new_loc , scale = new_scale )
		
		if plot:
			## Plot
			nrow,ncol = 2,1
			fig = plt.figure( figsize = ( 5 * ncol , 3 * nrow ) )
			fig.suptitle( "Normal Law" )
			
			ax = fig.add_subplot( nrow , ncol , 1 )
			ax.plot( t , X , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Y , color = "blue" , linestyle = ""  , marker = "." , alpha = 0.5 )
			
			ax = fig.add_subplot( nrow , ncol , 2 )
			ax.plot( t , new_loc , color = "red"  , linestyle = "-" , marker = ""  )
			ax.plot( t , Yn , color = "blue" , linestyle = ""  , marker = "." , alpha = 0.5 )
			plt.show()
		print( "Test predict for NormalLaw (Done)" )
	except:
		print( "Test predict for NormalLaw (Fail)" )
##}}}


def test_MultivarNormalLaw_stationary( plot = True ):##{{{
	print( "Test MultivariateNormalLaw stationary..." , end = "\r" )
	
	try:
		mean = np.array( [5,5] )
		cov  = skd.make_spd_matrix(2)
		size = 2000
		Y    = np.random.multivariate_normal( mean = mean , cov = cov , size = size )
		
		
		## Fit
		method = "moments"
		mnorm = sd.MultivariateNormalLaw( method = method )
		mnorm.fit(Y)
		Yb = np.random.multivariate_normal( mean = mnorm.mean[0,:] , cov = mnorm.cov[0,:,:] , size = size )
		cove = mnorm.cov[0,:,:]
		
		if plot:
			## Kernel
			xymin   = min(Y.min(),Yb.min())
			xymax   = max(Y.max(),Yb.max())
			XX,YY = np.mgrid[xymin:xymax:100j,xymin:xymax:100j]
			ZY  = kernel( Y  , XX , YY )
			ZYb = kernel( Yb , XX , YY )
			
			
			## Plot
			cmap0 = plt.cm.Blues
			cmap1 = plt.cm.Reds
			
			fig = plt.figure( figsize = (8,8) )
			
			ax = fig.add_subplot( 1 , 1, 1 )
			ax.contourf( XX , YY , ZY  , cmap = cmap0 )
			ax.contour(  XX , YY , ZYb , cmap = cmap1 )
			ax.plot( Y[:,0] , Y[:,1] , linestyle = "" , marker = "." , color = "black" )
			
			plt.tight_layout()
			plt.show()
		print( "Test MultivariateNormalLaw stationary (Done)" )
	except:
		print( "Test MultivariateNormalLaw stationary (Fail)" )
##}}}

def test_MultivarNormalLaw_nonstationary( plot = True ):##{{{
	print( "Test MultivariateNormalLaw non-stationary..." , end = "\r" )
	
	try:
		## Data
		size = 2000
		t    = np.linspace( 0 , 1 , size )
		X    = t**2 + 0.1
		mean = np.array( [X,-X] ).T
		cov0 = skd.make_spd_matrix(2)
		cov = np.array( [ cov0 * x for x in X ] )
		Y    = np.array( [ np.random.multivariate_normal( mean = mean[i,:] , cov = cov[i,:,] , size = 1 ) for i in range(size) ] ).squeeze()
		
		
		## Fit
		method = "moments"
		mnorm = sd.MultivariateNormalLaw( method = method )
		mnorm.fit( Y , mean_cov = X , cov_cov = X )
		
		## Generate dataset at some time step
		ndata = 4
		ltime = np.array( (size-1) * np.linspace( 0 , 1 , ndata ) , dtype = np.int )
		Yn = np.array( [ np.random.multivariate_normal( mean = mean[i,:]       , cov = cov[i,:,]       , size = size ) for i in ltime ] ).squeeze()
		Ye = np.array( [ np.random.multivariate_normal( mean = mnorm.mean[i,:] , cov = mnorm.cov[i,:,] , size = size ) for i in ltime ] ).squeeze()
		
		if plot:
			## Kernels
			Zn,Ze = [],[]
			XX = [None for _ in range(ndata)]
			YY = [None for _ in range(ndata)]
			
			
			for i in range(ndata):
				xymin = min( Yn[i,:,:].min() , Ye[i,:,:].min() )
				xymax = max( Yn[i,:,:].max() , Ye[i,:,:].max() )
				XX[i],YY[i] = np.mgrid[xymin:xymax:100j,xymin:xymax:100j]
				Zn.append( kernel( Yn[i,:,:] , XX[i] , YY[i] ) )
				Ze.append( kernel( Ye[i,:,:] , XX[i] , YY[i] ) )
			
			
			## Plot
			cmap0 = plt.cm.Blues
			cmap1 = plt.cm.inferno
			
			nrow,ncol = 2,2
			
			fig = plt.figure()
			
			for i in range(ndata):
				ax = fig.add_subplot( nrow , ncol , i + 1 , aspect = "equal" )
				ax.contourf( XX[i] , YY[i] , Zn[i] , cmap = cmap0 , extent = [xymin,xymax,xymin,xymax] )
				ax.contour(  XX[i] , YY[i] , Ze[i] , cmap = cmap1 )
				ax.plot( Yn[i,:,0] , Yn[i,:,1] , linestyle = "" , marker = "." , color = "black" , markersize = 1 )
				ax.plot( Ye[i,:,0] , Ye[i,:,1] , linestyle = "" , marker = "." , color = "red" , markersize = 1 )
			
			plt.tight_layout()
			plt.show()
		print( "Test MultivariateNormalLaw non-stationary (Done)" )
	except:
		print( "Test MultivariateNormalLaw non-stationary (Fail)" )
##}}}


def run_all_tests( plot = False ):##{{{
	test_np(plot)
	test_normal(plot)
	test_exp(plot)
	test_gamma(plot)
	test_gpd(plot)
	test_gev(plot)
	test_qr(plot)
	test_predict(plot)
	test_MultivarNormalLaw_stationary(plot)
	test_MultivarNormalLaw_nonstationary(plot)
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


