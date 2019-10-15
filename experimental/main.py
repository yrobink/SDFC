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

import SDFC.tools as sdt
import scipy.optimize as sco
import scipy.linalg as scl

from SDFC.tools.__LinkFct      import IdLinkFct


###############
## Fonctions ##
###############


def mean( Y , X = None , linkFct = IdLinkFct() , return_coef = False ):##{{{
	"""
		SDFC.NonParametric.mean
		=======================
		
		Estimate the mean
		
		Parameters
		----------
		Y       : np.array
			Dataset to fit the mean
		X       : np.array or None
			Covariate(s)
		linkFct : class based on SDFC.tools.LinkFct
			Link function, default is identity
		return_coef : bool
			If true, return coefficients with covariates, else return mean fitted
		
		Returns
		-------
		The mean or the coefficients of the regression
	"""
	out,coef = None,None
	if X is None:
		out = np.mean(Y)
		coef = linkFct.inverse(out)
	else:
		size = X.shape[0]
		if X.ndim == 1:
			X = X.reshape( (size,1) )
		design = np.hstack( ( np.ones( (size,1) ) , X ) )
		coef,_,_,_ = scl.lstsq( design , linkFct.inverse(Y) )
		out = linkFct( design @ coef )
	return coef if return_coef else out
##}}}

def var( Y , X = None , m = None , linkFct = IdLinkFct() , return_coef = False ):##{{{
	"""
		SDFC.NonParametric.var
		======================
		
		Estimate variance given a covariate (or not)
		
		Parameters
		----------
		Y       : np.array
			Dataset to fit the variance
		X       : np.array or None
			Covariate(s)
		m       : np.array or float or None
			mean of Y. If None, m = np.mean(Y)
		linkFct : class based on SDFC.tools.LinkFct
			Link function, default is identity
		return_coef : bool
			If true, return coefficients with covariates, else return variance fitted
		
		Returns
		-------
		The variance
	"""
	out,coef = None,None
	if X is None:
		out  = np.var(Y)
		coef = linkFct.inverse(out)
	else:
		m = np.mean( Y , axis = 0 ) if m is None else np.array( [m] ).reshape(-1,1)
		Yres = ( Y - m )**2
		size = X.shape[0]
		if X.ndim == 1:
			X = X.reshape(-1,1)
		design = np.hstack( ( np.ones( (size,1) ) , X ) )
		coef,_,_,_ = scl.lstsq( design , linkFct.inverse( Yres ) )
		out = np.abs( linkFct( design @ coef ) )
	
	return coef if return_coef else out
##}}}

def std( Y , X = None , m = None , linkFct = IdLinkFct() , return_coef = False ):##{{{
	"""
		SDFC.NonParametric.std
		======================
		
		Estimate standard deviation given a covariate (or not)
		
		Parameters
		----------
		Y       : np.array
			Dataset to fit the standard deviation
		X       : np.array or None
			Covariate(s)
		m       : np.array or float or None
			mean of Y. If None, m = np.mean(Y)
		linkFct : class based on SDFC.tools.LinkFct
			Link function, default is identity
		return_coef : bool
			If true, return coefficients with covariates, else return standard deviation fitted
		
		Returns
		-------
		The standard deviation
	"""
	out = np.sqrt( var( Y , X , m , linkFct ) )
	if return_coef:
		if X is None:
			coef = linkFct.inverse(out)
		else:
			size = X.shape[0]
			if X.ndim == 1:
				X = X.reshape( (size,1) )
			design = np.hstack( ( np.ones( (size,1) ) , X ) )
			coef,_,_,_ = scl.lstsq( design , linkFct.inverse( out ) )
			out = linkFct( design @ coef )
		return coef
	
	return out
##}}}

#############
## Classes ##
#############



class LawParam:##{{{
	"""
	Internal class, do not use
	"""
	def __init__( self , linkFct = IdLinkFct() , kind = "Parameter" ):##{{{
		self.linkFct    = linkFct
		self.coef_      = None
		self.design_    = None
		self.size       = 0
		self.dim        = 0
		self.shape      = None
		self._transform = None
		self._value     = None
		self._not_fixed = None
		self.kind       = kind
	##}}}
	
	def copy(self):##{{{
		p            = LawParam( self.linkFct , self.kind )
		p.coef_      = np.copy(self.coef_)   if self.coef_   is not None else None
		p.design_    = np.copy(self.design_) if self.design_ is not None else None
		p.size       = self.size
		p.dim        = self.dim
		p.shape      = self.shape
		p._transform = self._transform
		p._value     = np.copy(self._value)  if self._value  is not None else None
		p._not_fixed = self._not_fixed
		return p
	##}}}
	
	def init( self , X = None , fix_values = None , size = None , dim = 1 , transform = lambda x : x ): ##{{{
		self._not_fixed = fix_values is None
		self.dim        = dim
		self._transform = transform
		if fix_values is not None:
			if dim == 1:
				fix_values = np.array( [fix_values] ).reshape(-1,dim)
				if fix_values.size == 1:
					fix_values = np.repeat( fix_values[0,0] , X.size ).reshape(-1,dim)
			self._value = self.linkFct.inverse(fix_values)
		else:
			if X is None and size is not None:
				self.coef_   = np.zeros( (1,dim) )
				self.design_ = np.ones( (size,1) )
			elif X is not None:
				size = X.shape[0]
				if X.ndim == 1:
					X = X.reshape( (size,1) )
				self.design_ = np.hstack( ( np.ones( (size,1) ) , X ) )
				self.coef_   = np.zeros( (X.shape[1]+1,dim) )
			else:
				print( "Error" )
			
			if np.linalg.matrix_rank(self.design_) < self.design_.shape[1]:
				print( "SDFC.LawParam: singular design matrix" )
				self.coef_   = np.zeros( (1,dim) )
				self.design = np.ones( (size,1) )
			self.size  = np.prod(self.coef_.shape)
			self.shape = self.coef_.shape
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.LawParam\n"
		val += "-------------\n"
		val += "* kind     : {}\n".format(self.kind)
		val += "* fixed    : {}\n".format(not self._not_fixed)
		val += "* link_fct : {}\n".format(str(self.linkFct))
		val += "* coef_    : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def not_fixed(self):##{{{
		return self._not_fixed
	##}}}
	
	def set_coef( self , coef ):##{{{
		if self._not_fixed:
			self.coef_ = coef.reshape( self.shape )
	##}}}
	
	def set_intercept( self , intercept ):##{{{
		if self._not_fixed:
			self.coef_[0,:] = intercept
	##}}}
	
	def design_wo1( self ):##{{{
		"""
		Return design matrix without intercept
		"""
		return None if self.shape[0] == 1 else self.design_[:,1:]
	##}}}
	
	def update(self):##{{{
		if self._not_fixed:
			self._value = self._transform( self.design_ @ self.coef_ )
	##}}}
	
	def value( self ):##{{{
		return self._value
	##}}}
	
	def valueLf( self ):##{{{
		out = self.linkFct( self._value )
		return out
	##}}}
	
	def valueGrLf( self ):##{{{
		return self.linkFct.gradient( self._value )
	##}}}

##}}}


########################
## New Params classes ##
########################

class AbstractParam:##{{{
	def __init__( self , kind , n_samples , **kwargs ):
		self.kind       = kind
		self.link    = IdLinkFct() if kwargs.get("l_" + self.kind) is None else kwargs.get("l_" + self.kind)
		self.n_samples  = n_samples
		self.n_features = 0
		self.coef_      = None
		self.fit_       = None
	
	def is_fix(self):
		pass
	
	def value( self ):
		return self.link(self.fit_)
	
	def gradient( self ):
		return self.link.gradient(self.fit_)
	
##}}}

class CovariateParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.design_ = None
		
		## Build design matrix
		X = kwargs.get( "c_" + self.kind )
		if X.ndim == 1: X = X.reshape(-1,1)
		self.n_features = X.shape[1] + 1
		self.design_ = np.hstack( ( np.ones((self.n_samples,1)) , X ) )
		self.coef_   = np.zeros(self.n_features)
		
		if np.linalg.matrix_rank(self.design_) < self.n_features:
			print( "SDFC.LawParam: singular design matrix" )
			self.coef_   = np.zeros(1)
			self.design_ = np.ones((self.n_samples,1))
		
		self.update()
	
	def is_fix(self):
		return False
	
	def update( self ):
		self.fit_ = (self.design_ @ self.coef_).reshape(-1,1)
	
	def set_coef( self , coef ):
		self.coef_ = coef.squeeze()
		self.update()
	
	def design_wo1(self):
		return self.design_[:,1:]
##}}}

class StationaryParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.coef_      = np.zeros(1)
		self.n_features = 1
		self.fit_       = np.zeros(self.n_samples)
	
	def is_fix(self):
		return False
	
	def update( self ):
		self.fit_ = np.repeat( self.coef_ , self.n_samples ).squeeze()
	
	def set_coef( self , coef ):
		self.coef_ = np.array([coef]).squeeze()
		self.update()
	
	def design_wo1(self):
		return None
##}}}

class FixParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.fit_       = np.array( [self.link.inverse( kwargs.get( "f_" + self.kind ) )] ).squeeze()
		if self.fit_.size == 1: self.fit_ = np.repeat( self.fit_ , self.n_samples ).squeeze()
	
	def is_fix(self):
		return True
	
##}}}

class LawParams:##{{{
	def __init__( self , kinds , **kwargs ):
		self.kinds = kinds
		self._dparams = {}
	
	def add_params( self , n_samples , **kwargs ):
		for kind in self.kinds:
			k_param = self.filter( kind , **kwargs )
			config = self.infer_configuration(**k_param)
			if self.is_covariate(config):  self._dparams[kind] = CovariateParam(  kind , n_samples , **k_param )
			if self.is_stationary(config): self._dparams[kind] = StationaryParam( kind , n_samples , **k_param )
			if self.is_fix(config):        self._dparams[kind] = FixParam(        kind , n_samples , **k_param )
	
	def merge_coef( self ):
		coef = np.array([])
		for k in self._dparams:
			if not self._dparams[k].is_fix():
				coef = np.hstack( (coef,self._dparams[k].coef_.squeeze()) )
		return coef

	def split_coef( self , coef ):
		tcoef = []
		a,b = 0,0
		for k in self._dparams:
			if self._dparams[k].is_fix():
				tcoef.append(None)
			else:
				b = a + self._dparams[k].n_features
				tcoef.append( coef[a:b] )
				a = b
		
		return tcoef
	
	def update_coef( self , coef ):
		lcoef = self.split_coef(coef)
		i = 0
		for k in self._dparams:
			if not self._dparams[k].is_fix():
				self._dparams[k].set_coef(lcoef[i])
				i += 1
	
	def infer_configuration( self , **kwargs ):
		has_c = False
		has_s = False
		has_f = False
		has_l = False
		for k in kwargs:
			if k[:2] == "c_": has_c = True
			if k[:2] == "f_": has_f = True
			if k[:2] == "l_": has_l = True
		
		if len(kwargs) == 0 or ( len(kwargs) == 1 and has_l):
			has_s = True
		return has_c,has_s,has_f,has_l
	
	def is_covariate( self , c ):
		return c[0] and not ( c[1] or c[2] )
	
	def is_stationary( self , c ):
		return c[1] and not ( c[0] or c[2] )
	
	def is_fix( self , c ):
		return c[2] and not ( c[0] or c[1] )
	
	def filter( self , kind , **kwargs ):
		out = {}
		for k in kwargs:
			if kind in k:
				out[k] = kwargs[k]
		return out
##}}}


##############
## Base law ##
##############

class AbstractLaw:##{{{
	"""
	SDFC.AbstractLaw
	================
	
	Base class of SDFC laws
	
	Attributes
	----------
	method : string
		method used to fit
	coef_  : numpy.ndarray
		Coefficients fitted
	coefs_bootstrap: numpy.ndarray
		coef_ for each bootstrap
	confidence interval: numpy.ndarray[ shape = (2,coef_.size) ]
		Confidence interval, first line is the alpha/2 quantile, and second line the 1 - alpha/2 quantile
	alpha          : float
		Level of confidence interval
	"""
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ):##{{{
		"""
		Initialization of AbstractLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		
		"""
		self.method    = method
		self.params    = {}
		self._lparams  = None
		self.coef_     = None
		
		self._Y        = None
		self._size     = None
		
		self.n_bootstrap         = n_bootstrap
		self.coefs_bootstrap     = None
		self.confidence_interval = None
		self.alpha               = alpha
		
		self._debug = []
		
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.AbstractLaw\n"
		val += "----------------\n"
		val += "Base class, if your read this message, you have done a mistake"
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def _to_str( self ):##{{{
		tab = tt.Texttable( max_width = 0 )
		
		## Header
		header = [ str(type(self)).split(".")[-1][:-2] + " ({})".format(self.method) , "Link function" , "Is fix" , "coef" ]
		if self.confidence_interval is not None:
			header += [ "Quantile {}".format(self.alpha/2) , "Quantile {}".format( 1 - self.alpha / 2 ) ]
		tab.header( header )
		
		## Loop on params
		a = 0
		for p in self._lparams:
			coef = None if p.coef_ is None else str(p.coef_.squeeze().round(3).tolist())
			row = [ p.kind , str(p.linkFct).split(".")[-1] , str(not p._not_fixed) , coef ]
			if self.confidence_interval is not None:
				if p._not_fixed:
					b = a + p.size
					row += [ self.confidence_interval[i,a:b].squeeze().round(3).tolist() for i in range(2) ]
					a = b
				else:
					row += [ str(None) , str(None) ]
			tab.add_row( row )
		return tab.draw() + "\n"
	##}}}
	
	def _predict_param( self , param , cov = None ): ##{{{
		if cov is None:
			return param.valueLf()
		if param.size == 1:
			return np.repeat( param.valueLf()[0] , cov.size )
		if cov.ndim == 1:
			cov = cov.reshape( (cov.size,1) )
		return param.linkFct( param.coef_[0] + np.dot( cov , param.coef_[1:] ) )
	##}}}
	
	
	def _update_fit( self , coef ):
		pass

	def _negloglikelihood(self):
		pass
	
	def _optim_function( self , coef ):##{{{
		self._update_fit(coef)
		return self._negloglikelihood()
	##}}}
	
##}}}


##############
## New laws ##
##############

class NormalLaw(AbstractLaw):##{{{
	"""
	SDFC.NormalLaw
	==============
	
	Fit parameters of a Normal law, possibly with co-variable
	
	Attributes
	----------
	method : string
		method used to fit
	loc    : numpy.ndarray
		Location fitted
	scale  : numpy.ndarray
		Scale fitted
	coef_  : numpy.ndarray
		Coefficients fitted
	n_bootstrap: integer
		Numbers of bootstrap for confidence interval
	coefs_bootstrap: numpy.ndarray
		coef_ for each bootstrap
	confidence interval: numpy.ndarray[ shape = (2,coef_.size) ]
		Confidence interval, first line is the alpha/2 quantile, and second line the 1 - alpha/2 quantile
	alpha          : float
		Level of confidence interval
	"""
	
	def __init__( self , method = "MLE" , link_fct_loc = IdLinkFct() , link_fct_scale = IdLinkFct() , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of NormalLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		link_fct_loc   : a class herited from SDFC.tools.LinkFct
			Link function for loc, default is IdLinkFct()
		link_fct_scale : a class herited from SDFC.tools.LinkFct
			Link function for scale, default is SDFC.tools.IdLinkFct(). Interesting option is SDFC.tools.ExpLinkFct().
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , method , n_bootstrap , alpha )
		
		self.loc       = None
		self.scale     = None
		
		self._loc   = LawParam( linkFct = link_fct_loc   , kind = "loc"   )
		self._scale = LawParam( linkFct = link_fct_scale , kind = "scale" )
		self._lparams = [self._loc,self._scale]
	##}}}
	
	def __str__(self):##{{{
		return self._to_str()
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	
	def fit( self , Y , **kwargs ): #loc_cov = None , scale_cov = None , floc = None , fscale = None ): ##{{{
		"""
		Fit function for NormalLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		loc_cov   : None or numpy.ndarray
			Co-variates of loc in columns.
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		floc      : None or numpy.ndarray
			If not None, fix the value of loc parameter (so not fitted)
		fscale    : None or numpy.ndarray
			If not None, fix the value of scale parameter (so not fitted)
		"""
		self._size = Y.size
		
		######################
		lkinds = ["loc","scale"]
		self.params = LawParams( kinds = ["loc","scale"] )
		self.params.add_params( n_samples = Y.size , **kwargs )
		self._Y = Y.reshape(-1,1)
		self.method = "MLE"
		self.__fit()
		print(self.params._dparams["loc"].coef_)
		print(self.params._dparams["scale"].coef_)
		self.coef_ = self.params.merge_coef()
		######################
		
		
#		if self.n_bootstrap > 0:
#			self.coefs_bootstrap = []
#			if loc_cov is not None and loc_cov.ndim == 1:
#				loc_cov = loc_cov.reshape( (loc_cov.size,1) )
#			if scale_cov is not None and scale_cov.ndim == 1:
#				scale_cov = scale_cov.reshape( (scale_cov.size,1) )
#			if np.isscalar(floc):
#				floc = np.array([floc]).ravel()
#			if np.isscalar(fscale):
#				fscale = np.array([fscale]).ravel()
#			
#			
#			for i in range(self.n_bootstrap):
#				idx = np.random.choice( self._size , self._size )
#				Y_bs         = Y[idx]
#				loc_cov_bs   = None   if loc_cov   is None else loc_cov[idx,:]
#				scale_cov_bs = None   if scale_cov is None else scale_cov[idx,:]
#				floc_bs      = floc   if floc      is None or floc.size == 1   else floc[idx]
#				fscale_bs    = fscale if fscale    is None or fscale.size == 1 else fscale[idx]
#				
#				self._fit( Y_bs , loc_cov_bs , scale_cov_bs , floc_bs , fscale_bs )
#				self.coefs_bootstrap.append( self.coef_ )
#			
#			self.coefs_bootstrap = np.array( self.coefs_bootstrap )
#			self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
#		
#		self._fit( Y , loc_cov , scale_cov , floc , fscale )
		
	##}}}
	
	
	def __fit_moments(self):##{{{
		ploc   = self.params._dparams["loc"]
		pscale = self.params._dparams["scale"]
		
		## Fit loc
		if not ploc.is_fix():
			ploc.set_coef( mean( self._Y , ploc.design_wo1() , return_coef = True , linkFct = ploc.link ) )
		
		## Fit scale
		if not pscale.is_fix():
			pscale.set_coef( std( self._Y , pscale.design_wo1() , m = ploc.value() , return_coef = True , linkFct = pscale.link ) )
		
		self.loc   = ploc.value()
		self.scale = pscale.value()
	##}}}
	
	def __fit_mle(self):##{{{
		self.__fit_moments()
		coef_init = self.params.merge_coef()
		self.optim_result = sco.minimize( self._optim_function , coef_init , jac = self._gradient_optim_function , method = "BFGS" )
		self._update_fit( self.optim_result.x )
	##}}}
	
	def __fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self.__fit_moments()
		else:
			self.__fit_mle()
	##}}}
	
	
	def bootstrap_law( self , i ):##{{{
		"""
		Return a NormalLaw with coef from bootstrap
		
		Arguments
		---------
		i : integer
			Number of bootstrap
		
		Return
		------
		law : SDFC.NormalLaw
			A NormalLaw, None if n_bootstrap = 0
		"""
		if self.n_bootstrap == 0:
			return None
		law = NormalLaw( self.method , alpha = self.alpha )
		law._loc   = self._loc.copy()
		law._scale = self._scale.copy()
		loc,scale = self._split_param( self.coefs_bootstrap[i,:] )
		law._lparams = [law._loc,law._scale]
		law._loc.set_coef( loc )
		law._scale.set_coef( scale )
		law.coef_ = law._concat_param()
		law._update_param( law.coef_ )
		return law
	##}}}
	
	def predict_loc( self , loc_cov = None ):##{{{
		"""
		Return location parameter with a new co-variates
		
		Arguments
		---------
		loc_cov : np.array or None
			Covariate
		
		Return
		------
		loc : np.array
			Location parameters
		"""
		return self._predict_param( self._loc , loc_cov )
	##}}}
	
	def predict_scale( self , scale_cov = None ):##{{{
		"""
		Return scale parameter with a new co-variates
		
		Arguments
		---------
		scale_cov : np.array or None
			Covariate
		
		Return
		------
		scale : np.array
			Location parameters
		"""
		return self._predict_param( self._scale , scale_cov )
	##}}}
	
	
	
	def _negloglikelihood( self ): ##{{{
		scale2 = np.power( self.scale , 2 )
		return np.Inf if not np.all( self.scale > 0 ) else np.sum( np.log( scale2 ) ) / 2. + np.sum( np.power( self._Y - self.loc , 2 ) / scale2 ) / 2.
	##}}}
	
	def _update_fit( self , coef ):##{{{
		self.params.update_coef(coef)
		
		self.loc   = self.params._dparams["loc"].value()
		self.scale = self.params._dparams["scale"].value()
	##}}}
	
	def _gradient_optim_function( self , coef ): ##{{{
		self._update_fit(coef)
		grad = np.array( [] )
		Yc = self._Y - self.loc
		
		ploc = self.params._dparams["loc"]
		if not ploc.is_fix():
			grad_loc   = - ploc.design_.T @ (Yc / self.scale**2 * ploc.gradient() )
			grad = np.hstack( (grad,grad_loc.squeeze()) )
		
		pscale = self.params._dparams["scale"]
		if not pscale.is_fix():
			grad_scale = pscale.design_.T @ ( ( 1. / self.scale - Yc**2 / self.scale**3 ) * pscale.gradient() )
			grad = np.hstack( (grad,grad_scale.squeeze()) )
		return grad
	##}}}
##}}}





##########
## main ##
##########

if __name__ == "__main__":
	
	## Dataset
	size  = 2000
	t,X_loc,X_scale,_ = sdt.Dataset.covariates(size)
	loc   = 1. + 0.8 * X_loc
	scale = 0.08 * X_scale
	Y     = np.random.normal( loc = loc , scale = scale )
	
	
	## Fit
	norm = NormalLaw()
	norm.fit( Y , c_loc = X_loc , c_scale = X_scale )#, scale_cov = X_scale )
	
	print("Done")
