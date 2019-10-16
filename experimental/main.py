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
import texttable as tt

###############
## Fonctions ##
###############


def mean( Y , c_Y = None , link = IdLinkFct() , value = True ):##{{{
	"""
		SDFC.NonParametric.mean
		=======================
		
		Estimate the mean
		
		Parameters
		----------
		Y     : np.array
			Dataset to fit the mean
		c_Y   : np.array or None
			Covariate(s)
		link  : class based on SDFC.tools.LinkFct
			Link function, default is identity
		value : bool
			If true return value fitted, else return coefficients of fit
		
		Returns
		-------
		The mean or the coefficients of the regression
	"""
	out,coef = None,None
	if c_Y is None:
		out = np.mean(Y)
		coef = link.inverse(out)
	else:
		size = c_Y.shape[0]
		if c_Y.ndim == 1:
			c_Y = c_Y.reshape(-1,1)
		design = np.hstack( ( np.ones((Y.size,1)) , c_Y ) )
		coef,_,_,_ = scl.lstsq( design , link.inverse(Y) )
		out = link( design @ coef )
	return out if value else coef
##}}}

def var( Y , c_Y = None , m_Y = None , link = IdLinkFct() , value = True ):##{{{
	"""
		SDFC.NonParametric.var
		======================
		
		Estimate variance given a covariate (or not)
		
		Parameters
		----------
		Y     : np.array
			Dataset to fit the mean
		c_Y   : np.array or None
			Covariate(s)
		m_Y   : np.array or float or None
			mean of Y. If None, m = np.mean(Y)
		link  : class based on SDFC.tools.LinkFct
			Link function, default is identity
		value : bool
			If true return value fitted, else return coefficients of fit
		
		Returns
		-------
		The variance
	"""
	out,coef = None,None
	if c_Y is None:
		out  = np.var(Y)
		coef = link.inverse(out)
	else:
		m_Y = np.mean( Y , axis = 0 ) if m_Y is None else np.array( [m_Y] ).reshape(-1,1)
		Yres = ( Y - m_Y )**2
		if c_Y.ndim == 1: c_Y = c_Y.reshape(-1,1)
		design = np.hstack( ( np.ones((Y.size,1)) , c_Y ) )
		coef,_,_,_ = scl.lstsq( design , link.inverse( Yres ) )
		out = np.abs( link( design @ coef ) )
	
	return out if value else coef
##}}}

def std( Y , c_Y = None , m_Y = None , link = IdLinkFct() , value = True ):##{{{
	"""
		SDFC.NonParametric.std
		======================
		
		Estimate standard deviation given a covariate (or not)
		
		Parameters
		----------
		Y     : np.array
			Dataset to fit the mean
		c_Y   : np.array or None
			Covariate(s)
		m_Y   : np.array or float or None
			mean of Y. If None, m = np.mean(Y)
		link  : class based on SDFC.tools.LinkFct
			Link function, default is identity
		value : bool
			If true return value fitted, else return coefficients of fit
		
		Returns
		-------
		The standard deviation
	"""
	out = np.sqrt( var( Y , c_Y , m_Y , link ) )
	if not value:
		if c_Y is None:
			coef = link.inverse(out)
		else:
			if c_Y.ndim == 1: c_Y = c_Y.reshape(-1,1)
			design = np.hstack( ( np.ones((Y.size,1)) , c_Y ) )
			coef,_,_,_ = scl.lstsq( design , link.inverse( out ) )
			out = link( design @ coef )
		return coef
	
	return out
##}}}


####################
## Params classes ##
####################

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
	
	@property
	def value( self ):
		return self.link(self.fit_)
	
	def gradient( self ):
		return self.link.gradient(self.fit_)
	
	def set_coef( self , coef ):
		pass
	
##}}}

class CovariateParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , resample , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.design_ = None
		
		## Build design matrix
		X = kwargs.get( "c_" + self.kind )
		if X.ndim == 1: X = X.reshape(-1,1)
		self.n_features = X.shape[1] + 1
		if resample is None:
			self.design_ = np.hstack( ( np.ones((self.n_samples,1)) , X ) )
		else:
			self.design_ = np.hstack( ( np.ones((self.n_samples,1)) , X[resample,:] ) )
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
	def __init__( self , kind , n_samples , resample , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.coef_      = np.zeros(1)
		self.n_features = 1
		self.fit_       = np.zeros(self.n_samples)
		self.design_    = np.ones( (self.n_samples,1) )
	
	def is_fix(self):
		return False
	
	def update( self ):
		self.fit_ = np.repeat( self.coef_ , self.n_samples ).reshape(-1,1)
	
	def set_coef( self , coef ):
		self.coef_ = np.array([coef]).squeeze()
		self.update()
	
	def design_wo1(self):
		return None
##}}}

class FixParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , resample , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.fit_       = np.array( [self.link.inverse( kwargs.get( "f_" + self.kind ) )] ).reshape(-1,1)
		if self.fit_.size == 1: self.fit_ = np.repeat( self.fit_ , self.n_samples ).reshape(-1,1)
	
	def is_fix(self):
		return True
	
##}}}

class LawParams:##{{{
	def __init__( self , kinds , **kwargs ):
		self.kinds = kinds
		self._dparams = {}
	
	def add_params( self , n_samples , resample , **kwargs ):
		for kind in self.kinds:
			k_param = self.filter( kind , **kwargs )
			config = self.infer_configuration(**k_param)
			if self.is_covariate(config):  self._dparams[kind] = CovariateParam(  kind , n_samples , resample , **k_param )
			if self.is_stationary(config): self._dparams[kind] = StationaryParam( kind , n_samples , resample , **k_param )
			if self.is_fix(config):        self._dparams[kind] = FixParam(        kind , n_samples , resample , **k_param )
	
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
		for c,k in zip(lcoef,self._dparams):
			self._dparams[k].set_coef(c)
	
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
	
	def __init__( self , method , n_bootstrap , alpha ):##{{{
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
		self.method    = method.lower()
		self.params    = {}
		self.coef_     = None
		
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
		header = [ str(type(self)).split(".")[-1][:-2] + " ({})".format(self.method) , "Link" , "Type" , "coef" ]
		if self.confidence_interval is not None:
			header += [ "Quantile {}".format(self.alpha/2) , "Quantile {}".format( 1 - self.alpha / 2 ) ]
		tab.header( header )
		
		## Loop on params
		a = 0
		for k in self.params._dparams:
			p = self.params._dparams[k]
			coef = None if p.coef_ is None else str(p.coef_.squeeze().round(3).tolist())
			row = [ p.kind , str(p.link).split(".")[-1] , str(type(p)).split(".")[-1].split("Param")[0] , coef ]
			if self.confidence_interval is not None:
				if not p.is_fix():
					b = a + p.n_features
					row += [ self.confidence_interval[i,a:b].squeeze().round(3).tolist() for i in range(2) ]
					a = b
				else:
					row += [ str(None) , str(None) ]
			tab.add_row( row )
		return tab.draw() + "\n"
	##}}}
	
	
	def _predict_covariate( self , kind , c_p ): ##{{{
		p = self.params._dparams[kind]
		
		if isinstance(p,CovariateParam) and c_p is not None:
			return p.coef_[0] + c_p.T @ p.coef_[1:]
		return p.value
	##}}}
	
	def _update_coef(func):##{{{
		def wrapper(*args):
			args[0].params.update_coef(args[1])
			return func(*args)
		return wrapper
	##}}}
	
	def _fit_mle( self ):##{{{
		self.optim_result = sco.minimize( self._negloglikelihood , self.params.merge_coef() , jac = self._gradient_nlll , method = "BFGS" )
		self.params.update_coef( self.optim_result.x )
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
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of NormalLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , method , n_bootstrap , alpha )
	##}}}
	
	def __str__(self):##{{{
		return self._to_str()
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	
	@property
	def loc(self):##{{{
		return self.params._dparams["loc"].value
	##}}}
	
	@property
	def scale(self):##{{{
		return self.params._dparams["scale"].value
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
	
	def predict_loc( self , c_loc = None ):##{{{
		"""
		Return location parameter with a new co-variates
		
		Arguments
		---------
		c_loc : np.array or None
			Covariate
		
		Return
		------
		loc : np.array
			Location parameters, if c_loc is None return self.loc
		"""
		return self._predict_covariate( "loc" , c_loc )
	##}}}
	
	def predict_scale( self , c_scale  = None ):##{{{
		"""
		Return scale parameter with a new co-variates
		
		Arguments
		---------
		c_scale : np.array or None
			Covariate
		
		Return
		------
		scale : np.array
			Scale parameters, if c_scale is None return self.scale
		"""
		return self._predict_covariate( "scale" , c_scale )
	##}}}
	
	
	def fit( self , Y , **kwargs ): ##{{{
		"""
		Fit function for NormalLaw
		
		Arguments
		---------
		
		Y       : numpy.ndarray
			Data to fit
		c_loc   : numpy.ndarray or None
			Covariate of loc
		f_loc   : numpy.ndarray or None
			Fix value of loc
		l_loc   : SDFC.tools.LinkFct (optional)
			Link function of loc
		c_scale : numpy.ndarray or None
			Covariate of scale
		f_scale : numpy.ndarray or None
			Fix value of scale
		l_scale : SDFC.tools.LinkFct (optional)
			Link function of scale
		"""
		
		## Bootstrap part
		if self.n_bootstrap > 0:
			self.coefs_bootstrap = []
			for _ in range(self.n_bootstrap):
				idx = np.random.choice( Y.size , Y.size , replace = True )
				self.params = LawParams( kinds = ["loc","scale"] )
				self.params.add_params( n_samples = Y.size , resample = idx , **kwargs )
				self._Y = Y.reshape(-1,1)[idx,:]
				self._fit()
				self.coefs_bootstrap.append( self.coef_ )
			self.coefs_bootstrap = np.array( self.coefs_bootstrap )
			self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
		
		## Fit part
		self.params = LawParams( kinds = ["loc","scale"] )
		self.params.add_params( n_samples = Y.size , resample = None , **kwargs )
		self._Y = Y.reshape(-1,1)
		self._fit()
		del self._Y
		
		
	##}}}
	
	def _fit_moments(self):##{{{
		ploc   = self.params._dparams["loc"]
		pscale = self.params._dparams["scale"]
		
		## Fit loc
		if not ploc.is_fix():
			cm = mean( self._Y , ploc.design_wo1() , value = False , link = ploc.link )
			ploc.set_coef( cm )
		
		## Fit scale
		if not pscale.is_fix():
			pscale.set_coef( std( self._Y , pscale.design_wo1() , m_Y = self.loc , value = False , link = pscale.link ) )
	##}}}
	
	def _fit_mle(self):##{{{
		self._fit_moments()
		AbstractLaw._fit_mle(self)
	##}}}
	
	def _fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self._fit_moments()
		else:
			self._fit_mle()
		self.coef_ = self.params.merge_coef()
	##}}}
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		scale2 = np.power( self.scale , 2 )
		return np.Inf if not np.all( self.scale > 0 ) else np.sum( np.log( scale2 ) ) / 2. + np.sum( np.power( self._Y - self.loc , 2 ) / scale2 ) / 2.
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
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
	loc   = 1. + 0.8 * X_loc# - 0.5 * X_loc**2
	scale = 0.08 * X_scale
	Y     = np.random.normal( loc = loc , scale = scale )
	
	
	## Fit
	norm = NormalLaw( method = "MLE" , n_bootstrap = 100 )
	norm.fit( Y , c_loc = X_loc , c_scale = X_scale , l_scale = sdt.IdLinkFct() )
	print(norm)
	
	print("Done")
