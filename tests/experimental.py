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
import scipy.optimize as sco
import scipy.linalg as scl
import pandas as pd
import xarray as xr
import SDFC as sd
import SDFC.NonParametric as sdnp

## Plot libraries ##
##==================

import matplotlib as mpl
import matplotlib.pyplot as plt
import CTP.plot as cplt

#mpl.rcParams['font.size'] = 30

##}}}
###############

###############
## Fonctions ##
###############


#############
## Classes ##
#############

## Univ link function
##===================

class UnivariateLink:##{{{
	def __init__(self): pass

	def __call__( self , x ):
		return self.transform(x)
##}}}

class Identity(UnivariateLink): ##{{{
	
	def __init__( self ):
		UnivariateLink.__init__(self)
	
	def transform( self , x ):
		return x
	
	def inverse( self , x ):
		return x
	
	def jacobian( self , x ):
		return np.ones_like(x)
##}}}

class Exponential(UnivariateLink): ##{{{
	
	def __init__( self ):
		UnivariateLink.__init__(self)
	
	def transform( self , x ):
		return np.exp(x)
	
	def inverse( self , x ):
		return np.log(x)
	
	def jacobian( self , x ):
		return np.exp(x)
##}}}


## Global link function
##=====================

class GlobalLink:##{{{
	def __init__( self , *args , **kwargs ):
		self._special_fit_allowed = False
		self._n_features = 0
	
	def transform( self , coef , X ):
		pass
	
	def jacobian( self , coef , X ):
		pass
	
	def pseudo_inverse( self , params , lin_coef , X ):
		pass
	
	@property
	def n_features(self):
		return self._n_features
	
##}}}

class FixedParams(GlobalLink):##{{{
	def __init__( self , value ):
		self.value_ = value
	
	def transform( self , *args , **kwargs ):
		return self.value_
##}}}

class TransformLinear(GlobalLink): ##{{{
	
	def __init__( self , *args , **kwargs ):##{{{
		GlobalLink.__init__( self , *args , **kwargs )
		self._link     = kwargs.get("link")
		if self._link is None:
			self._link = Identity()
	##}}}
	
	def _linear_transform( self , coef , X ): ##{{{
		return coef[0] + X @ coef[1:]
	##}}}
	
	def transform( self , coef , X ):##{{{
		out = self._link.transform( self._linear_transform(coef,X) )
		return out
	##}}}
	
	def jacobian( self , coef , X ): ##{{{
		jac = np.zeros( (coef.size,X.shape[0]) )
		jac[0,:] = 1
		jac[1:,:] = X.T
		return self._link.jacobian( self._linear_transform( coef , X ) ) * jac
	##}}}
	
##}}}

class TensorLink(GlobalLink):##{{{
	
	def __init__( self , l_p , s_p , *args , **kwargs ):##{{{
		GlobalLink.__init__( self , *args , **kwargs )
		self._l_p = [ l if l is not None else TransformLinear() for l in l_p ]
		for i in range(len(self._l_p)):
			if isinstance(self._l_p[i],UnivariateLink):
				self._l_p[i] = TransformLinear( link = self._l_p[i] )
		self._s_p = s_p
		self._n_features = np.sum(self._s_p)
		self._special_fit_allowed = np.all( [isinstance(l,TransformLinear) for l in self._l_p] )
	##}}}
	
	def transform( self , coef , X ): ##{{{
		list_p = []
		ib,ie = 0,0
		for s,l,x in zip(self._s_p,self._l_p,X):
			ie += s
			list_p.append( l.transform( coef[ib:ie] , x ) )
			ib += s
		return list_p
	##}}}
	
	def jacobian( self , coef , X ): ##{{{
		list_jac = []
		ib,ie = 0,0
		jac = np.zeros( (len(self._s_p),coef.size,X[0].shape[0]) )
		for i,(s,l,x) in enumerate(zip(self._s_p,self._l_p,X)):
			ie += s
			jac[i,ib:ie,:] = l.jacobian( coef[ib:ie] , x )
			ib += s
		return jac
	##}}}
	
	def pseudo_inverse( self , params , lin_coef , X ):##{{{
		pass
	##}}}
	
##}}}

## Non parametric functions
##=========================

## Statistical distribution classes
##=================================

class AbstractLaw:##{{{
	
	class _Info(object):##{{{
		def __init__(self):
			self.mle_optim_result   = None
			self.cov_from_optim_mle = False
		
		@property
		def cov_(self):
			if self.cov_from_optim_mle:
				return self.mle_optim_result.hess_inv
			
			return None
	##}}}
	
	## Init functions
	##===============
	
	def __init__( self , name_params : list , method : str ): ##{{{
		self._method      = method.lower()
		self._name_params = name_params
		self._c_global    = None
		self._l_global    = None
		self._coef        = None
		self.info_        = AbstractLaw._Info()
	##}}}
	
	def _init_link( self , **kwargs ): ##{{{
		if kwargs.get("l_global") is not None:
			self._l_global = kwargs.get("l_global")
			self._c_global = kwargs.get("c_global")
		else:
			l_p = []
			c_p = []
			s_p = [] ## size
			for p in self._name_params:
				if kwargs.get("l_{}".format(p)) is not None:
					l_p.append(kwargs.get("l_{}".format(p)))
				elif kwargs.get("f_{}".format(p)) is not None:
					l_p.append( FixedParams( kwargs.get("f_{}".format(p)) ) )
				else:
					l_p.append(None)
				if kwargs.get("c_{}".format(p)) is not None:
					c = kwargs.get("c_{}".format(p)).squeeze()
					if c.ndim == 1: c = c.reshape(-1,1)
					c_p.append(c)
					s_p.append( 1 + c.shape[1] )
				else:
					c_p.append(None)
			self._l_global = TensorLink( l_p , s_p )
			self._c_global = c_p
	##}}}
	
	def _update_coef(func):##{{{
		def wrapper(*args):
			self  = args[0]
			coef_ = args[1]
			self.coef_ = coef_
			self._set_params( *self._l_global.transform( coef_ , self._c_global ) )
			return func(*args)
		return wrapper
	##}}}
	
	
	## Properties
	##===========
	
	@property
	def method(self):##{{{
		return self._method
	##}}}
	
	@property
	def coef_(self):##{{{
		return self._coef
	##}}}
	
	@coef_.setter
	def coef_( self , coef_ ): ##{{{
		self._coef = coef_
		self._set_params( *self._l_global.transform( coef_ , self._c_global ) )
	##}}}
	
	@property
	def cov_(self):##{{{
		return self.info_.cov_
	##}}}
	
	## Fit functions
	##==============
	
	def _random_valid_point(self):##{{{
		"""
		Try to find a valid point in the neighborhood of self.coef_
		"""
		coef_ = self.coef_.copy()
		cov_  = 0.1 * np.identity(coef_.size)
		
		p_coef = coef_.copy()
		n_it   = 1
		while not np.isfinite(self._negloglikelihood(p_coef)) or not np.all(np.isfinite(self._gradient_nlll(p_coef))):
			if n_it % 100 == 0: cov_ *= 2
			p_coef = np.random.multivariate_normal( coef_ , cov_ )
			n_it += 1
		self.coef_ = p_coef
	##}}}
	
	def _fit_MLE(self): ##{{{
		
		self._init_MLE()
		self._random_valid_point()
		
		self.info_.mle_optim_result = sco.minimize( self._negloglikelihood , self.coef_ , jac = self._gradient_nlll , method = "BFGS" )
		self.info_.cov_from_optim_mle = True
		self.coef_ = self.info_.mle_optim_result.x
		
	##}}}
	
	def _fit_Bayesian(self):##{{{
		pass
	##}}}
	
	def fit( self , Y , **kwargs ): ##{{{
		
		## Add Y
		self._Y = Y.reshape(-1,1)
		
		## Init link functions
		self._init_link(**kwargs)
		
		## Now fit
		if self._method not in ["mle","bayesian"] and self._l_global._special_fit_allowed:
			self._special_fit()
		elif self._method == "mle" :
			self._fit_MLE()
		else:
			self._fit_Bayesian()
		
	##}}}
	
##}}}

##===================================

class Normal(AbstractLaw):##{{{
	
	def __init__( self , method = "MLE" ):##{{{
		AbstractLaw.__init__( self , ["loc","scale"] , method )
		self._loc   = None
		self._scale = None
	##}}}
	
	## Properties
	##===========
	
	@property
	def loc(self):##{{{
		return self._loc
	##}}}
	
	@property
	def scale(self):##{{{
		return self._scale
	##}}}
	
	def _set_params( self , loc , scale ):##{{{
		self._loc,self._scale = loc,scale
	##}}}
	
	
	## Fit methods
	##============
	
	def _fit_moments( self ): ##{{{
		
		coefs = np.zeros(self._l_global.n_features)
		
		## Find loc
		##=========
		X_loc = self._c_global[0]
		coefs[:self._l_global._s_p[0]] = sdnp.mean( self._Y , X_loc , self._l_global._l_p[0]._link , False ).squeeze()
		self.coef_ = coefs
		
		## Find scale
		##===========
		X_scale = self._c_global[1]
		coefs[self._l_global._s_p[0]:] = sdnp.std( self._Y , X_scale , self.loc , self._l_global._l_p[1]._link , False ).squeeze()
		self.coef_ = coefs
	##}}}
	
	def _special_fit( self ):##{{{
		if self.method == "moments":
			self._fit_moments()
	##}}}
	
	def _init_MLE( self ): ##{{{
		if self._l_global._special_fit_allowed:
			self._fit_moments()
		else:
			pass
	##}}}
	
	def _negloglikelihood( self , coef ): ##{{{
		self.coef_ = coef
		scale2 = np.power( self.scale , 2 )
		shape = self._Y.shape
		return np.Inf if not np.all( self.scale > 0 ) else np.sum( np.log( scale2 ) ) / 2. + np.sum( np.power( self._Y - self.loc.reshape(shape) , 2 ) / scale2.reshape(shape) ) / 2.
	##}}}
	
	def _gradient_nlll( self , coef ): ##{{{
		self.coef_ = coef
		## Parameters
		shape = self._Y.shape
		loc   = self.loc.reshape(shape)
		scale = self.scale.reshape(shape)
		Z     = ( self._Y - loc ) / scale
		
		## Compute gradient
		T0 = - Z / scale
		T1 = - Y * Z / scale**2 + loc * Z / scale**2 + 1 / scale
		jac = self._l_global.jacobian( coef , self._c_global )
		jac[0,:,:] *= T0.T
		jac[1,:,:] *= T1.T
		
		return jac.sum( axis = (0,2) )
	##}}}
	
##}}}


## GEV part
##=========

class GEVPrLink(GlobalLink):##{{{
	def __init__( self , *args , **kwargs ):
		GlobalLink.__init__( self , *args , **kwargs )
	
	def transform( self , coef , X ):
		E = np.exp( coef[3] / coef[0] * X[:,0] )
		loc   = coef[0] * E
		scale = coef[1] * E
		shape = coef[2] + np.zeros_like(X[:,0])
		return loc,scale,shape
	
	def jacobian( self , coef , X ):
		E = np.exp( coef[3] / coef[0] * X[:,0] )
		jac = np.zeros( (3 , 4 , X[:,0].size) )
		jac[0,0,:] = E - coef[3] * X[:,0] / coef[0] * E
		jac[1,0,:] = - coef[1] * coef[3] * X[:,0] / coef[0]**2 * E
		jac[1,1,:] = E
		jac[2,2,:] = 1
		jac[0,3,:] = X[:,0] * E
		jac[1,3,:] = coef[1] * X[:,0] * E / coef[0]
		
		return jac
	
	def pseudo_inverse( self , params , lin_coef , X ):
		
		coef = np.zeros(4)
		design = np.stack( (np.ones_like(X),X) , -1 ).squeeze()
		
		idxloc   = np.isfinite(np.log(params[:,0]))
		idxscale = np.isfinite(np.log(params[:,1]))
		resloc   = scl.lstsq( design[idxloc,:]   , np.log(params[idxloc,0]) )
		resscale = scl.lstsq( design[idxscale,:] , np.log(params[idxscale,1]) )
		coef[0] = np.exp(resloc[0][0])
		coef[1] = np.exp(resscale[0][0])
		coef[2] = params[:,2].mean()
		
		alphaloc   = resloc[0][1] * coef[0]
		alphascale = resscale[0][1] * coef[0]
		coef[3]    = ( alphaloc + alphascale ) / 2
		
		return coef
##}}}

class GEV:##{{{
	
	def __init__( self , link , restart_fit = 0 ):##{{{
		self.link  = link
		self.optim = None
		self.coef_ = None
		self.restart_fit = restart_fit
		self.n_restart   = 0
	##}}}
	
	def negloglikelihood( self , coef ): ##{{{
		self.coef_ = coef
		loc,scale,shape = self.link.transform( self.coef_ , self._X )
		return - np.sum( sc.genextreme.logpdf( self._Y , loc = loc , scale = scale , c = -shape ) )
	##}}}

	def grad_nlll( self , coef ): ##{{{
		self.coef_ = coef
		jac = self.link.jacobian( self.coef_ , self._X )
		loc,scale,shape = self.link.transform( self.coef_ , self._X )
		
		
		Z  = ( self._Y - loc ) / scale
		ZZ = 1 + shape * Z
		ZZi = ZZ**( - 1 / shape )
		kappa = ( 1 + 1 / shape ) / ZZ - ZZi / (shape * ZZ)
		
		T0 = - shape * kappa / scale
		T1 = 1 / scale - shape * Z / scale * kappa
		T2 = np.log(ZZ) * ( ZZi - 1 ) / shape**2 + Z * kappa
		
		jac[0,:] *= T0
		jac[1,:] *= T1
		jac[2,:] *= T2
		
		return jac.sum( axis = (0,2) )
	##}}}
	
	def init_mle( self ):##{{{
		sdgev = sd.GEV( method = "lmoments-experimental" )
		sdgev.fit( self._Y , c_loc = np.exp(self._X) , c_scale = np.exp(self._X) )
		params = np.stack( (sdgev.loc,sdgev.scale,sdgev.shape) ).squeeze().T
		minit   = self.link.pseudo_inverse( params , sdgev.coef_ , self._X )
		mcov    = np.identity(4)
		
		## Perturbation if the init point is not well defined
		init = minit.copy()
		nit = 0
		while not np.isfinite(self.negloglikelihood(init)) or not np.isfinite(self.grad_nlll(init).all()):
			init = np.random.multivariate_normal( mean = minit , cov = mcov )
			nit += 1
			if nit % 100 == 0: mcov *= 2
		
		self.coef_ = init
	##}}}
	
	def fit( self , Y , X ):##{{{
		self._Y = Y.squeeze()
		self._X = X.reshape(-1,1) if X.ndim == 1 else X
		self.init_mle()
		self.optim = sco.minimize( self.negloglikelihood , self.coef_ , jac = self.grad_nlll )
		
		restart_init = self.restart_fit
		if not self.optim.success and self.restart_fit > 0:
			self.restart_fit -= 1
			self.fit( Y , X )
		self.n_restart   = restart_init - self.restart_fit
		self.restart_fit = restart_init
		
	##}}}
##}}}

def test_GEV(): ##{{{
	n_samples = 2000
	t,X,_,_ = sd.tools.Dataset.covariates( n_samples )
	X = X.reshape(-1,1)
	l_global = GEVPrLink()
	
	
	coef_true = np.array([0.5,1.5,-0.2,1.3])
	loc,scale,shape = l_global.transform( coef_true , X )
	Y = sc.genextreme.rvs( loc = loc , scale = scale , c = -shape )
	
	gev = GEV( l_global , 10 )
	gev.fit( Y , X )
	
	print(coef_true)
	print(gev.optim.x)
	print(gev.n_restart)
##}}}


##########
## main ##
##########

if __name__ == "__main__":
	np.seterr( all = "ignore" )
	
	## Build data
	##===========
	n_samples = 2000
	t,X_loc,X_scale,_ = sd.tools.Dataset.covariates( n_samples )
	X_loc   = X_loc.reshape(-1,1)
	X_scale = X_scale.reshape(-1,1)
	
	loc   = 1. + 3 * X_loc
	scale = np.exp( 2 * X_scale )
#	scale = np.zeros(n_samples) + 0.5
	
	Y = np.random.normal( loc = loc , scale = scale )
	
	## Fit
	##====
	l_global = TensorLink( [Identity() , Exponential()] , [2,2] )
	c_global = [X_loc,X_scale]
	norm = Normal( method = "mle" )
#	norm.fit( Y , c_global = c_global , l_global = l_global )
#	norm.fit( Y , c_loc = X_loc , c_scale = X_scale , l_scale = Exponential() )
	norm.fit( Y , f_loc = loc , c_scale = X_scale , l_scale = Exponential() )
	print(norm.info_.mle_optim_result)
	
	## And plot it
	##============
	if True:
		nrow,ncol = 1,1
		fig = plt.figure( figsize = cplt.figsize(nrow,ncol) )
		
		ax = fig.add_subplot( nrow , ncol , 1 )
		ax.plot( t , Y , color = "blue" , linestyle = "" , marker = "." )
		ax.plot( t , norm.loc , color = "red" )
		ax.plot( t , norm.loc + norm.scale , color = "red" , linestyle = "--" )
		ax.plot( t , norm.loc - norm.scale , color = "red" , linestyle = "--" )
#		ax.plot( t , X_loc )
#		ax.plot( t , np.exp(-X_scale) )
		
		fig.set_tight_layout(True)
		plt.show()
	
	print("Done")
