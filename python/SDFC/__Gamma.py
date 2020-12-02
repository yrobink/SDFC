# -*- coding: utf-8 -*-

## Copyright(c) 2020 Yoann Robin
## 
## This file is part of SDFC.
## 
## SDFC is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## SDFC is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with SDFC.  If not, see <https://www.gnu.org/licenses/>.


##############
## Packages ##
##############

import numpy          as np
import scipy.stats    as sc
import scipy.linalg   as scl
import scipy.optimize as sco
import scipy.special  as scp

from SDFC.__AbstractLaw        import AbstractLaw
from SDFC.NonParametric.__mean import mean
from SDFC.NonParametric.__var  import var


#############
## Classes ##
#############

class Gamma(AbstractLaw):
	"""
	Class to fit a Gamma law with covariates, available methods are:
	
	moments  : use empirical estimator of mean and standard deviation to find loc and scale, possibly with least square
			   regression if covariates are given
	bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of n_mcmc_iteration sample draw from
	           the posterior P(coef_ | Y)
	mle      : Maximum likelihood estimation
	
	Parameters
	==========
	scale : scale parameter
	shape : shape parameter
	"""
	__doc__ += AbstractLaw.__doc__
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of Gamma law
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		"""
		AbstractLaw.__init__( self , ["scale","shape"] , method , n_bootstrap , alpha )
	##}}}
	
	def __str__(self):##{{{
		return self._to_str()
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	@property
	def scale(self):##{{{
		return self.params._dparams["scale"].value
	##}}}
	
	@property
	def shape(self):##{{{
		return self.params._dparams["shape"].value
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
	
	def predict_shape( self , c_shape  = None ):##{{{
		"""
		Return scale parameter with a new co-variates
		
		Arguments
		---------
		c_shape : np.array or None
			Covariate
		
		Return
		------
		shape : np.array
			Shape parameters, if c_scale is None return self.shape
		"""
		return self._predict_covariate( "shape" , c_shape )
	##}}}
	
	
	def _fit_moments(self):##{{{
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		n_samples = pscale.n_samples
		
		if not pscale.is_fix() and not pshape.is_fix():
			mX = np.ones( (n_samples,1) )
			vX = np.ones( (n_samples,1) )
			for i in range(1,pscale.n_features):
				for j in range(pshape.n_features):
					mX = np.hstack( (mX,np.reshape( pscale.design_[:,i]    * pshape.design_[:,j] , (n_samples,1) ) ) )
					vX = np.hstack( (vX,np.reshape( pscale.design_[:,i]**2 * pshape.design_[:,j] , (n_samples,1) ) ) )
			m = mean( self._Y , mX[:,1:] )
			v = var(  self._Y , vX[:,1:] )
			
			idx  = np.logical_or( np.abs(m) < 1e-8 , v < 1e-8 )
			cidx = np.logical_not(idx)
			scale = np.zeros_like(m)
			shape = np.zeros_like(m)
			scale[cidx] = v[cidx] / m[cidx]
			shape[cidx] = m[cidx]**2 / v[cidx]
			
			if np.any(idx):
				scale[idx] = scale[cidx].min()
				shape[idx] = shape[cidx].min()
			
			self.params.update_coef( mean( scale , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
			self.params.update_coef( mean( shape , pshape.design_wo1() , value = False , link = pshape.link ) , "shape" )
		elif pscale.is_fix():
			
			m = mean( self._Y  , pshape.design_wo1() * self.scale )
			v = var(  self._Y  , pshape.design_wo1() * self.scale**2 )
			
			shape = m**2 / v
			self.params.update_coef( mean( shape , pshape.design_wo1() , value = False , link = pshape.link ) , "shape" )
			
		elif pshape.is_fix():
			m = mean( self._Y  , pscale.design_wo1()    * self.shape )
			v = var(  self._Y  , pscale.design_wo1()**2 * self.shape )
			
			scale = v / m
			self.params.update_coef( mean( scale , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
	##}}}
	
	def _initialization_mle(self):##{{{
		self._fit_moments()
	##}}}
	
	def _fit( self ): ##{{{
		if self.method == "moments":
			self._fit_moments()
	##}}}
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		if not np.all(self.scale > 0) or not np.all(self.shape > 0) or not np.all(self._Y > 0):
			return np.Inf
		
		return np.sum( self._Y / self.scale + scp.loggamma(self.shape) + self.shape * np.log(self.scale) - (self.shape-1) * np.log(self._Y) )
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		grad = np.array( [] )
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		if np.all(self.scale > 0) and np.all(self.shape > 0) and np.all(self._Y > 0):
			if not pscale.is_fix():
				grad_scale = pscale.design_.T @ ( ( self.shape / self.scale - self._Y / self.scale**2 ) * pscale.gradient() )
				grad = np.hstack( (grad,grad_scale.squeeze()) )
			if not pshape.is_fix():
				grad_shape = pshape.design_.T @ ( ( scp.digamma(self.shape) + np.log(self.scale) - np.log(self._Y) ) * pshape.gradient() )
				grad = np.hstack( (grad,grad_shape.squeeze()) )
		else:
			grad = np.zeros( coef.size ) + np.nan
		return grad
	##}}}





