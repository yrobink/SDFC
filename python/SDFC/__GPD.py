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

import numpy         as np
import scipy.special as scs

from SDFC.__AbstractLaw            import AbstractLaw
from SDFC.NonParametric.__mean     import mean
from SDFC.NonParametric.__std      import std
from SDFC.NonParametric.__lmoments import lmoments


#############
## Classes ##
#############

class GPD(AbstractLaw):
	"""
	Class to fit a GPD law with covariates, available methods are:
	
	moments  : use empirical estimator of mean and standard deviation to find loc and scale, possibly with least square
			   regression if covariates are given
	lmoments : Use L-Moments estimation, only in stationary context
	lmoments_experimental: Use non-stationary L-Moments with Quantile Regression, experimental and not published, only
	           used to find an initialization of MLE
	bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of n_mcmc_iteration sample draw from
	           the posterior P(coef_ | Y)
	mle      : Maximum likelihood estimation
	
	WARNING: For this class, f_loc must be always given, because we fit a pareto beyond the loc parameter!
	========
	
	Parameters
	==========
	loc   : location parameter, must be pass with f_loc
	scale : scale parameter
	shape : shape parameter
	"""
	__doc__ += AbstractLaw.__doc__
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of Normal law
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , ["loc","scale","shape"] , method , n_bootstrap , alpha )
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
		
		idx = (self._Y > self.loc)
		Y   = (self._Y[idx] - self.loc[idx]).reshape(-1,1)
		if not pscale.is_fix():
			c_scale = pscale.design_wo1()
			if c_scale is not None:
				c_scale = c_scale.reshape(-1,pscale.n_features-1)[idx]
			self.params.update_coef( std( Y , c_scale , m_Y = 0 , value = False , link = pscale.link ) , "scale" )
		
		if not pshape.is_fix():
			self.params.set_intercept( -1e-8 , "shape" )
	##}}}
	
	def _fit_lmoments( self ): ##{{{
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		idx = (self._Y > self.loc)
		Y   = (self._Y[idx] - self.loc[idx]).reshape(-1,1)
		
		## L-moments
		lmom = lmoments( Y )
		if not pscale.is_fix() and not pshape.is_fix():
			itau     = lmom[0] / lmom[1]
			scale_lm = lmom[0] * ( itau - 1 )
			scale_lm = scale_lm if scale_lm > 0 else 1e-8
			shape_lm = 2 - itau
			self.params.set_intercept( scale_lm , "scale" )
			self.params.set_intercept( shape_lm , "shape" )
		elif not pscale.is_fix():
			scale = lmom[0] * ( 1 - self.shape )
			scale[ np.logical_not(scale > 0) ] = 1e-8
			self.params.update_coef( mean( scale , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
		elif not pshape.is_fix():
			Y /= self.scale[idx].reshape(-1,1)
			lmom = lmoments(Y)
			itau     = lmom[0] / lmom[1]
			self.params.set_intercept( 2 - itau , "shape" )
	##}}}
	
	def _fit_lmoments_experimental( self ): ##{{{
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		idx = (self._Y > self.loc)
		Y   = (self._Y[idx] - self.loc[idx]).reshape(-1,1)
		
		## First step, find lmoments
		c_Y = self.params.merge_covariate()
		if c_Y is None:
			self._fit_lmoments()
			return
		
		c_Y = c_Y[idx.squeeze(),:]
		lmom = lmoments( Y , c_Y )
		
		if not pscale.is_fix() and not pshape.is_fix():
			itau  = lmom[:,0] / lmom[:,1]
			scale = lmom[:,0] * ( itau - 1 )
			shape = 2 - itau
			
			scale_design = pscale.design_wo1()
			if scale_design is not None: scale_design = scale_design[idx.squeeze(),:]
			self.params.update_coef( mean( scale , scale_design , link = pscale.link , value = False ) , "scale" )
			
			shape_design = pshape.design_wo1()
			if shape_design is not None: shape_design = shape_design[idx.squeeze(),:]
			self.params.update_coef( mean( shape , shape_design , link = pshape.link , value = False ) , "shape" )
		elif not pscale.is_fix():
			scale = lmom[:,0].reshape(-1,1) * ( 1 - self.shape )
			scale_design = pscale.design_wo1()
			if scale_design is not None: scale_design = scale_design[idx.squeeze(),:]
			self.params.update_coef( mean( scale , scale_design , link = pscale.link , value = False ) , "scale" )
		elif not pshape.is_fix():
			Y    /= self.scale[idx].reshape(-1,1)
			lmom  = lmoments( Y , pshape.design_wo1() )
			shape = 2 - lmom[:,0] / lmom[:,1]
			shape_design = pshape.design_wo1()
			if shape_design is not None: shape_design = shape_design[idx.squeeze(),:]
			self.params.update_coef( mean( shape , shape_design , link = pshape.link , value = False ) , "shape" )
	##}}}
	
	def _initialization_mle(self):##{{{
		self._fit_lmoments_experimental()
		nlll = self._negloglikelihood(self.coef_)
		grad = self._gradient_nlll(self.coef_)
		
		f_scale = 1
		f_shape = 1
		while ( not nlll < np.inf ) or np.any(np.isnan(grad)):
			pscale = self.params._dparams["scale"]
			pshape = self.params._dparams["shape"]
			
			if pshape.is_fix() and not pscale.is_fix():
				coef_ = np.zeros(pscale.n_features)
				coef_[0] = pscale.link.inverse( 1. * f_scale )
				self.params.update_coef( coef_ , "scale" )
			elif not pshape.is_fix():
				coef_ = np.zeros(pshape.n_features)
				coef_[0] = pshape.link.inverse( 1e-1 / f_shape )
				self.params.update_coef( coef_ , "shape" )
			else:
				self._fit_quantiles()
			f_scale *= 2
			f_shape *= 2
			nlll = self._negloglikelihood(self.coef_)
			grad = self._gradient_nlll(self.coef_)
#		self._fit_moments()
#		nlll_mom = self._negloglikelihood(self.coef_)
#		grad_mom = np.any(np.isnan(self._gradient_nlll(self.coef_)))
#		
#		self._fit_lmoments()
#		nlll_lmo = self._negloglikelihood(self.coef_)
#		grad_lmo = np.any(np.isnan(self._gradient_nlll(self.coef_)))
#		
#		if grad_mom and grad_lmo and ( nlll_mom < np.inf or nlll_lmo < np.inf ):
#			if nlll_mom < nlll_lmo:
#				self._fit_moments()
#			else:
#				self._fit_lmoments()
#		elif grad_mom and nlll_mom < np.inf:
#			self._fit_moments()
#		elif grad_lmo and nlll_lmo < np.inf:
#			self._fit_lmoments()
#		else:
#			pscale = self.params._dparams["scale"]
#			pshape = self.params._dparams["shape"]
#			if not pscale.is_fix():
#				self.params.set_intercept( 1.   , "scale" )
#			if not pshape.is_fix():
#				self.params.set_intercept( 0.01 , "shape" )
		
	##}}}
	
	def _fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self._fit_moments()
		elif self.method == "lmoments":
			self._fit_lmoments()
		elif self.method == "lmoments-experimental":
			self._fit_lmoments_experimental()
	##}}}
	
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		## Impossible scale
		if not np.all( self.scale > 0 ):
			return np.inf
		
		## Remove exponential case
		shape = self.shape
		zero_shape = ( np.abs(shape) < 1e-10 )
		if np.any(zero_shape):
			shape[zero_shape] = -1e-10
		
		
		##
		idx   = (self._Y > self.loc).squeeze()
		loc   = self.loc[idx,:]
		scale = self.scale[idx,:]
		shape = shape[idx,:]
		Z = 1. + shape * ( self._Y[idx,:] - loc ) / scale
		
		if not np.all(Z > 0):
			return np.inf
		res = np.sum( np.log( scale ) + np.log(Z) * ( 1 + 1. / shape ) )
		return res
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		
		## Remove exponential case
		shape = self.shape
		zero_shape = ( np.abs(shape) < 1e-10 )
		if np.any(zero_shape):
			shape[zero_shape] = -1e-10
		
		
		##
		idx   = (self._Y > self.loc).squeeze()
		Y      = self._Y[idx,:]
		loc   = self.loc[idx,:]
		scale = self.scale[idx,:]
		shape = shape[idx,:]
		
		Z        = ( Y - loc ) / scale
		ZZ       = 1. + shape * Z
		exponent = 1. + 1. / shape
		
		grad = np.array([])
		
		pscale = self.params._dparams["scale"]
		if not pscale.is_fix():
			gr_scale   = pscale.gradient()[idx,:]
			A = gr_scale * ( - exponent * shape * Z / ZZ / scale + 1. / scale )
			B = pscale.design_[idx,:].T
			grad_scale =  B @ A
			grad       = np.hstack( (grad,grad_scale.squeeze()) )
		
		pshape = self.params._dparams["shape"]
		if not pshape.is_fix():
			gr_shape   = pshape.gradient()[idx,:].reshape(-1,1)
			grad_shape = pshape.design_[idx,:].T @ ( gr_shape * ( - np.log(ZZ) / shape**2 + exponent * Z / ZZ ) ) if np.all( ZZ > 0 ) else np.repeat(np.nan,pshape.n_features)
			grad       = np.hstack( (grad,grad_shape.squeeze()) )
		return grad
	##}}}

