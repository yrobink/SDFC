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

import numpy as np

from SDFC.__AbstractLaw        import AbstractLaw
from SDFC.NonParametric.__mean import mean


#############
## Classes ##
#############

class Exponential(AbstractLaw):
	"""
	Class to fit an Exponential law with covariates, available methods are:
	
	moments  : use empirical estimator of mean and standard deviation to find loc and scale, possibly with least square
			   regression if covariates are given
	bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of n_mcmc_iteration sample draw from
	           the posterior P(coef_ | Y)
	mle      : Maximum likelihood estimation
	
	Parameters
	==========
	scale : scale parameter
	"""
	__doc__ += AbstractLaw.__doc__
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of Exponential law
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments", "lmoments", "quantiles" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , ["scale"] , method , n_bootstrap , alpha )
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
	
	
	def _fit_moments(self):##{{{
		
		pscale = self.params._dparams["scale"]
		
		## Fit loc
		if not pscale.is_fix():
			self.params.update_coef( mean( self._Y , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
	##}}}
	
	def _initialization_mle(self):##{{{
		self._fit_moments()
	##}}}
	
	def _fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self._fit_moments()
	##}}}
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		if not np.all(self.scale > 0):
			return np.Inf
		
		return np.sum( np.log(self.scale) - self._Y / self.scale )
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		grad_scale = np.zeros_like(coef) + np.nan
		
		pscale = self.params._dparams["scale"]
		if np.all(self.scale > 0):
			grad_scale = pscale.design_.T @ ( ( 1. / self.scale - self._Y / self.scale**2 ) * pscale.gradient() )
		
		return grad_scale.squeeze()
	##}}}


