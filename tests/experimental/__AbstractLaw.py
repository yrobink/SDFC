
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
import scipy.optimize as sco

from .core.__LHS import LHS
from .core.__RHS import RHS


###############
## Class(es) ##
###############

class AbstractLaw:
	
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
	
	def __init__( self , names : list , method : str ): ##{{{
		self._method = method.lower()
		self._lhs    = LHS(names,0)
		self._rhs    = RHS(self._lhs)
		self.info_   = AbstractLaw._Info()
	##}}}
	
	
	## Properties
	##===========
	
	@property
	def method(self):##{{{
		return self._method
	##}}}
	
	@property
	def coef_(self):##{{{
		return self._rhs.coef_
	##}}}
	
	@coef_.setter
	def coef_( self , coef_ ): ##{{{
		self._rhs.coef_ = coef_
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
		
		## Init LHS/RHS
		self._lhs.n_samples = Y.size
		self._rhs.build(**kwargs)
		
		## Now fit
		if self._method not in ["mle","bayesian"] and self._rhs.l_global._special_fit_allowed:
			self._special_fit()
		elif self._method == "mle" :
			self._fit_MLE()
		else:
			self._fit_Bayesian()
		
	##}}}


