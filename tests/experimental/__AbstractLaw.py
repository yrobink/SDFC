
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
from .__Link        import TensorLink
from .__Link        import FixedParams


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
			s_p = []
			n_samples = self._Y.size
			for p in self._name_params:
				if kwargs.get("c_{}".format(p)) is not None:
					c = kwargs.get("c_{}".format(p)).squeeze()
					if c.ndim == 1: c = c.reshape(-1,1)
					c_p.append(c)
					s_p.append( 1 + c.shape[1] )
				else:
					c_p.append(None)
					if kwargs.get("f_{}".format(p)) is not None:
						s_p.append(0)
					else:
						s_p.append(1)
				
				if kwargs.get("l_{}".format(p)) is not None:
					l_p.append(kwargs.get("l_{}".format(p)))
				elif kwargs.get("f_{}".format(p)) is not None:
					l_p.append( FixedParams( kwargs.get("f_{}".format(p)) , n_samples = n_samples , n_features = 0 ) )
				else:
					l_p.append(None)
				
			self._l_global = TensorLink( l_p , s_p , n_features = np.sum(s_p) , n_samples = n_samples )
			self._c_global = c_p
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


