
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


import numpy as np
from .__AbstractLaw import AbstractLaw
from .__Link        import FixedParams

import SDFC.NonParametric as sdnp

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
		self._loc,self._scale = loc.squeeze(),scale.squeeze()
	##}}}
	
	
	## Fit methods
	##============
	
	def _fit_moments( self ): ##{{{
		
		coefs = np.zeros(self._l_global.n_features)
		
		## Find loc
		##=========
		if not isinstance(self._l_global._l_p[0],FixedParams):
			X_loc = self._c_global[0]
			a = sdnp.mean( self._Y , X_loc , self._l_global._l_p[0]._link , False ).squeeze()
			coefs[:self._l_global._s_p[0]] = a
			self.coef_ = coefs
		
		## Find scale
		##===========
		if not isinstance(self._l_global._l_p[1],FixedParams):
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
			self.coef_ = self._l_global.valid_point( self )
	##}}}
	
	def _negloglikelihood( self , coef ): ##{{{
		self.coef_ = coef
		shape = self._Y.shape
		scale2 = np.power( self.scale , 2 )
		if not np.isscalar(scale2): scale2 = scale2.reshape(shape)
		return np.Inf if not np.all( self.scale > 0 ) else np.sum( np.log( scale2 ) ) / 2. + np.sum( np.power( self._Y - self.loc.reshape(shape) , 2 ) / scale2 ) / 2.
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
		T1 = - self._Y * Z / scale**2 + loc * Z / scale**2 + 1 / scale
		jac = self._l_global.jacobian( coef , self._c_global )
		p = 0
		if not isinstance(self._l_global._l_p[0],FixedParams):
			jac[p,:,:] *= T0
			p += 1
		if not isinstance(self._l_global._l_p[1],FixedParams):
			jac[p,:,:] *= T1
		
		return jac.sum( axis = (0,1) )
	##}}}
	
##}}}

