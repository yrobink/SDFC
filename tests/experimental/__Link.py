
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
	
	def __init__( self , *args , **kwargs ):##{{{
		self._special_fit_allowed = False
		self._n_features = kwargs.get("n_features")
		self._n_samples  = kwargs.get("n_samples")
	##}}}
	
	def transform( self , coef , X ):##{{{
		pass
	##}}}
	
	def jacobian( self , coef , X ):##{{{
		pass
	##}}}
	
	@property
	def n_features(self):##{{{
		return self._n_features
	##}}}
	
	@property##{{{
	def n_samples(self):
		return self._n_samples
	##}}}
	
##}}}

class FixedParams(GlobalLink):##{{{
	
	def __init__( self , value , *args , **kwargs ):##{{{
		GlobalLink.__init__( self , *args , **kwargs )
		self.value_ = np.array([value])
		if self.value_.size == 1:
			self.value_ = np.repeat( value , self.n_samples )
	##}}}
	
	def transform( self , *args , **kwargs ):##{{{
		return self.value_
	##}}}
##}}}

class TransformLinear(GlobalLink): ##{{{
	
	def __init__( self , *args , **kwargs ):##{{{
		GlobalLink.__init__( self , *args , **kwargs )
		self._link     = kwargs.get("link")
		if self._link is None: self._link = Identity()
	##}}}
	
	def _linear_transform( self , coef , X ): ##{{{
		if X is None: return np.repeat( coef[0] , self.n_samples ).reshape(self.n_samples,1)
		return (coef[0] + X @ coef[1:]).reshape(self.n_samples,1)
	##}}}
	
	def transform( self , coef , X ):##{{{
		out = self._link.transform( self._linear_transform(coef,X) )
		return out
	##}}}
	
	def jacobian( self , coef , X ): ##{{{
		jac = np.zeros( (self.n_samples,self.n_features) )
		jac[:,0]  = 1
		jac[:,1:] = X
		return self._link.jacobian( self._linear_transform( coef , X ) ) * jac
	##}}}
	
##}}}

class TensorLink(GlobalLink):##{{{
	
	def __init__( self , l_p , s_p , *args , **kwargs ):##{{{
		GlobalLink.__init__( self , *args , **kwargs )
		self._l_p = [ l if l is not None else TransformLinear(n_features=s,n_samples=self.n_samples) for s,l in zip(s_p,l_p) ]
		for i in range(len(self._l_p)):
			if isinstance(self._l_p[i],UnivariateLink):
				self._l_p[i] = TransformLinear( *args , link = self._l_p[i] , n_features = s_p[i] , n_samples = self.n_samples )
		self._s_p = s_p
		self._special_fit_allowed = np.all( [isinstance(l,(TransformLinear,FixedParams)) for l in self._l_p] )
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
		c_size = 1
		for x in X:
			if x is not None:
				c_size = x.shape[0]
		jac = np.zeros( (np.nonzero(self._s_p)[0].size,self.n_samples,self.n_features) )
		i = 0
		for s,l,x in zip(self._s_p,self._l_p,X):
			if s > 0:
				ie += s
				jac[i,:,ib:ie] = l.jacobian( coef[ib:ie] , x )
				ib += s
				i += 1
		return jac
	##}}}
	
##}}}

