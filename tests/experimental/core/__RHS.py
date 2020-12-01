
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

from experimental.__Link import UnivariateLink
from experimental.__Link import Identity

###############
## Class(es) ##
###############


class GlobalLink2:##{{{
	
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

class FixedParams2(GlobalLink2):##{{{
	
	def __init__( self , value , *args , **kwargs ):##{{{
		kwargs["n_features"] = 0
		GlobalLink2.__init__( self , *args , **kwargs )
		self.value_ = np.array([value])
		if self.value_.size == 1:
			self.value_ = np.repeat( value , self.n_samples )
	##}}}
	
	def transform( self , *args , **kwargs ):##{{{
		return self.value_
	##}}}
	
##}}}

class TransformLinear2(GlobalLink2): ##{{{
	
	def __init__( self , *args , **kwargs ):##{{{
		GlobalLink2.__init__( self , *args , **kwargs )
		self._l     = kwargs.get("l")
		self._c     = kwargs.get("c")
		if self._l is None: self._l = Identity()
		self.design_ = np.ones( (self.n_samples,1) )
		if self._c is not None:
			self.design_ = np.hstack( (self.design_,self._c) )
		
	##}}}
	
	def _linear_transform( self , coef , X ): ##{{{
		return self.design_ @ coef
	##}}}
	
	def transform( self , coef , X ):##{{{
		out = self._l.transform( self._linear_transform(coef,X) )
		return out
	##}}}
	
	def jacobian( self , coef , X ): ##{{{
#		jac = np.zeros( (self.n_samples,self.n_features) )
#		jac[:,0]  = 1
#		jac[:,1:] = X
		return self._l.jacobian( self._linear_transform( coef , X ).reshape(-1,1) ) * self.design_
	##}}}
	
##}}}

class TensorLink2(GlobalLink2):##{{{
	
	def __init__( self , l_p , s_p , *args , **kwargs ):##{{{
		GlobalLink2.__init__( self , *args , **kwargs )
		self._l_p = l_p
		self._s_p = s_p
		self._special_fit_allowed = np.all( [isinstance(l,(TransformLinear2,FixedParams2)) for l in self._l_p] )
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



class LHS:
	def __init__( self , names : list , n_samples : int ):
		self.names     = names
		self.n_lhs     = len(self.names)
		self.n_samples = n_samples
		self._values   = { n : None for n in self.names }
		self.jacobian_ = None
		self._fixed    = { n : False for n in self.names }
	
	def is_fixed( self , name ):
		return self._fixed.get(name)
	
	@property
	def values_( self ):
		return self._values
	
	@values_.setter
	def values_( self , values ):
		for n,v in zip(self.names,values):
			self._values[n] = v


class RHS:
	"""
	Class that contains coef_ to fit, covariates and link function. This class
	contains:
	- self.lhs_ : the Left Hand Side part, updated by this class
	- self.c_global : list of covariate
	- self.l_global : global link function
	- self.s_global : length of coef per params of lhs
	
	The parameter self.coef_ is a property, and when it is set the LHS is
	updated accordingly
	
	"""
	def __init__( self , lhs_ : LHS ): ##{{{
		"""
		d
		"""
		self.lhs_     = lhs_
		self.c_global = None
		self.l_global = None
		self.s_global = None
		self._coef_   = None
	##}}}
	
	def build( self , **kwargs ):##{{{
		"""
		Here five kinds of arguments can be passed:
		- c_<param> : covariate of the param,
		- l_<param> : link function of the param,
		- f_<param> : fixed values of the LHS
		- c_global  : list of all covariates, sorted by lhs order
		- l_global  : global link function generated the LHS
		If c_global is set, all arguments (except l_global) are ignored
		"""
		
		## If global covariate and link functions are defined, just set it
		##================================================================
		if kwargs.get("c_global") is not None:
			self.c_global = kwargs["c_global"]
			self.l_global = kwargs.get("l_global")
			if self.l_global is None:
				self.l_global = Id
			return
		
		## Else loop on lhs to find global parameters
		##===========================================
		self.c_global = []
		self.s_global = []
		l_global      = [] ## This list will be passed to TensorLink2
		for lhs in self.lhs_.names:
			## Start with covariate
			if kwargs.get("c_{}".format(lhs)) is not None:
				c = kwargs["c_{}".format(lhs)].squeeze()
				if c.ndim == 1: c = c.reshape(-1,1)
				self.c_global.append(c)
				self.s_global.append(1 + c.shape[1])
			## No covariate, two choices : lhs is 1d or fixed
			else:
				self.c_global.append(None)
				if kwargs.get("f_{}".format(lhs)) is not None:
					self.s_global.append(0)
					self.lhs_._fixed["lhs"] = True
				else:
					self.s_global.append(1)
			
			## Now the link functions
			if kwargs.get("f_{}".format(lhs)) is not None:
				l_global.append( FixedParams2( kwargs["f_{}".format(lhs)] , n_samples = self.lhs_.n_samples ) )
			else:
				l = kwargs.get("l_{}".format(lhs))
				if l is None or issubclass(l.__class__,UnivariateLink):
					l = TransformLinear2( c = self.c_global[-1] , l = l , n_samples = self.lhs_.n_samples )
				l_global.append(l)
		
		self.l_global = TensorLink2( l_global , self.s_global , n_features = np.sum(self.s_global) , n_samples = self.lhs_.n_samples )
	##}}}
	
	## Properties
	## {{{
	
	@property
	def coef_(self):
		return self._coef_
	
	@coef_.setter
	def coef_( self , coef_ ):
		self._coef_         = np.array( [coef_] ).squeeze()
		self.lhs_.values_   = self.l_global.transform( self.coef_ , self.c_global )
		self.lhs_.jacobian_ = self.l_global.jacobian( self.coef_ , self.c_global )
	
	##}}}

