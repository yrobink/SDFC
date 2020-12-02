
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


###############
## Class(es) ##
###############


class UnivariateLink:##{{{
	def __init__(self): pass

	def __call__( self , x ):
		return self.transform(x)
##}}}

class ULIdentity(UnivariateLink): ##{{{
	
	def __init__( self ):
		UnivariateLink.__init__(self)
	
	def transform( self , x ):
		return x
	
	def inverse( self , x ):
		return x
	
	def jacobian( self , x ):
		return np.ones_like(x)
##}}}

class ULExponential(UnivariateLink): ##{{{
	
	def __init__( self ):
		UnivariateLink.__init__(self)
	
	def transform( self , x ):
		return np.exp(x)
	
	def inverse( self , x ):
		return np.log(x)
	
	def jacobian( self , x ):
		return np.exp(x)
##}}}


class AbstractLink:##{{{
	"""
	SDFC.tools.AbstractLink
	=======================
	
	Abstract base class for link function
	
	"""
	def __init__( self ):
		pass
	
	def __str__(self):
		return "SDFC.tools.LinkFct"
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		"""
		Evaluation of link function
		"""
		pass
	
	def gradient( self , x ):
		"""
		Gradient of link function
		"""
		pass
	
	def inverse( self , x ):
		"""
		Inverse of link function
		"""
		pass
##}}}

class ChainLink(AbstractLink): ##{{{
	"""
	SDFC.tools.ChainLink
	====================
	
	This class is used to chain two link functions, i.e. $link1 \circ link0$
	
	"""
	def __init__( self , link1 , link0 ):
		self.link1 = link1
		self.link0 = link0
	
	def __str__(self):
		return "SDFC.tools.ChainLink between {} and {}".format(str(self.link1,self.link0))
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return self.link1( self.link0( x ) )
	
	def gradient( self , x ):
		return self.link0.gradient(x) * self.link1.gradient( self.link0(x) )
	
	def inverse( self , x ):
		return self.link0.inverse( self.link1.inverse( x ) )
##}}}

class IdLink(AbstractLink): ##{{{
	"""
	SDFC.tools.IdLink
	=================
	
	Identity link function
	
	"""
	
	def __init__(self):
		pass
	
	def __str__(self):
		return "SDFC.tools.IdLink"
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return x
	
	def gradient( self , x ):
		return np.ones( x.shape )
	
	def inverse( self , x ):
		return x
##}}}

class InverseLink(AbstractLink): ##{{{
	"""
	SDFC.tools.InverseLink
	======================
	
	Inverse link function, i.e.:
		f(x) = 1/x
		f^{-1}(x) = 1/x
		df(x) = - 1 / x**2
	
	"""
	def __init__(self):
		pass
	
	def __str__(self):
		return "SDFC.tools.InverseLink"
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return 1. / x
	
	def gradient( self , x ):
		return - 1. / x**2
	
	def inverse( self , x ):
		return 1. / x
##}}}

class ExpLink(AbstractLink):##{{{
	"""
	SDFC.tools.ExpLink
	==================
	
	Exponential link function, i.e.:
		f(x) = exp(s*x) + b
		f^{-1}(x) = log(x-b) / s
		df(x) = s*exp(s*x)
	This function is used to bound a variable into level b, by upper if s > 0 or lower if s < 0.
	"""
	def __init__( self , b = 0 , s = 1 ):
		self.b = b
		self.s = s
	
	def __str__(self):
		return "SDFC.tools.ExpLink ({}{})".format( ">" if self.s > 0 else "<" , self.b )
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return np.exp(self.s * x) + self.b
	
	def gradient( self , x ):
		return self.s * np.exp(self.s * x)
	
	def inverse( self , x ):
		return np.log(x - self.b) / self.s
##}}}

class LogitLink(AbstractLink):##{{{
	"""
	SDFC.tools.LogitLink
	====================
	
	Logit link function, i.e.:
		f(x) = ( b - a ) / ( 1 + exp(-sx) ) + b
		f^{-1}(x) = - log( (b-a) / (x-b) - 1) / s
		df(x) = s * (b-a) * exp( -sx ) / ( 1 + exp(-sx) )**2
	
	This function constrain the parameters estimated between a and b, and the parameters s control the growth between a and b.
	
	
	"""
	def __init__( self , a = 0 , b = 1 , s = 1 ):
		self.a      = a
		self.b      = b
		self.s      = s
	
	def __str__(self):
		return "SDFC.tools.LogitLink, ({}< x{} < {})".format(self.a,self.s,self.b)
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return (self.b - self.a) / ( 1. + np.exp( - x * self.s ) ) + self.a
	
	def gradient( self , x ):
		e = np.exp( - self.s * x )
		return self.s * (self.b - self.a) * e / ( 1 + e )**2
	
	def inverse( self , x ):
		x = np.array( [x] ).ravel()
		idx_lo = x < self.a
		idx_up = x > self.b
		x[idx_lo] = self.a + 1e-3
		x[idx_up] = self.b - 1e-3
		return - np.log( (self.b - self.a) / ( x - self.a ) - 1 ) / self.s
##}}}

class SemiBoundedLink(AbstractLink):##{{{
	"""
	SDFC.tools.SemiBoundedLink
	==========================
	A simple semi bounded function, use it if ExpLink has too many overflow (Here the gradient is not well defined at x = b, and the inverse also).
	Values:
		- x < b : f(x) = sx
		- x > b : f(x) = b
	"""
	def __init__( self , b = 0 , s = - 1 ):
		"""
		Parameters
		----------
		b : float, optional
			Bound. The default is 0.
		s : float, optional
			Slope. The default is -1.
		
		Returns
		-------
		None.
		"""
		self.b = b
		self.s = s
	
	def __str__(self):
		return "SDFC.tools.SemiBoundedLink"
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		if np.isscalar(x):
			return self.s * x if x < self.b else self.b
		return np.where( x < self.b , self.s * x , self.b )
	
	def gradient( self , x ):
		return np.where( x < self.b , self.s , 0. )
	
	def inverse( self , x ):
		return self.__call__(x)
##}}}

class BoundedLink(AbstractLink):##{{{
	"""
	SDFC.tools.BoundedLink
	======================
	A simple bounded function, use it if LogitLink has too many overflow (Here the gradient is not well defined at x = a and x = b, and the inverse also).
	Values:
		- x < a : f(x) = a
		- x > b : f(x) = b
		- a <= x <= b : f(x) = x
	"""
	def __init__( self , a , b ):
		self.a = a
		self.b = b
	
	def __str__(self):
		return "SDFC.tools.BoundedLink ({}<x<{})".format(self.a,self.b)
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		if np.isscalar(x):
			if x < self.a:
				return self.a
			elif x > self.b:
				return self.b
			else:
				return x
		return np.where( x > self.a , np.where( x < self.b , x , self.b ) , self.a )
	
	def gradient( self , x ):
		return np.where( (self.a < x) &  (x < self.b) , 1. , 0. )
	
	def inverse( self , x ):
		return self.__call__(x)
##}}}

