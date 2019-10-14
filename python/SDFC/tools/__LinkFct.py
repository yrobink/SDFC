# -*- coding: utf-8 -*-

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SDFC (Statistical    ##
## Distribution Fit with Covariates) library. This library makes it possible    ##
## to regress the parameters of some statistical law with co-variates.          ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SDFC (Statistical Distribution Fit with Covariates). Cette librairie         ##
## permet de calculer de regresser les parametres de lois statistiques selon    ##
## plusieurs co-variables                                                       ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################

###############
## Libraries ##
###############

import numpy as np


###########
## Class ##
###########

class LinkFct:##{{{
	"""
	SDFC.tools.LinkFct
	==================
	
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

class ChainLinkFct(LinkFct): ##{{{
	"""
	SDFC.tools.ChainLinkFct
	=======================
	
	This class is used to chain two link functions, i.e. $linkFct1 \circ linkFct0$
	
	"""
	def __init__( self , linkFct1 , linkFct0 ):
		self.linkFct1 = linkFct1
		self.linkFct0 = linkFct0
	
	def __str__(self):
		return "SDFC.tools.ChainLinkFct between {} and {}".format(str(self.linkFct1,self.linkFct2))
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return self.linkFct1( self.linkFct0( x ) )
	
	def gradient( self , x ):
		return self.linkFct0.gradient(x) * self.linkFct1.gradient( self.linkFct0(x) )
	
	def inverse( self , x ):
		return self.linkFct0.inverse( self.linkFct1.inverse( x ) )
##}}}

class IdLinkFct(LinkFct): ##{{{
	"""
	SDFC.tools.IdLinkFct
	====================
	
	Identity link function
	
	"""
	
	def __init__(self):
		pass
	
	def __str__(self):
		return "SDFC.tools.IdLinkFct"
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return x
	
	def gradient( self , x ):
		return np.ones( x.shape )
	
	def inverse( self , x ):
		return x
##}}}

class InverseLinkFct(LinkFct): ##{{{
	"""
	SDFC.tools.InverseLinkFct
	=========================
	
	Inverse link function, i.e.:
		f(x) = 1/x
		f^{-1}(x) = 1/x
		df(x) = - 1 / x**2
	
	"""
	def __init__(self):
		pass
	
	def __str__(self):
		return "SDFC.tools.InverseLinkFct"
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return 1. / x
	
	def gradient( self , x ):
		return - 1. / x**2
	
	def inverse( self , x ):
		return 1. / x
##}}}

class ExpLinkFct(LinkFct):##{{{
	"""
	SDFC.tools.ExpLinkFct
	=====================
	
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
		return "SDFC.tools.ExpLinkFct, bounds:{}, sign:{}".format(self.b,self.s)
	
	def __repr__(self):
		return self.__str__()
	
	def __call__( self , x ):
		return np.exp(self.s * x) + self.b
	
	def gradient( self , x ):
		return self.s * np.exp(self.s * x)
	
	def inverse( self , x ):
		return np.log(x - self.b) / self.s
##}}}

class LogitLinkFct(LinkFct):##{{{
	"""
	SDFC.tools.LogitLinkFct
	=======================
	
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
		return "SDFC.tools.LogitLinkFct, bounds: ({},{}), speed: {}".format(self.a,self.b,self.s)
	
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

class SemiBoundedLinkFct(LinkFct):##{{{
	"""
	SDFC.tools.SemiBoundedLinkFct
	=============================
	A simple semi bounded function, use it if ExpLinkFct has too many overflow (Here the gradient is not well defined at x = b, and the inverse also).
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
		return "SDFC.tools.SemiBoundedLinkFct, bounds : {}, slope : {}".format(self.b,self.s)
	
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

class BoundedLinkFct(LinkFct):##{{{
	"""
	SDFC.tools.BoundedLinkFct
	=========================
	A simple bounded function, use it if LogitLinkFct has too many overflow (Here the gradient is not well defined at x = a and x = b, and the inverse also).
	Values:
		- x < a : f(x) = a
		- x > b : f(x) = b
		- a <= x <= b : f(x) = x
	"""
	def __init__( self , a , b ):
		self.a = a
		self.b = b
	
	def __str__(self):
		return "SDFC.tools.BoundedLinkFct, bounds : ({},{})".format(self.a,self.b)
	
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

