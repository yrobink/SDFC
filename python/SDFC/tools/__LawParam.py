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
from SDFC.tools.__LinkFct import IdLinkFct


###########
## Class ##
###########

class LawParam:##{{{
	"""
	Internal class, do not use
	"""
	def __init__( self , linkFct = IdLinkFct() , kind = "Parameter" ):##{{{
		self.linkFct    = linkFct
		self.coef_      = None
		self.design_    = None
		self.size       = 0
		self.dim        = 0
		self.shape      = None
		self._transform = None
		self._value     = None
		self._not_fixed = None
		self.kind       = kind
	##}}}
	
	def copy(self):##{{{
		p            = LawParam( self.linkFct , self.kind )
		p.coef_      = np.copy(self.coef_)   if self.coef_   is not None else None
		p.design_    = np.copy(self.design_) if self.design_ is not None else None
		p.size       = self.size
		p.dim        = self.dim
		p.shape      = self.shape
		p._transform = self._transform
		p._value     = np.copy(self._value)  if self._value  is not None else None
		p._not_fixed = self._not_fixed
		return p
	##}}}
	
	def init( self , X = None , fix_values = None , size = None , dim = 1 , transform = lambda x : x ): ##{{{
		self._not_fixed = fix_values is None
		self.dim        = dim
		self._transform = transform
		if fix_values is not None:
			if dim == 1:
				fix_values = np.array( [fix_values] ).reshape(-1,dim)
				if fix_values.size == 1:
					fix_values = np.repeat( fix_values[0,0] , X.size ).reshape(-1,dim)
			self._value = self.linkFct.inverse(fix_values)
		else:
			if X is None and size is not None:
				self.coef_   = np.zeros( (1,dim) )
				self.design_ = np.ones( (size,1) )
			elif X is not None:
				size = X.shape[0]
				if X.ndim == 1:
					X = X.reshape( (size,1) )
				self.design_ = np.hstack( ( np.ones( (size,1) ) , X ) )
				self.coef_   = np.zeros( (X.shape[1]+1,dim) )
			else:
				print( "Error" )
			
			if np.linalg.matrix_rank(self.design_) < self.design_.shape[1]:
				print( "SDFC.LawParam: singular design matrix" )
				self.coef_   = np.zeros( (1,dim) )
				self.design = np.ones( (size,1) )
			self.size  = np.prod(self.coef_.shape)
			self.shape = self.coef_.shape
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.LawParam\n"
		val += "-------------\n"
		val += "* kind     : {}\n".format(self.kind)
		val += "* fixed    : {}\n".format(not self._not_fixed)
		val += "* link_fct : {}\n".format(str(self.linkFct))
		val += "* coef_    : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def not_fixed(self):##{{{
		return self._not_fixed
	##}}}
	
	def set_coef( self , coef ):##{{{
		if self._not_fixed:
			self.coef_ = coef.reshape( self.shape )
	##}}}
	
	def set_intercept( self , intercept ):##{{{
		if self._not_fixed:
			self.coef_[0,:] = intercept
	##}}}
	
	def design_wo1( self ):##{{{
		"""
		Return design matrix without intercept
		"""
		return None if self.shape[0] == 1 else self.design_[:,1:]
	##}}}
	
	def update(self):##{{{
		if self._not_fixed:
			self._value = self._transform( self.design_ @ self.coef_ )
	##}}}
	
	def value( self ):##{{{
		return self._value
	##}}}
	
	def valueLf( self ):##{{{
		out = self.linkFct( self._value )
		return out
	##}}}
	
	def valueGrLf( self ):##{{{
		return self.linkFct.gradient( self._value )
	##}}}

##}}}


class LawParam_deprecated:##{{{
	"""
	Internal class, do not use
	"""
	def __init__( self , linkFct = IdLinkFct() , kind = "Parameter" ):##{{{
		self.linkFct    = linkFct
		self.coef_      = None
		self.design_    = None
		self.size       = 0
		self._value     = None
		self._not_fixed = None
		self.kind       = kind
	##}}}
	
	def copy(self):##{{{
		p = LawParam( self.linkFct , self.kind )
		p.coef_ = np.copy(self.coef_) if self.coef_ is not None else None
		p.design_ = np.copy(self.design_) if self.design_ is not None else None
		p.size = self.size
		p._value = np.copy(self._value) if self._value is not None else None
		p._not_fixed = self._not_fixed
		return p
	##}}}
	
	def init( self , X = None , fix_values = None , size = None ): ##{{{
		self._not_fixed = fix_values is None
		if fix_values is not None:
			fix_values = np.array( [fix_values] ).ravel()
			if fix_values.size == 1 and size is not None:
				self._value = self.linkFct.inverse( np.repeat( fix_values , size ).ravel() )
			elif fix_values.size > 1:
				self._value = self.linkFct.inverse( fix_values )
			else:
				print( "Error" )
		else:
			if X is None and size is not None:
				self.coef_   = np.zeros( (1,) )
				self.design_ = np.ones( (size,1) )
			elif X is not None:
				size = X.shape[0]
				if X.ndim == 1:
					X = X.reshape( (size,1) )
				self.design_ = np.hstack( ( np.ones( (size,1) ) , X ) )
				self.coef_   = np.zeros( (X.shape[1]+1) )
			else:
				print( "Error" )
			
			if np.linalg.matrix_rank(self.design_) < self.design_.shape[1]:
				print( "SDFC.LawParam: singular design matrix" )
				self.coef_   = np.zeros( (1,) )
				self.design = np.ones( (size,1) )
			self.size = self.coef_.size
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.LawParam\n"
		val += "-------------\n"
		val += "* kind     : {}\n".format(self.kind)
		val += "* fixed    : {}\n".format(not self._not_fixed)
		val += "* link_fct : {}\n".format(str(self.linkFct))
		val += "* coef_    : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def not_fixed(self):##{{{
		return self._not_fixed
	##}}}
	
	def set_coef( self , coef ):##{{{
		if self._not_fixed:
			self.coef_ = coef.ravel()
	##}}}
	
	def set_intercept( self , intercept ):##{{{
		if self._not_fixed:
			self.coef_[0] = intercept
	##}}}
	
	def design_wo1( self ):##{{{
		"""
		Return design matrix without intercept
		"""
		return None if self.size == 1 else self.design_[:,1:]
	##}}}
	
	def update(self):##{{{
		if self._not_fixed:
#			print("LawParam.design.shape={}".format(self.design_.shape))
#			print("LawParam.coef.shape={}".format(self.coef_.shape))
#			print("LawParam.coef={}".format(self.coef_))
			self._value = np.ravel( self.design_ @ self.coef_ )
	##}}}
	
	def value( self ):##{{{
		return np.ravel( self._value )
	##}}}
	
	def valueLf( self ):##{{{
		return np.ravel( self.linkFct( self._value ) )
	##}}}
	
	def valueGrLf( self ):##{{{
		return self.linkFct.gradient( self._value )
	##}}}

##}}}


