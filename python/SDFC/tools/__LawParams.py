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
from SDFC.tools.__Link import IdLink


###########
## Class ##
###########

class AbstractParam:##{{{
	def __init__( self , kind , n_samples , **kwargs ):
		self.kind       = kind
		self.link    = IdLink() if kwargs.get("l_" + self.kind) is None else kwargs.get("l_" + self.kind)
		self.n_samples  = n_samples
		self.n_features = 0
		self.coef_      = None
		self.fit_       = None
	
	def is_fix(self):
		pass
	
	@property
	def value( self ):
		return self.link(self.fit_)
	
	def gradient( self ):
		return self.link.gradient(self.fit_)
	
	def set_coef( self , coef ):
		pass
	
##}}}

class CovariateParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , resample , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.design_ = None
		
		## Build design matrix
		X = kwargs.get( "c_" + self.kind )
		if X.ndim == 1: X = X.reshape(-1,1)
		self.n_features = X.shape[1] + 1
		if resample is None:
			self.design_ = np.hstack( ( np.ones((self.n_samples,1)) , X ) )
		else:
			self.design_ = np.hstack( ( np.ones((self.n_samples,1)) , X[resample,:] ) )
		self.coef_   = np.zeros(self.n_features)
		
		if np.linalg.matrix_rank(self.design_) < self.n_features:
			print( "SDFC.LawParam: singular design matrix" )
			self.coef_   = np.zeros(1)
			self.design_ = np.ones((self.n_samples,1))
		
		self.update()
	
	def is_fix(self):
		return False
	
	def update( self ):
		self.fit_ = (self.design_ @ self.coef_).reshape(-1,1)
	
	def set_coef( self , coef ):
		self.coef_ = coef.squeeze()
		self.update()
	
	def design_wo1(self):
		return self.design_[:,1:]
##}}}

class StationaryParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , resample , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.coef_      = np.zeros(1)
		self.n_features = 1
		self.fit_       = np.zeros(self.n_samples)
		self.design_    = np.ones( (self.n_samples,1) )
	
	def is_fix(self):
		return False
	
	def update( self ):
		self.fit_ = np.repeat( self.coef_ , self.n_samples ).reshape(-1,1)
	
	def set_coef( self , coef ):
		self.coef_ = np.array([coef]).squeeze()
		self.update()
	
	def design_wo1(self):
		return None
##}}}

class FixParam(AbstractParam):##{{{
	def __init__( self , kind , n_samples , resample , **kwargs ):
		AbstractParam.__init__( self , kind , n_samples , **kwargs )
		self.fit_       = np.array( [self.link.inverse( kwargs.get( "f_" + self.kind ) )] ).reshape(-1,1)
		if self.fit_.size == 1: self.fit_ = np.repeat( self.fit_ , self.n_samples ).reshape(-1,1)
	
	def is_fix(self):
		return True
	
##}}}

class LawParams:##{{{
	
	def __init__( self , kinds , **kwargs ):##{{{
		self.kinds = kinds
		self._dparams = {}
		self.coef_ = None
	##}}}
	
	def add_params( self , n_samples , resample , **kwargs ):##{{{
		for kind in self.kinds:
			k_param = self.filter( kind , **kwargs )
			config = self.infer_configuration(**k_param)
			if self.is_covariate(config):  self._dparams[kind] = CovariateParam(  kind , n_samples , resample , **k_param )
			if self.is_stationary(config): self._dparams[kind] = StationaryParam( kind , n_samples , resample , **k_param )
			if self.is_fix(config):        self._dparams[kind] = FixParam(        kind , n_samples , resample , **k_param )
	##}}}
	
	def merge_coef( self ):##{{{
		self.coef_ = np.array([])
		for k in self._dparams:
			if not self._dparams[k].is_fix():
				self.coef_ = np.hstack( (self.coef_,self._dparams[k].coef_.squeeze()) )
	##}}}
	
	def split_coef( self , coef ):##{{{
		tcoef = []
		a,b = 0,0
		for k in self._dparams:
			if self._dparams[k].is_fix():
				tcoef.append(None)
			else:
				b = a + self._dparams[k].n_features
				tcoef.append( coef[a:b] )
				a = b
		return tcoef
	##}}}
	
	def update_coef( self , coef , kind = None ):##{{{
		if kind is None:
			lcoef = self.split_coef(coef)
			for c,k in zip(lcoef,self._dparams):
				self._dparams[k].set_coef(c)
		else:
			self._dparams[kind].set_coef(coef)
		self.merge_coef()
	##}}}
	
	def infer_configuration( self , **kwargs ):##{{{
		has_c = False
		has_s = False
		has_f = False
		has_l = False
		for k in kwargs:
			if k[:2] == "c_": has_c = True
			if k[:2] == "f_": has_f = True
			if k[:2] == "l_": has_l = True
		
		if len(kwargs) == 0 or ( len(kwargs) == 1 and has_l):
			has_s = True
		return has_c,has_s,has_f,has_l
	
	def is_covariate( self , c ):
		return c[0] and not ( c[1] or c[2] )
	
	def is_stationary( self , c ):
		return c[1] and not ( c[0] or c[2] )
	
	def is_fix( self , c ):
		return c[2] and not ( c[0] or c[1] )
	
	def filter( self , kind , **kwargs ):
		out = {}
		for k in kwargs:
			if kind in k:
				out[k] = kwargs[k]
		return out
	##}}}

##}}}

