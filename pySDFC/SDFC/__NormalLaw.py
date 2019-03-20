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

import numpy          as np
import scipy.stats    as sc
import scipy.linalg   as scl
import scipy.optimize as sco


#############
## Classes ##
#############

class NormalLaw:
	"""
	SDFC.NormalLaw
	==============
	
	Fit Normal law, possibly with co-variable
	
	"""
	
	def __init__( self , use_phi = False , method = "BFGS" , verbose = False ): ##{{{
		"""
		Initialization of NormalLaw
		
		Parameters
		----------
		use_phi : bool
			If True, the exponential link function is used to fit the the scale parameter, default False
		method  : string
			Method called to minimize the negloglikelihood function, default "BFGS"
		verbose : bool
			If True, warning and error are printed
		
		Attributes
		----------
		
		loc          : numpy.ndarray
			Location parameter(s)
		scale        : numpy.ndarray
			Scale parameter(s)
		loc_design   : numpy.ndarray
			Design matrix for loc
		scale_design : numpy.ndarray
			Design matrix for scale
		nloc         : integer
			Number of co-variate for loc + 1 (intercept)
		nscale       : integer
			Number of co-variate for scale + 1 (intercept)
		ncov         : integer
			nloc + nscale
		optim_result : scipy.optimize.OptimizeResult
			Result of minimization of likelihood
		loc_coef_    : numpy.ndarray
			coefficient fitted for loc
		scale_coef_  : numpy.ndarray
			coefficient fitted for scale
		
		"""
		self._use_phi = use_phi
		self._method  = method
		self.verbose = verbose
		
		self._Y            = None
		self._size         = None
		self._Nlog2pi      = None
		self.loc          = None
		self.scale        = None
		self.loc_design   = None
		self.scale_design = None
		self.nloc         = None
		self.nscale       = None
		self.ncov         = None
		self.optim_result = None
		self.loc_coef_    = None
		self.scale_coef_  = None
	##}}}
	
	def fit( self , Y , loc_cov = None , scale_cov = None ): ##{{{
		"""
		Fit function for NormalLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		loc_cov   : None or numpy.ndarray
			Co-variates of loc in columns.
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		"""
		self._Y    = Y.ravel()
		self._size = Y.size
		self._Nlog2pi = self._size * np.log( 2 * np.pi ) / 2.
		self.loc   = np.zeros(self._size)
		self.scale = np.zeros(self._size)
		
		## Design matrix
		loc_cov   = loc_cov   if ( loc_cov   is None or loc_cov.ndim > 1 )   else loc_cov.reshape( (self._size,1) )
		scale_cov = scale_cov if ( scale_cov is None or scale_cov.ndim > 1 ) else scale_cov.reshape( (self._size,1) )
		self.loc_design   = np.hstack( (np.ones( (self._size,1) ) , loc_cov) )   if loc_cov   is not None else np.zeros( (self._size,1) ) + 1.
		self.scale_design = np.hstack( (np.ones( (self._size,1) ) , scale_cov) ) if scale_cov is not None else np.zeros( (self._size,1) ) + 1.
		
		if np.linalg.matrix_rank(self.loc_design) < self.loc_design.shape[1]:
			if self.verbose:
				print( "SFDC.NormalLaw: singular design matrix for loc, co-variable coefficients are set to 0" )
			self.loc_design = np.ones( (self._size,1) )
		
		if np.linalg.matrix_rank(self.scale_design) < self.scale_design.shape[1]:
			if self.verbose:
				print( "SFDC.NormalLaw: singular design matrix for scale, co-variable coefficients are set to 0" )
			self.scale_design = np.ones( (self._size,1) )
		
		
		self.nloc   = self.loc_design.shape[1]
		self.nscale = self.scale_design.shape[1]
		self.ncov   = self.nloc + self.nscale
		
		## Initial condition
		init_loc,init_scale = self._find_init()
		param_init = np.zeros( (self.ncov) )
		param_init[:self.nloc] = init_loc
		param_init[self.nloc:] = init_scale
		
		## Optimization
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = self._method )
		
		## Set result
		self.loc_coef_   = self.optim_result.x[:self.nloc]
		self.scale_coef_ = self.optim_result.x[self.nloc:]
		self._update_param( self.optim_result.x )
	##}}}
	
	def _link( self , x ): ##{{{
		return np.exp(x) if self._use_phi else x
	##}}}
	
	def _link_inv( self , x ): ##{{{
		return np.log(x) if self._use_phi else x
	##}}}
	
	def _update_param( self , param ):##{{{
		self.loc   = np.dot( self.loc_design   , param[:self.nloc] )
		self.scale = self._link( np.dot( self.scale_design , param[self.nloc:] ) )
	##}}}

	def _find_init( self ):##{{{
		
		## Estimate location
		init_loc,_,_,_ = scl.lstsq( self.loc_design , self._Y )
		
		## Estimate scale
		init_scale = np.zeros( self.nscale )
		init_scale[0] = self._link_inv(np.std(self._Y))
		
		return init_loc,init_scale
	##}}}
	
	def _negloglikelihood( self ): ##{{{
		scale2 = np.power( self.scale , 2 )
		return np.Inf if not np.all( self.scale > 0 ) else self._Nlog2pi + np.sum( np.log( scale2 ) ) / 2. + np.sum( np.power( self._Y - self.loc , 2 ) / scale2 ) / 2.
	##}}}
	
	def _optim_function( self , param ):##{{{
		self._update_param(param)
		return self._negloglikelihood()
	##}}}
	
	def _gradient_optim_function( self , param ): ##{{{
		self._update_param(param)
		
		Yc = self._Y - self.loc
		grad_loc   = - self.loc_design.transpose() @ (Yc / self.scale**2)
		grad_phi   = 1. if self._use_phi else self.scale
		grad_scale = self.scale_design.transpose() @ ( 1. / grad_phi - (Yc / self.scale)**2 / grad_phi )
		
#		grad_scale = np.zeros(self.nscale)
#		if not self._use_phi:
#			grad_scale = self.scale_design.transpose() @ ( 1. / self.scale - (Yc**2 / self.scale**3) )
#		else:
#			for i in range(self.nscale):
#				grad_scale[i] = np.sum(self.scale_design[:,i] - self.scale_design[:,i] * (Yc / self.scale)**2)
#		grad_loc   = np.zeros( self.nloc   )
#		grad_scale = np.zeros( self.nscale )
#		grad_loc[0] = - np.sum( Yc / self.scale**2 )
#		for i in range(self.nloc):
#			grad_loc[i] = - np.sum( self.loc_design[:,i] * Yc / self.scale**2 )
		
#		grad_scale[0] = np.sum( 1. / self.scale ) - np.sum( Yc**2 / self.scale**3 )
#		for i in range(self.nscale):
#			grad_scale[i] = np.sum( self.scale_design[:,i] / self.scale ) - np.sum( self.scale_design[:,i] * Yc**2 / self.scale**3 )
		
		
		return np.hstack( (grad_loc,grad_scale) )
	##}}}

