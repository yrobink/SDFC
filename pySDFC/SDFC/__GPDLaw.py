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

##################################################################################
## This software is based on extRemes library (R), see:                         ##
## https://cran.r-project.org/web/packages/extRemes/index.html                  ##
##################################################################################

##################################################################################
## Ce programme est basé sur la librairie R extRemes, voir:                     ##
## https://cran.r-project.org/web/packages/extRemes/index.html                  ##
##################################################################################


###############
## Libraries ##
###############

import numpy                  as np
import scipy.stats            as sc
import scipy.linalg           as scl
import scipy.optimize         as sco
import scipy.spatial.distance as ssd

from SDFC.__lmoments import lmoments


#############
## Classes ##
#############

class GPDLaw:
	"""
	SDFC.GPDLaw
	===========
	
	Fit generalized pareto law, possibly with co-variable
	
	"""
	
	def __init__( self , use_phi = False , method = "BFGS" , verbose = False ): ##{{{
		"""
		Initialization of GPDLaw
		
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
		shape        : numpy.ndarray
			Shape parameter(s)
		scale_design : numpy.ndarray
			Design matrix for scale
		shape_design : numpy.ndarray
			Design matrix for shape
		nscale       : integer
			Number of co-variate for scale + 1 (intercept)
		nshape       : integer
			Number of co-variate for shape + 1 (intercept)
		ncov         : integer
			nscale + nshape
		optim_result : scipy.optimize.OptimizeResult
			Result of minimization of likelihood
		scale_coef_  : numpy.ndarray
			coefficient fitted for scale
		shape_coef_  : numpy.ndarray
			coefficient fitted for shape
		
		"""
		self._use_phi = use_phi
		self._method  = method
		self.verbose = verbose
		
		self._Y           = None
		self._size        = None
		self.loc          = None
		self.scale        = None
		self.shape        = None
		self.scale_design = None
		self.shape_design = None
		self.nscale       = None
		self.nshape       = None
		self.ncov         = None
		self.optim_result = None
		self.scale_coef_  = None
		self.shape_coef_  = None
	##}}}
	
	def fit( self , Y , loc , scale_cov = None , shape_cov = None ): ##{{{
		"""
		Fit function for GPDLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		loc       : numpy.ndarray or float
			Loc parameter (use QuantileRegression to find the loc parameters)
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		shape_cov : None or numpy.ndarray
			Co-variates of shape in columns.
		"""
		self._Y    = Y.ravel()
		self._size = Y.size
		self.loc   = loc.ravel() if loc.size == self._size else np.repeat( loc[0] , self._size )
		self.scale = np.zeros(self._size)
		self.shape = np.zeros(self._size)
		
		## Design matrix
		scale_cov = scale_cov if ( scale_cov is None or scale_cov.ndim > 1 ) else scale_cov.reshape( (self._size,1) )
		shape_cov = shape_cov if ( shape_cov is None or shape_cov.ndim > 1 ) else shape_cov.reshape( (self._size,1) )
		self.scale_design = np.hstack( (np.ones( (self._size,1) ) , scale_cov) ) if scale_cov is not None else np.zeros( (self._size,1) ) + 1.
		self.shape_design = np.hstack( (np.ones( (self._size,1) ) , shape_cov) ) if shape_cov is not None else np.zeros( (self._size,1) ) + 1.
		
		if np.linalg.matrix_rank(self.scale_design) < self.scale_design.shape[1]:
			if self.verbose:
				print( "SDFC.GPDLaw: singular design matrix for scale, co-variable coefficients are set to 0" )
			self.scale_design = np.ones( (self._size,1) )
		
		if np.linalg.matrix_rank(self.shape_design) < self.shape_design.shape[1]:
			if self.verbose:
				print( "SDFC.GPDLaw: singular design matrix for shape, co-variable coefficients are set to 0" )
			self.shape_design = np.ones( (self._size,1) )
		
		self.nscale = self.scale_design.shape[1]
		self.nshape = self.shape_design.shape[1]
		self.ncov   = self.nscale + self.nshape
		
		## Initial condition
		init_scale,init_shape = self._find_init()
		param_init = np.zeros( (self.ncov) )
		param_init[0] = init_scale
		param_init[self.nscale] = init_shape
		
		## Optimization
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = self._method )
		
		## Set result
		self.scale_coef_ = self.optim_result.x[:self.nscale]
		self.shape_coef_ = self.optim_result.x[self.nscale:]
		self._update_param( self.optim_result.x )
	##}}}
	
	def _link( self , x ): ##{{{
		return np.exp(x) if self._use_phi else x
	##}}}
	
	def _link_inv( self , x ): ##{{{
		return np.log(x) if self._use_phi else x
	##}}}
	
	def _update_param( self , param ):##{{{
		## Extract coefficients from param
		scale_coef = param[:self.nscale]
		shape_coef = param[self.nscale:]
		
		## Set scale and shape
		self.scale = self._link( self.scale_design @ scale_coef )
		self.shape = self.shape_design @ shape_coef
	##}}}

	def _find_init( self ):##{{{
		## LMoments initial condition
		idx_excess = (self._Y > self.loc)
		excess = self._Y[idx_excess] - self.loc[idx_excess]
#		lmo1   = np.mean(excess)
#		lmo2   = np.sum( ssd.squareform( ssd.pdist( excess.reshape( (excess.size,1) ) , "cityblock" )) ) / ( 2 * excess.size * (excess.size - 1) )
		lmo1   = lmoments( excess , 1 )
		lmo2   = lmoments( excess , 2 )
		itau     = lmo1 / lmo2
		scale_lm = lmo1 * ( itau - 1 )
		scale_lm = scale_lm if scale_lm > 0 else 1e-8
		shape_lm = - ( itau - 2 )
		
		## Check negloglikelihood
		self.scale = np.repeat( scale_lm , self._size )
		self.shape = np.repeat( shape_lm , self._size )
		eval_lm = self._negloglikelihood()
		
		## MOMS initial condition
		scale_mm = np.sqrt( np.var(excess) )
		scale_mm = scale_mm if scale_mm > 0 else 1e-8
		shape_mm = -1e-8
		
		## Check negloglikelihood
		self.scale = np.repeat( scale_mm , self._size )
		self.shape = np.repeat( shape_mm , self._size )
		eval_mm = self._negloglikelihood()
		
		## Keep best
		init_scale,init_shape = self._link_inv(1.),1e-2
		if np.isfinite(eval_mm) or np.isfinite(eval_lm):
			init_scale = self._link_inv( scale_mm if eval_mm < eval_lm else scale_lm )
			init_shape = shape_mm if eval_mm < eval_lm else shape_lm
		return init_scale,init_shape
	##}}}
	
	def _negloglikelihood( self ): ##{{{
		## Impossible scale
		if np.any( self.scale <= 0 ):
			return np.inf
		
		## Fuck exponential case
		zero_shape = ( np.abs(self.shape) < 1e-10 )
		if np.any(zero_shape):
			self.shape[zero_shape] = -1e-10
		
		##
		idx_excess = (self._Y > self.loc)
		loc   = self.loc[idx_excess]
		scale = self.scale[idx_excess]
		shape = self.shape[idx_excess]
		Z = 1. + shape * ( self._Y[idx_excess] - loc ) / scale
		
		if np.any(Z <= 0):
			return np.inf
		
		res = np.sum( np.log( scale ) + np.log(Z) * ( 1 + 1. / shape ) )
		
		return res if np.isfinite(res) else np.inf
	##}}}
	
	def _optim_function( self , param ):##{{{
		self._update_param(param)
		return self._negloglikelihood()
	##}}}
	
	def _gradient_optim_function( self , param ): ##{{{
		
		self._update_param(param)
		
		idx_excess = ( self._Y > self.loc )
		Y      = self._Y[idx_excess]
		loc    = self.loc[idx_excess]
		scale  = self.scale[idx_excess]
		shape  = self.shape[idx_excess]
		Y_zero   = Y - loc
		Z        = 1. + shape * Y_zero / scale
		exponent = 1. + 1. / shape
		
		grad_phi   = scale if self._use_phi else 1.
		grad_scale = self.scale_design[idx_excess,:].transpose() @ (grad_phi / scale - ( grad_phi * exponent * shape * Y_zero / (scale**2) ) / Z)
		grad_shape = self.shape_design[idx_excess,:].transpose() @ ( - np.log(Z) / (shape**2) + exponent * Y_zero / scale / Z) if np.all(Z>0) else np.repeat(np.nan,self.nshape)
		
#		phi_grad = scale if self._use_phi else 1.
#		coef = phi_grad / scale - ( phi_grad * exponent * shape * Y_zero / (scale**2) ) / Z
#		grad_scale = np.apply_along_axis( lambda X : np.sum(coef * X) , 0 , self.scale_design[idx_excess,:] )
#		if np.all(Z > 0):
#			coef = - np.log(Z) / (shape**2) + exponent * Y_zero / scale / Z
#		else:
#			coef = np.nan
#		grad_shape = np.apply_along_axis( lambda X : np.sum(coef * X) , 0 , self.shape_design[idx_excess,:] )
		
		
		return np.hstack( (grad_scale,grad_shape) )
	##}}}

