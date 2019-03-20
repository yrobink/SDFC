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

import numpy          as np
import scipy.special  as scs
import scipy.linalg   as scl
import scipy.optimize as sco

from SDFC.__lmoments import lmoments
from SDFC.__QuantileRegression import QuantileRegression

#############
## Classes ##
#############


class GEVLaw:##{{{
	"""
	SDFC.GEVLaw
	===========
	
	Fit generalized extreme value law, possibly with co-variable
	
	"""
	
	def __init__( self , use_phi = False , method = "BFGS" , verbose = False ): ##{{{
		"""
		Initialization of GEVLaw
		
		Parameters
		----------
		use_phi     : bool
			If True, the exponential link function is used to fit the the scale parameter, default False
		method      : string
			Method called to minimize the negloglikelihood function, default "BFGS"
		verbose     : bool
			If True, warning and error are printed
		
		Attributes
		----------
		
		loc          : numpy.ndarray
			Location parameter(s)
		scale        : numpy.ndarray
			Scale parameter(s)
		shape        : numpy.ndarray
			Shape parameter(s)
		loc_design   : numpy.ndarray
			Design matrix for loc
		scale_design : numpy.ndarray
			Design matrix for scale
		shape_design : numpy.ndarray
			Design matrix for shape
		nloc         : integer
			Number of co-variate for loc + 1 (intercept)
		nscale       : integer
			Number of co-variate for scale + 1 (intercept)
		nshape       : integer
			Number of co-variate for shape + 1 (intercept)
		ncov         : integer
			nloc + nscale + nshape
		optim_result : scipy.optimize.OptimizeResult
			Result of minimization of likelihood
		scale_coef_  : numpy.ndarray
			coefficient fitted for scale
		shape_coef_  : numpy.ndarray
			coefficient fitted for shape
		
		Notes
		-----
		
		To find first estimates, two methods are available. "extRemes" use the same method as extRemes R package. "quantiles"
		use quantile regression. For example, the "loc" parameter is the quantile exp(-1), and the scale depend linearly of the
		centered quantile function.
		
		
		"""
		self._use_phi     = use_phi
		self._method      = method
		self.verbose      = verbose
		
		self._Y           = None
		self._size        = None
		self.loc          = None
		self.scale        = None
		self.shape        = None
		self.loc_design   = None
		self.scale_design = None
		self.shape_design = None
		self.nloc         = None
		self.nscale       = None
		self.nshape       = None
		self.ncov         = None
		self.optim_result = None
		self.loc_coef_    = None
		self.scale_coef_  = None
		self.shape_coef_  = None
	##}}}
	
	def fit( self , Y , loc_cov = None , scale_cov = None , shape_cov = None ): ##{{{
		"""
		Fit function for GEVLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		loc_cov   : None or numpy.ndarray
			Co-variates of loc in columns.
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		shape_cov : None or numpy.ndarray
			Co-variates of shape in columns.
		"""
		self._Y    = Y.ravel()
		self._size = Y.size
		self.loc   = np.zeros(self._size)
		self.scale = np.zeros(self._size)
		self.shape = np.zeros(self._size)
		
		## Design matrix
		loc_cov   = loc_cov   if ( loc_cov   is None or loc_cov.ndim > 1 )   else loc_cov.reshape( (self._size,1) )
		scale_cov = scale_cov if ( scale_cov is None or scale_cov.ndim > 1 ) else scale_cov.reshape( (self._size,1) )
		shape_cov = shape_cov if ( shape_cov is None or shape_cov.ndim > 1 ) else shape_cov.reshape( (self._size,1) )
		self.loc_design   = np.hstack( (np.ones( (self._size,1) ) , loc_cov) )   if loc_cov   is not None else np.zeros( (self._size,1) ) + 1.
		self.scale_design = np.hstack( (np.ones( (self._size,1) ) , scale_cov) ) if scale_cov is not None else np.zeros( (self._size,1) ) + 1.
		self.shape_design = np.hstack( (np.ones( (self._size,1) ) , shape_cov) ) if shape_cov is not None else np.zeros( (self._size,1) ) + 1.
		
		if np.linalg.matrix_rank(self.loc_design) < self.loc_design.shape[1]:
			if self.verbose:
				print( "SDFC.GEVLaw: singular design matrix for loc, co-variable coefficients are set to 0" )
			self.loc_design = np.ones( (self._size,1) )
		
		if np.linalg.matrix_rank(self.scale_design) < self.scale_design.shape[1]:
			if self.verbose:
				print( "SDFC.GEVLaw: singular design matrix for scale, co-variable coefficients are set to 0" )
			self.scale_design = np.ones( (self._size,1) )
		
		if np.linalg.matrix_rank(self.shape_design) < self.shape_design.shape[1]:
			if self.verbose:
				print( "SDFC.GEVLaw: singular design matrix for shape, co-variable coefficients are set to 0" )
			self.shape_design = np.ones( (self._size,1) )
		
		self.nloc   = self.loc_design.shape[1]
		self.nscale = self.scale_design.shape[1]
		self.nshape = self.shape_design.shape[1]
		self.ncov   = self.nloc + self.nscale + self.nshape
		
		## Initial condition
		param_init = self._find_init()
		
		## Optimization
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = self._method )
		
		## Set result
		self.loc_coef_   = self.optim_result.x[:self.nloc]
		self.scale_coef_ = self.optim_result.x[self.nloc:(self.nloc+self.nscale)]
		self.shape_coef_ = self.optim_result.x[(self.nloc+self.nscale):]
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
		loc_coef   = param[:self.nloc]
		scale_coef = param[self.nloc:(self.nloc+self.nscale)]
		shape_coef = param[(self.nloc+self.nscale):]
		
		## Set scale and shape
		self.loc   = self.loc_design @ loc_coef
		self.scale = self._link( self.scale_design @ scale_coef )
		self.shape = self.shape_design @ shape_coef
	##}}}
	
	def _find_init_extRemes( self ): ##{{{
		
		## LMoments
		lmom1 = lmoments( self._Y , 1 )
		lmom2 = lmoments( self._Y , 2 )
		lmom3 = lmoments( self._Y , 3 )
		
		tau3  = lmom3 / lmom2
		co    = 2. / ( 3. + tau3 ) - np.log(2) / np.log(3)
		kappa = 7.8590 * co + 2.9554 * co**2
		g     = scs.gamma( 1. + kappa )
		
		init_loc   = np.zeros( self.nloc )
		init_scale = np.zeros( self.nscale )
		init_shape = np.zeros( self.nshape )
		
		init_scale[0] = lmom2 * kappa / ( (1 - np.power( 2 , - kappa )) * g )
		init_loc[0]   = lmom1 - init_scale[0] * (1 - g) / kappa
		init_shape[0] = - kappa
		init_scale[0] = self._link_inv(init_scale[0])
		
		
		param   = np.hstack( (init_loc,init_scale,init_shape) )
		test_ll = self._optim_function(param)
		
		## Moments
		m     = np.mean(self._Y)
		s     = np.sqrt(6) * np.std(self._Y) / np.pi
		
		init_loc2   = np.zeros( self.nloc )
		init_scale2 = np.zeros( self.nscale )
		init_shape2 = np.zeros( self.nshape )
		
		init_loc2[0]   = m - 0.57722 * s
		init_scale2[0] = np.log(s)
		init_shape2[0] = 1e-8
		init_scale2[0] = self._link_inv(init_scale2[0])
		
		param2   = np.hstack( (init_loc2,init_scale2,init_shape2) )
		test_ll2 = self._optim_function(param2)
		
		if not np.isfinite(test_ll) and not np.isfinite(test_ll2):
			init_loc3   = np.zeros( self.nloc )
			init_scale3 = np.zeros( self.nscale )
			init_shape3 = np.zeros( self.nshape )
			init_loc3[0]   = 0
			init_scale3[0] = 1
			init_shape3[0] = 0.1
			param3 = np.hstack( (init_loc3,init_scale3,init_shape3) )
			test_ll3 = self._optim_function(param3)
			return param3,test_ll3
		
		return (param,test_ll) if test_ll < test_ll2 else (param2,test_ll2)
	##}}}
	
	def _find_init_quantiles( self ): ##{{{
		
		init_loc   = np.zeros( self.nloc )
		init_scale = np.zeros( self.nscale )
		init_shape = np.zeros( self.nshape )
		
		## Fit loc
		loc = None
		if self.nloc == 1:
			init_loc[0] = sc.rv_histogram( np.histogram( self._Y , 100 ) ).ppf( np.exp(-1) )
			loc = np.repeat( init_loc[0] , self._Y.size )
		else:
			reg = QuantileRegression( [np.exp(-1)] )
			reg.fit( self._Y , self.loc_design[:,1:] )
			init_loc = reg.coef_.ravel()
			loc = reg.predict().ravel()
		
		## Fit shape	
		lmom1 = lmoments( self._Y - loc , 1 )
		lmom2 = lmoments( self._Y - loc , 2 )
		lmom3 = lmoments( self._Y - loc , 3 )
		
		tau3  = lmom3 / lmom2
		co    = 2. / ( 3. + tau3 ) - np.log(2) / np.log(3)
		kappa = 7.8590 * co + 2.9554 * co**2
		init_shape[0] = - kappa
		
		## Fit scale
		qscale = np.array([0.25,0.5,0.75])
		if self.nscale == 1:
			coef = - kappa / ( np.power( -np.log(qscale) , kappa ) - 1 )
			init_scale[0] = self._link_inv(np.mean( sc.rv_histogram( np.histogram( self._Y - loc , 100 ) ).ppf( qscale ) * coef ))
		else:
			reg = QuantileRegression( qscale )
			reg.fit( self._Y - loc , self.scale_design[:,1:] )
			fshape = np.repeat( -kappa , self._Y.size )
			coef = np.array( [ fshape / ( np.power( -np.log(x) , - fshape ) - 1 ) for x in qscale ] ).T
			fscale = np.mean( reg.predict() * coef , axis = 1 )
			init_scale,_,_,_ = scl.lstsq( self.scale_design , self._link_inv(fscale) )
		
		
		param   = np.hstack( (init_loc,init_scale,init_shape) )
		while not np.isfinite(self._optim_function(param)) or not np.all(np.isfinite(self._gradient_optim_function(param))):
			param[(self.nloc+self.nscale)] *= 0.95
		
		test_ll = self._optim_function(param)
		
		return param,test_ll
	##}}}
	
	def _find_init( self ):##{{{
		
		param_e,ll_e = self._find_init_extRemes()
		param_q,ll_q = self._find_init_quantiles()
		
		return param_e if ll_e < ll_q else param_q
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
		Z = 1 + self.shape * ( self._Y - self.loc ) / self.scale
		
		if np.any(Z <= 0):
			return np.inf
		
		res = np.sum( ( 1. + 1. / self.shape ) * np.log(Z) + np.power( Z , - 1. / self.shape ) + np.log(self.scale) )
		
		
		return res if np.isfinite(res) else np.inf
	##}}}
	
	def _optim_function( self , param ):##{{{
		self._update_param(param)
		return self._negloglikelihood()
	##}}}
	
	def _logZafun( self, Z , alpha ):##{{{
		return alpha * np.log( 1. + self.shape * Z )
	##}}}
	
	def _Zafun( self , Z , alpha ):##{{{
		return np.exp( self._logZafun( Z , alpha ) )
	##}}}
	
	def _gradient_optim_function( self , param ): ##{{{
		
		self._update_param(param)
		
		## Impossible
		if np.any( 1. + self.shape * ( self._Y - self.loc ) / self.scale <= 0 ):
			return np.zeros( self.ncov ) + np.nan
		
		## Usefull values
		Z      = ( self._Y - self.loc ) / self.scale
		Za1    = self._Zafun( Z , 1. )
		ishape = 1. / self.shape
		Zamsi  = self._Zafun( Z , - ishape ) ## Za of Minus Shape Inverse
		
		## Vectors
		phi_vect   = 1. if self._use_phi else 1. / self.scale
		loc_vect   = (Zamsi - 1 - self.shape) / ( self.scale * Za1 )
		scale_vect = phi_vect * ( 1. + Z * ( Zamsi - 1 - self.shape ) / Za1 )
		shape_vect = ( Zamsi - 1. ) * np.log(Za1) * ishape**2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1
		
		## Gradients
		grad_loc   = np.dot( self.loc_design.transpose()   , loc_vect   )
		grad_scale = np.dot( self.scale_design.transpose() , scale_vect )
		grad_shape = np.dot( self.shape_design.transpose() , shape_vect )
		
		return np.hstack( (grad_loc,grad_scale,grad_shape) )
	##}}}
##}}}




