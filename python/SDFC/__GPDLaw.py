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

import numpy                  as np
import scipy.optimize         as sco

from SDFC.__AbstractLaw            import AbstractLaw
from SDFC.NonParametric.__mean     import mean
from SDFC.NonParametric.__lmoments import lmoments
from SDFC.NonParametric.__std      import std
from SDFC.tools.__LawParam         import LawParam
from SDFC.tools.__LinkFct          import IdLinkFct


###########
## Class ##
###########

class GPDLaw(AbstractLaw):
	"""
	SDFC.GPDLaw
	===========
	
	Fit parameters of a Generalized Pareto Distribution, possibly with co-variable
	
	Attributes
	----------
	
	method : string
		method used to fit
	loc    : numpy.ndarray
		Location, given by user in fit function
	scale  : numpy.ndarray
		Scale fitted
	shape  : numpy.ndarray
		Shape fitted
	coef_  : numpy.ndarray
		Coefficients fitted
	coefs_bootstrap: numpy.ndarray
		coef_ for each bootstrap
	confidence interval: numpy.ndarray[ shape = (2,coef_.size) ]
		Confidence interval, first line is the alpha/2 quantile, and second line the 1 - alpha/2 quantile
	alpha          : float
		Level of confidence interval
	"""
	
	def __init__( self , method = "MLE" , link_fct_scale = IdLinkFct() , link_fct_shape = IdLinkFct() , n_bootstrap = 0 , alpha = 0.05 ):##{{{
		"""
		Initialization of GPDLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments", "lmoments" and "MLE" (Maximum Likelihood estimation)
		link_fct_scale : a class herited from SDFC.tools.LinkFct
			Link function for scale, default is SDFC.tools.IdLinkFct(). Interesting option is SDFC.tools.ExpLinkFct().
		link_fct_shape : a class herited from SDFC.tools.LinkFct
			Link function for shape, default is IdLinkFct(). Interesting option is SDFC.tools.LogitLinkFct( -0.5 , 0.5 ).
		"""
		AbstractLaw.__init__( self , method , n_bootstrap , alpha )
		
		self.loc       = None
		
		self.scale     = None
		self.shape     = None
		
		self._scale = LawParam( linkFct = link_fct_scale , kind = "scale" )
		self._shape = LawParam( linkFct = link_fct_shape , kind = "shape" )
		
		self._lparams = [self._scale,self._shape]
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.GPDLaw\n"
		val += "-----------\n"
		val += "* fit_method : {}\n".format(self.method)
		val += "* link_scale : {}\n".format(str(self._scale.linkFct))
		val += "* link_shape : {}\n".format(str(self._shape.linkFct))
		val += "* coef_      : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def fit( self , Y , loc , scale_cov = None , shape_cov = None , fscale = None , fshape = None ):##{{{
		"""
		Fit function for GPDLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		loc       : numpy.ndarray
			Location, i.e. threshold
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		shape_cov : None or numpy.ndarray
			Co-variates of shape in columns.
		fscale    : None or numpy.ndarray
			If not None, fix the value of scale parameter (so not fitted)
		fshape    : None or numpy.ndarray
			If not None, fix the value of shape parameter (so not fitted)
		"""
		
		self._size = Y.size
		
		if self.n_bootstrap > 0:
			self.coefs_bootstrap = []
			if scale_cov is not None and scale_cov.ndim == 1:
				scale_cov = scale_cov.reshape( (scale_cov.size,1) )
			if shape_cov is not None and shape_cov.ndim == 1:
				shape_cov = shape_cov.reshape( (shape_cov.size,1) )
			if np.isscalar(fscale):
				fscale = np.array([fscale]).ravel()
			if np.isscalar(fshape):
				fshape = np.array([fshape]).ravel()
			
			loc    = np.array( [loc] ).ravel()
			if loc.size == 1: loc = np.repeat( loc , self._size )
			
			for i in range(self.n_bootstrap):
				idx = np.random.choice( self._size , self._size )
				Y_bs         = Y[idx]
				loc_bs       = loc[idx] 
				scale_cov_bs = None   if scale_cov is None else scale_cov[idx,:]
				shape_cov_bs = None   if shape_cov is None else shape_cov[idx,:]
				fscale_bs    = fscale if fscale    is None or fscale.size == 1 else fscale[idx]
				fshape_bs    = fshape if fshape    is None or fshape.size == 1 else fshape[idx]
				
				self._fit( Y_bs , loc_bs , scale_cov_bs , shape_cov_bs , fscale_bs , fshape_bs )
				self.coefs_bootstrap.append( self.coef_ )
			
			self.coefs_bootstrap = np.array( self.coefs_bootstrap )
			self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
		
		self._fit( Y , loc , scale_cov , shape_cov , fscale , fshape )
		
	##}}}
	
	def bootstrap_law( self , i ):##{{{
		"""
		Return a GPDLaw with coef from bootstrap
		
		Arguments
		---------
		i : integer
			Number of bootstrap
		
		Return
		------
		law : SDFC.GPDLaw
			A GPDLaw, None if n_bootstrap = 0
		"""
		if n_bootstrap == 0:
			return None
		law = GPDLaw( self.method , alpha = self.alpha )
		law._scale = self._scale.copy()
		law._shape = self._shape.copy()
		scale,shape = self._split_param( self.coefs_bootstrap[i,:] )
		law._scale.set_coef( scale )
		law._shape.set_coef( shape )
		law.coef_ = law._concat_param()
		law._update_param( law.coef_ )
		return law
	##}}}
	
	def predict_scale( self , scale_cov = None ):##{{{
		"""
		Return scale parameter with a new co-variates
		
		Arguments
		---------
		scale_cov : np.array or None
			Covariate
		
		Return
		------
		scale : np.array
			Location parameters
		"""
		return self._predict_param( self._scale , scale_cov )
	##}}}
	
	def predict_shape( self , shape_cov = None ):##{{{
		"""
		Return shape parameter with a new co-variates
		
		Arguments
		---------
		shape_cov : np.array or None
			Covariate
		
		Return
		------
		shape : np.array
			Location parameters
		"""
		return self._predict_param( self._shape , shape_cov )
	##}}}
	
	def _fit( self , Y , loc , scale_cov , shape_cov , fscale , fshape ):##{{{
		
		self._Y    = np.ravel(Y)
		self._size = Y.size
		
		self.loc    = np.array( [loc] ).ravel()
		if self.loc.size == 1: self.loc = np.repeat( loc , self._size )
		self._scale.init( X = scale_cov , fix_values = fscale , size = self._size )
		self._shape.init( X = shape_cov , fix_values = fshape , size = self._size )
		self._lparams = [self._scale,self._shape]
		
		if self.method == "moments":
			self._fit_moments()
		elif self.method == "lmoments":
			self._fit_lmoments()
		else:
			self._fit_mle()
		self.coef_ = self._concat_param()
		
	##}}}
	
	def _fit_moments( self ):##{{{
		
		## Center data
		idx_excess = (self._Y > self.loc)
		excess = self._Y[idx_excess] - self.loc[idx_excess]
		
		if self._scale.not_fixed():
			sX = self._scale.design_wo1()
			self._scale.set_coef( std( excess , sX , m = 0 , return_coef = True , linkFct = self._scale.linkFct ) )
		self._shape.set_intercept( self._shape.linkFct.inverse( -1e-8 ) )
		
		self._scale.update()
		self._shape.update()
		self.scale = self._scale.valueLf()
		self.shape = self._shape.valueLf()
	##}}}
	
	def _fit_lmoments( self ): ##{{{
		
		## Center data
		idx_excess = (self._Y > self.loc)
		excess = self._Y[idx_excess] - self.loc[idx_excess]
		
		## L-moments
		if self._scale.not_fixed() and self._shape.not_fixed():
			lmo1   = lmoments( excess , 1 )
			lmo2   = lmoments( excess , 2 )
			itau     = lmo1 / lmo2
			scale_lm = lmo1 * ( itau - 1 )
			scale_lm = scale_lm if scale_lm > 0 else 1e-8
			shape_lm = 2 - itau
			self._scale.set_intercept( self._scale.linkFct.inverse( scale_lm ) )
			self._shape.set_intercept( self._shape.linkFct.inverse( shape_lm ) )
		elif self._scale.not_fixed():
			self.shape = self._shape.valueLf()
			lmo1   = lmoments( excess , 1 )
			self.scale = lmo1 * ( 1 - self.shape )
			self.scale[ np.logical_not(self.scale > 0) ] = 1e-8
			self._scale.set_coef( mean( self.scale , self._scale.design_wo1() , linkFct = self._scale.linkFct , return_coef = True ) )
		elif self._shape.not_fixed():
			self.scale = self._scale.valueLf()
			excess /= self.scale[idx_excess]
			lmo1   = lmoments( excess , 1 )
			lmo2   = lmoments( excess , 2 )
			itau     = lmo1 / lmo2
			self._shape.set_intercept( self._shape.linkFct.inverse( 2 - itau ) )
		
		
		self._scale.update()
		self._shape.update()
		self.scale = self._scale.valueLf()
		self.shape = self._shape.valueLf()
		
	##}}}
	
	def _fit_mle( self ):##{{{
		
		## Initial condition
		self._fit_moments()
		param_mm = self._concat_param()
		lle_mm = self._negloglikelihood()
		
		self._fit_lmoments()
		param_lm = self._concat_param()
		lle_lm = self._negloglikelihood()
		
		## Keep best
		if np.isfinite(lle_mm) or np.isfinite(lle_lm):
			if lle_mm < lle_lm:
				self._update_param(param_mm)
			else:
				self._update_param(param_lm)
		else:
			self._scale.set_intercept( self._scale.linkFct.inverse(1.)   )
			self._shape.set_intercept( self._shape.linkFct.inverse(1e-2) )
		
		
		param_init = self._concat_param()
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = "BFGS" )
		self._update_param( self.optim_result.x )
	##}}}
	
#	def _split_param( self , param ):##{{{
#		param_scale = None
#		param_shape = None
#		
#		if self._scale.not_fixed():
#			param_scale = param[:self._scale.size]
#			if self._shape.not_fixed():
#				param_shape = param[self._scale.size:]
#		elif self._shape.not_fixed():
#			param_shape = param[:self._shape.size]
#		
#		return param_scale,param_shape
#	##}}}
#	
#	def _concat_param( self ):##{{{
#		return self._gen_concat_param( [self._scale,self._shape] )
#	##}}}
	
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
	
	def _update_param( self , param ):##{{{
		
		param_scale,param_shape = self._split_param(param)
		
		## Extract coefficients from param
		self._scale.set_coef( param_scale )
		self._shape.set_coef( param_shape )
		self._scale.update()
		self._shape.update()
		
		## Set scale and shape
		self.scale = np.ravel( self._scale.valueLf() )
		self.shape = np.ravel( self._shape.valueLf() )
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
		Z        = ( Y - loc ) / scale
		ZZ       = 1. + shape * Z
		exponent = 1. + 1. / shape
		
		grad = np.array([])
		
		if self._scale.not_fixed():
			gr_scale   = self._scale.valueGrLf()[idx_excess]
			grad_scale = self._scale.design_[idx_excess,:].T @ ( gr_scale * ( - exponent * shape * Z / ZZ / scale + 1. / scale ) )
			grad       = np.hstack( (grad,grad_scale) )
		
		if self._shape.not_fixed():
			gr_shape   = self._shape.valueGrLf()[idx_excess]
			grad_shape = self._shape.design_[idx_excess,:].T @ ( gr_shape * ( - np.log(ZZ) / shape**2 + exponent * Z / ZZ ) ) if np.all( ZZ > 0 ) else np.repeat(np.nan,self._shape.size)
			grad       = np.hstack( (grad,grad_shape) )
		
		return grad
	##}}}




