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
import scipy.optimize as sco
import scipy.special  as scs

from SDFC.__AbstractLaw            import AbstractLaw
from SDFC.NonParametric.__mean     import mean
from SDFC.NonParametric.__quantile import quantile
from SDFC.NonParametric.__lmoments import lmoments
from SDFC.tools.__LawParam         import LawParam
from SDFC.tools.__LinkFct          import IdLinkFct



###########
## Class ##
###########

class GEVLaw(AbstractLaw):
	"""
	SDFC.GEVLaw
	===========
	
	Fit parameters of a Generalized Extreme Value law, possibly with co-variable
	
	Attributes
	----------
	method : string
		method used to fit
	loc    : numpy.ndarray
		Location fitted
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
	
	def __init__( self , method = "MLE" , link_fct_loc = IdLinkFct() , link_fct_scale = IdLinkFct() , link_fct_shape = IdLinkFct() , n_bootstrap = 0 , alpha = 0.05 ):##{{{
		"""
		Initialization of GEVLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments", "lmoments", "quantiles" and "MLE" (Maximum Likelihood estimation)
		link_fct_loc   : a class herited from SDFC.tools.LinkFct
			Link function for loc, default is IdLinkFct()
		link_fct_scale : a class herited from SDFC.tools.LinkFct
			Link function for scale, default is SDFC.tools.IdLinkFct(). Interesting option is SDFC.tools.ExpLinkFct().
		link_fct_shape : a class herited from SDFC.tools.LinkFct
			Link function for shape, default is IdLinkFct(). Interesting option is SDFC.tools.LogitLinkFct( -0.5 , 0.5 ).
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		
		"""
		AbstractLaw.__init__( self , method , n_bootstrap , alpha )
		
		self.loc    = None
		self.scale  = None
		self.shape  = None
		
		self._loc   = LawParam( linkFct = link_fct_loc   , kind = "loc"   )
		self._scale = LawParam( linkFct = link_fct_scale , kind = "scale" )
		self._shape = LawParam( linkFct = link_fct_shape , kind = "shape" )
		
		self._lparams = [self._loc,self._scale,self._shape]
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.GEVLaw\n"
		val += "-----------\n"
		val += "* fit_method : {}\n".format(self.method)
		val += "* link_loc   : {}\n".format(str(self._loc.linkFct))
		val += "* link_scale : {}\n".format(str(self._scale.linkFct))
		val += "* link_shape : {}\n".format(str(self._shape.linkFct))
		val += "* coef_      : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def fit( self , Y , loc_cov = None , scale_cov = None , shape_cov = None , floc = None , fscale = None , fshape = None ):##{{{
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
		floc      : None or numpy.ndarray (see Notes)
			If not None, fix the value of loc parameter (so not fitted)
		fscale    : None or numpy.ndarray (see Notes)
			If not None, fix the value of scale parameter (so not fitted)
		fshape    : None or numpy.ndarray (see Notes)
			If not None, fix the value of shape parameter (so not fitted)
		
		Notes: floc, fscale and fshape are NOT supported for method = "moments" and method = "lmoments"
		"""
		
		self._size = Y.size
		
		if self.n_bootstrap > 0:
			self.coefs_bootstrap = []
			if loc_cov is not None and loc_cov.ndim == 1:
				loc_cov = loc_cov.reshape( (loc_cov.size,1) )
			if scale_cov is not None and scale_cov.ndim == 1:
				scale_cov = scale_cov.reshape( (scale_cov.size,1) )
			if shape_cov is not None and shape_cov.ndim == 1:
				shape_cov = shape_cov.reshape( (shape_cov.size,1) )
			if np.isscalar(floc):
				floc = np.array([floc]).ravel()
			if np.isscalar(fscale):
				fscale = np.array([fscale]).ravel()
			if np.isscalar(fshape):
				fshape = np.array([fshape]).ravel()
			
			
			for i in range(self.n_bootstrap):
				idx = np.random.choice( self._size , self._size )
				Y_bs         = Y[idx]
				loc_cov_bs   = None   if loc_cov   is None else loc_cov[idx,:]
				scale_cov_bs = None   if scale_cov is None else scale_cov[idx,:]
				shape_cov_bs = None   if shape_cov is None else shape_cov[idx,:]
				floc_bs      = floc   if floc      is None or floc.size == 1   else floc[idx]
				fscale_bs    = fscale if fscale    is None or fscale.size == 1 else fscale[idx]
				fshape_bs    = fshape if fshape    is None or fshape.size == 1 else fshape[idx]
				
				self._fit( Y_bs , loc_cov_bs , scale_cov_bs , shape_cov_bs , floc_bs , fscale_bs , fshape_bs )
				self.coefs_bootstrap.append( self.coef_ )
			
			self.coefs_bootstrap = np.array( self.coefs_bootstrap )
			self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
		
		self._fit( Y , loc_cov , scale_cov , shape_cov , floc , fscale , fshape )
		
	##}}}
	
	def predict_loc( self , loc_cov = None ):##{{{
		return self._predict_param( self._loc , loc_cov )
	##}}}
	
	def predict_scale( self , scale_cov = None ):##{{{
		return self._predict_param( self._scale , scale_cov )
	##}}}
	
	def predict_shape( self , shape_cov = None ):##{{{
		return self._predict_param( self._shape , shape_cov )
	##}}}
	
	def _fit( self , Y , loc_cov = None , scale_cov = None , shape_cov = None , floc = None , fscale = None , fshape = None ):##{{{
		self._Y    = np.ravel(Y)
		self._size = Y.size
		
		self._loc.init(   X = loc_cov   , fix_values = floc   , size = self._size )
		self._scale.init( X = scale_cov , fix_values = fscale , size = self._size )
		self._shape.init( X = shape_cov , fix_values = fshape , size = self._size )
		
		
		if self.method == "moments":
			self._fit_moments()
		elif self.method == "lmoments":
			self._fit_lmoments()
		elif self.method == "quantiles":
			self._fit_quantiles()
		else:
			self._fit_mle()
		
		self.coef_ = self._concat_param()
	##}}}
	
	def _fit_moments( self ):##{{{
		
		m = np.mean(self._Y)
		s = np.sqrt(6) * np.std(self._Y) / np.pi
		
		iloc   = m - 0.57722 * s
		iscale = np.log(s)
		ishape = 1e-8
		
		self._loc.set_intercept(   self._loc.linkFct.inverse( iloc )     )
		self._scale.set_intercept( self._scale.linkFct.inverse( iscale ) )
		self._shape.set_intercept( self._shape.linkFct.inverse( ishape ) )
		
		
		self._loc.update()
		self._scale.update()
		self._shape.update()
		self.loc   = self._loc.valueLf()
		self.scale = self._scale.valueLf()
		self.shape = self._shape.valueLf()
		
	##}}}
	
	def _fit_lmoments( self ): ##{{{
		lmom1 = lmoments( self._Y , 1 )
		lmom2 = lmoments( self._Y , 2 )
		lmom3 = lmoments( self._Y , 3 )
		
		tau3  = lmom3 / lmom2
		co    = 2. / ( 3. + tau3 ) - np.log(2) / np.log(3)
		kappa = 7.8590 * co + 2.9554 * co**2
		g     = scs.gamma( 1. + kappa )
		
		iscale = lmom2 * kappa / ( (1 - np.power( 2 , - kappa )) * g )
		iloc   = lmom1 - iscale * (1 - g) / kappa
		ishape = - kappa
		
		self._loc.set_intercept(   self._loc.linkFct.inverse( iloc )     )
		self._scale.set_intercept( self._scale.linkFct.inverse( iscale ) )
		self._shape.set_intercept( self._shape.linkFct.inverse( ishape ) )
		
		self._loc.update()
		self._scale.update()
		self._shape.update()
		self.loc   = self._loc.valueLf()
		self.scale = self._scale.valueLf()
		self.shape = self._shape.valueLf()
		
	##}}}
	
	def _fit_quantiles( self ):##{{{
		## Fit loc
		if self._loc.not_fixed():
			if self._loc.size == 1:
				self._loc.coef_[0] = self._loc.linkFct.inverse( quantile( self._Y , [np.exp(-1)] , return_coef = True ) )
			else:
				loc = quantile( self._Y , [np.exp(-1)] , self._loc.design_wo1() )
				self._loc.coef_ = mean( loc , self._loc.design_wo1() , linkFct = self._loc.linkFct , return_coef = True ).ravel()
		self._loc.update()
		self.loc = self._loc.valueLf()
		
		
		## Fit scale
		if self._scale.not_fixed():
			qscale = np.array([0.25,0.5,0.75])
			coef   = -1. / np.log( - np.log(qscale) )
			if self._scale.size == 1:
				self._scale.coef_[0] = self._scale.linkFct.inverse( np.mean( np.percentile( self._Y - self.loc , 100 * qscale ) * coef ) )
			else:
				qreg = quantile( self._Y - self.loc , qscale , self._scale.design_wo1() )
				fscale = np.mean( qreg * coef , axis = 1 )
				fscale[np.logical_not(fscale > 0)] = 0.1
				self._scale.set_coef( mean( fscale , self._scale.design_wo1() , linkFct = self._scale.linkFct , return_coef = True ) )
		self._scale.update()
		self.scale = self._scale.valueLf()
		
		## Fit shape
		if self._shape.not_fixed():
			p0,p1 = 0.1,0.9
			if self._shape.size == 1:
				qval = np.percentile( (self._Y - self.loc) / self.scale , [100*p0,100*p1] )
				kappa = qval[0] / qval[1]
				llp0,llp1 = np.log( - np.log( p0 ) ) , np.log( - np.log( p1 ) )
				self._shape.coef_[0] = self._shape.linkFct.inverse( 2 * (llp0 - kappa * llp1 ) / ( llp0**2 - kappa * llp1**2 ) )
			else:
				qval = quantile( (self._Y - self.loc) / self.scale , [p0,p1] , self._shape.design_wo1() )
				kappa = qval[:,0] / qval[:,1]
				llp0,llp1 = np.log( - np.log( p0 ) ) , np.log( - np.log( p1 ) )
				shape = 2 * (llp0 - kappa * llp1 ) / ( llp0**2 - kappa * llp1**2 )
				
				self._shape.set_coef( mean( shape , self._shape.design_wo1() , linkFct = self._shape.linkFct , return_coef = True ) )
		self._shape.update()
		self.shape = self._shape.valueLf()
		
	##}}}
	
	def _fit_mle( self ):##{{{
		
		## Initial condition
#		self._fit_moments()
#		init_loc_mm   = self._loc.coef_.copy()
#		init_scale_mm = self._scale.coef_.copy()
#		init_shape_mm = self._shape.coef_.copy()
#		lle_mm = self._negloglikelihood()
#		
#		self._fit_lmoments()
#		init_loc_lm   = self._loc.coef_.copy()
#		init_scale_lm = self._scale.coef_.copy()
#		init_shape_lm = self._shape.coef_.copy()
#		lle_lm = self._negloglikelihood()
		
		self._fit_quantiles()
#		init_loc_q   = self._loc.coef_.copy()
#		init_scale_q = self._scale.coef_.copy()
#		init_shape_q = self._shape.coef_.copy()
#		lle_q  = self._negloglikelihood()
		
		
		## Keep best
#		if np.isfinite(lle_mm) or np.isfinite(lle_lm) or np.isfinite(lle_q):
#			if lle_mm < lle_lm and lle_mm < lle_q:
#				self._loc.set_coef( init_loc_mm )
#				self._scale.set_coef( init_scale_mm )
#				self._shape.set_coef( init_shape_mm )
#			elif lle_lm < lle_mm and lle_lm < lle_q:
#				self._loc.set_coef( init_loc_lm )
#				self._scale.set_coef( init_scale_lm )
#				self._shape.set_coef( init_shape_lm )
#			else:
#				self._loc.set_coef( init_loc_q )
#				self._scale.set_coef( init_scale_q )
#				self._shape.set_coef( init_shape_q )
#		else:
#			self._loc.set_intercept(   self._loc.linkFct.inverse(0.) )
#			self._scale.set_intercept( self._scale.linkFct.inverse(1.) )
#			self._shape.set_intercept( self._shape.linkFct.inverse(1e-2) )
		
		
		## Optimization
		param_init = self._concat_param()
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = "BFGS" )
		self._update_param( self.optim_result.x )
	##}}}
	
	def _split_param( self , param ):##{{{
		param_loc   = None
		param_scale = None
		param_shape = None
		
		
		if self._loc.not_fixed() and self._scale.not_fixed() and self._shape.not_fixed():
			param_loc   = param[:self._loc.size]
			param_scale = param[self._loc.size:(self._loc.size+self._scale.size)]
			param_shape = param[(self._loc.size+self._scale.size):]
		elif self._loc.not_fixed() and self._scale.not_fixed():
			param_loc   = param[:self._loc.size]
			param_scale = param[self._loc.size:]
			param_shape = None
		elif self._scale.not_fixed() and self._shape.not_fixed():
			param_loc   = None
			param_scale = param[:self._scale.size]
			param_shape = param[self._scale.size:]
		elif self._loc.not_fixed() and self._shape.not_fixed():
			param_loc   = param[:self._loc.size]
			param_scale = None
			param_shape = param[self._loc.size:]
		elif self._loc.not_fixed():
			param_loc   = param
			param_scale = None
			param_shape = None
			pass
		elif self._scale.not_fixed():
			param_loc   = None
			param_scale = param
			param_shape = None
			pass
		elif self._shape.not_fixed():
			param_loc   = None
			param_scale = None
			param_shape = param
		
		return param_loc,param_scale,param_shape
	##}}}
	
	def _concat_param( self ):##{{{
		param = None
		param_loc   = self._loc.coef_   if self._loc.not_fixed()   else np.array([])
		param_scale = self._scale.coef_ if self._scale.not_fixed() else np.array([])
		param_shape = self._shape.coef_ if self._shape.not_fixed() else np.array([])
		
		param = np.hstack( (param_loc,param_scale,param_shape) )

#		if self._loc.not_fixed() and self._scale.not_fixed() and self._shape.not_fixed():
#			param = np.hstack( (self._loc.coef_,self._scale.coef_,self._shape.coef_) )
#		elif self._loc.not_fixed() and self._scale.not_fixed():
#			param = np.hstack( (self._loc.coef_,self._scale.coef_) )
#		elif self._loc.not_fixed() and self._shape.not_fixed():
#			param = np.hstack( (self._loc.coef_,self._shape.coef_) )
#		elif self._scale.not_fixed() and self._shape.not_fixed():
#			param = np.hstack( (self._scale.coef_,self._shape.coef_) )
#		elif self._loc.not_fixed():
#			param = self._loc.coef_
#		elif self._scale.not_fixed():
#			param = self._scale.coef_
#		elif self._shape.not_fixed():
#			param = self._shape.coef_
		return param
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
	
	def _update_param( self , param ):##{{{
		
		param_loc,param_scale,param_shape = self._split_param( param )
		
		## Extract coefficients from param
		self._loc.set_coef( param_loc )
		self._scale.set_coef( param_scale )
		self._shape.set_coef( param_shape )
		self._loc.update()
		self._scale.update()
		self._shape.update()
		
		## Set scale and shape
		self.loc   = self._loc.valueLf()
		self.scale = self._scale.valueLf()
		self.shape = self._shape.valueLf()
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
			return np.zeros( param.size ) + np.nan
		
		## Usefull values
		Z      = ( self._Y - self.loc ) / self.scale
		Za1    = self._Zafun( Z , 1. )
		ishape = 1. / self.shape
		Zamsi  = self._Zafun( Z , - ishape ) ## Za of Minus Shape Inverse
		
		## Gradient
		grad = np.array( [] )
		
		if self._loc.not_fixed():
			loc_vect   = self._loc.valueGrLf()   * ( Zamsi - 1 - self.shape ) / ( self.scale * Za1 )
			grad_loc   = np.dot( self._loc.design_.T   , loc_vect   )
			grad = np.hstack( (grad,grad_loc) )
		if self._scale.not_fixed():
			scale_vect = self._scale.valueGrLf() * ( 1. + Z * ( Zamsi - 1 - self.shape ) / Za1 ) / self.scale
			grad_scale = np.dot( self._scale.design_.T , scale_vect )
			grad = np.hstack( (grad,grad_scale) )
		if self._shape.not_fixed():
			shape_vect = self._shape.valueGrLf() * ( ( Zamsi - 1. ) * np.log(Za1) * ishape**2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1 )
			grad_shape = np.dot( self._shape.design_.T , shape_vect )
			grad = np.hstack( (grad,grad_shape) )
		
		return grad
		
	##}}}
	



