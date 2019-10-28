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
import scipy.special  as scs
import scipy.optimize as sco

from SDFC.__AbstractLaw            import AbstractLaw
from SDFC.NonParametric.__mean     import mean
from SDFC.NonParametric.__quantile import quantile
from SDFC.NonParametric.__lmoments import lmoments


#############
## Classes ##
#############

class GEV(AbstractLaw):
	"""
	SDFC.GEV
	========
	
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
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of Normal law
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , ["loc","scale","shape"] , method , n_bootstrap , alpha )
	##}}}
	
	def __str__(self):##{{{
		return self._to_str()
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	
	@property
	def loc(self):##{{{
		return self.params._dparams["loc"].value
	##}}}
	
	@property
	def scale(self):##{{{
		return self.params._dparams["scale"].value
	##}}}
	
	@property
	def shape(self):##{{{
		return self.params._dparams["shape"].value
	##}}}
	
	def predict_loc( self , c_loc = None ):##{{{
		"""
		Return location parameter with a new co-variates
		
		Arguments
		---------
		c_loc : np.array or None
			Covariate
		
		Return
		------
		loc : np.array
			Location parameters, if c_loc is None return self.loc
		"""
		return self._predict_covariate( "loc" , c_loc )
	##}}}
	
	def predict_scale( self , c_scale  = None ):##{{{
		"""
		Return scale parameter with a new co-variates
		
		Arguments
		---------
		c_scale : np.array or None
			Covariate
		
		Return
		------
		scale : np.array
			Scale parameters, if c_scale is None return self.scale
		"""
		return self._predict_covariate( "scale" , c_scale )
	##}}}
	
	def predict_shape( self , c_shape  = None ):##{{{
		"""
		Return scale parameter with a new co-variates
		
		Arguments
		---------
		c_shape : np.array or None
			Covariate
		
		Return
		------
		shape : np.array
			Shape parameters, if c_scale is None return self.shape
		"""
		return self._predict_covariate( "shape" , c_shape )
	##}}}
	
	
	def _fit_moments(self):##{{{
		ploc   = self.params._dparams["loc"]
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		m = np.mean(self._Y)
		s = np.sqrt(6) * np.std(self._Y) / np.pi
		
		iloc   = m - 0.57722 * s
		iscale = np.log(s)
		ishape = 1e-8
		
		
		## Fit scale
		if not pscale.is_fix():
			self.params.set_intercept( iscale , "scale" )
		
		## Fit loc
		if not ploc.is_fix():
			if pscale.is_fix():
				iloc = m - 0.57722 * np.exp(pscale.value)
				self.params.update_coef( mean( iloc , ploc.design_wo1() , value = False , link = ploc.link ) , "loc" )
			else:
				self.params.set_intercept( iloc , "loc" )
		
		## Fit shape
		if not pshape.is_fix():
			self.params.set_intercept( ishape , "shape" )
		
	##}}}
	
	def _fit_lmoments( self ): ##{{{
		
		ploc   = self.params._dparams["loc"]
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		lmom = lmoments( self._Y )
		
		tau3  = lmom[2] / lmom[1]
		co    = 2. / ( 3. + tau3 ) - np.log(2) / np.log(3)
		kappa = 7.8590 * co + 2.9554 * co**2
		g     = scs.gamma( 1. + kappa )
		
		
		iscale = lmom[1] * kappa / ( (1 - np.power( 2 , - kappa )) * g )
		iloc   = lmom[0] - iscale * (1 - g) / kappa
		ishape = - kappa
		
		## Fit scale
		if not pscale.is_fix():
			self.params.set_intercept( iscale , "scale" )
		
		## Fit loc
		if not ploc.is_fix():
			if pscale.is_fix():
				iloc = lmom[0] - pscale.value * (1 - g) / kappa
				self.params.update_coef( mean( iloc , ploc.design_wo1() , value = False , link = ploc.link ) , "loc" )
			else:
				self.params.set_intercept( iloc , "loc" )
		
		## Fit shape
		if not pshape.is_fix():
			self.params.set_intercept( ishape , "shape" )
	##}}}
	
	def _fit_lmoments_experimental(self):##{{{
		
		ploc   = self.params._dparams["loc"]
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		## First step, find lmoments
		c_Y = self.params.merge_covariate()
		if c_Y is None:
			self._fit_lmoments()
			return
		lmom = lmoments( self._Y , c_Y )
		
		## Find shape
		def uni_shape_solver(tau):
			bl,bu=-1,1
			fct = lambda x : 3 / 2 + tau / 2 - ( 1 - 3**x ) / (1 - 2**x )
			while fct(bl) * fct(bu) > 0:
				bl *= 2
				bu *= 2
			opt = sco.root_scalar( fct , method = "brenth" , bracket = [bl , bu] )
			return opt.root
		shape_solver = np.vectorize(uni_shape_solver)
		tau3 = lmom[:,2] / lmom[:,1]
		shape = shape_solver(tau3)
		
		## Find scale
		gshape = scs.gamma( 1 - shape )
		scale = - lmom[:,1] * shape / ( gshape * ( 1 - 2**shape ) )
		
		## Find loc
		loc = lmom[:,0] - scale * ( gshape - 1 ) / shape
		
		
		if not ploc.is_fix():
			self.params.update_coef( mean( loc   , ploc.design_wo1()   , link = ploc.link   , value = False ) , "loc"   )
		if not pscale.is_fix():
			self.params.update_coef( mean( scale , pscale.design_wo1() , link = pscale.link , value = False ) , "scale" )
		if not pshape.is_fix():
			self.params.update_coef( mean( shape , pshape.design_wo1() , link = pshape.link , value = False ) , "shape" )
	##}}}
	
	def _fit_quantiles( self ):##{{{
		
		ploc   = self.params._dparams["loc"]
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		## Fit loc
		if not ploc.is_fix():
			loc = quantile( self._Y , [np.exp(-1)] , c_Y = ploc.design_wo1() , value = True )
			self.params.update_coef( mean( loc , ploc.design_wo1() , link = ploc.link , value = False ) , "loc" )
		
		## Fit scale
		if not pscale.is_fix():
			qscale = np.array([0.25,0.5,0.75])
			coef   = -1. / np.log( - np.log(qscale) )
			qreg = quantile( self._Y - self.loc , qscale , pscale.design_wo1() )
			fscale = np.mean( (qreg * coef).reshape(-1,qscale.size) , axis = 1 ).reshape(-1,1)
			fscale[np.logical_not(fscale > 0)] = 0.1
			self.params.update_coef( mean( fscale , pscale.design_wo1() , link = pscale.link , value = False ) , "scale" )
		
		## Fit shape
		if not pshape.is_fix():
			p0,p1 = 0.1,0.9
			qval = quantile( (self._Y - self.loc) / self.scale , [p0,p1] , pshape.design_wo1() ).reshape(-1,2)
			kappa = qval[:,0] / qval[:,1]
			llp0,llp1 = np.log( - np.log( p0 ) ) , np.log( - np.log( p1 ) )
			shape = ( 2 * (llp0 - kappa * llp1 ) / ( llp0**2 - kappa * llp1**2 ) ).reshape(-1,1)
			self.params.update_coef( mean( shape , pshape.design_wo1() , link = pshape.link , value = False ) , "shape" )
	##}}}
	
	def _initialization_mle(self):##{{{
		
		self._fit_lmoments_experimental()
		
		nlll = self._negloglikelihood(self.coef_)
		grad = self._gradient_nlll(self.coef_)
		
		f_scale = 1
		f_shape = 1
		while ( not nlll < np.inf ) or np.any(np.isnan(grad)):
			pscale = self.params._dparams["scale"]
			pshape = self.params._dparams["shape"]
			
			if pshape.is_fix() and not pscale.is_fix():
				coef_ = np.zeros(pscale.n_features)
				coef_[0] = 1. * f_scale
				self.params.update_coef( coef_ , "scale" )
			elif not pshape.is_fix():
				coef_ = np.zeros(pshape.n_features)
				coef_[0] = 1e-1 / f_shape
				self.params.update_coef( coef_ , "shape" )
			else:
				self._fit_quantiles()
			f_scale *= 2
			f_shape *= 2
			nlll = self._negloglikelihood(self.coef_)
			grad = self._gradient_nlll(self.coef_)
	##}}}
	
	def _fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self._fit_moments()
		elif self.method == "lmoments":
			self._fit_lmoments()
		elif self.method == "lmoments-experimental":
			self._fit_lmoments_experimental()
		elif self.method == "quantiles":
			self._fit_quantiles()
	##}}}
	
	
	def _logZafun( self, Z , alpha ):##{{{
		return alpha * np.log( 1. + self.shape * Z )
	##}}}
	
	def _Zafun( self , Z , alpha ):##{{{
		return np.exp( self._logZafun( Z , alpha ) )
	##}}}
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		## Impossible scale
		if not np.all( self.scale > 0 ):
			return np.inf
		
		## Fuck exponential case
		zero_shape = ( np.abs(self.shape) < 1e-10 )
		shape = self.shape
		if np.any(zero_shape):
			shape[zero_shape] = 1e-10
		
		##
		Z = 1 + shape * ( self._Y - self.loc ) / self.scale
		
		if not np.all(Z > 0):
			return np.inf
		
		res = np.sum( ( 1. + 1. / shape ) * np.log(Z) + np.power( Z , - 1. / shape ) + np.log(self.scale) )
		
		
		return res if np.isfinite(res) else np.inf
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		

		zero_shape = ( np.abs(self.shape) < 1e-10 )
		shape = self.shape
		if np.any(zero_shape):
			shape[zero_shape] = 1e-10
		
		## Impossible
		if not np.all(self.scale > 0) or not np.all( 1. + shape * ( self._Y - self.loc ) / self.scale > 0 ):
			return np.zeros( coef.size ) + np.nan
		
		## Usefull values
		Z      = ( self._Y - self.loc ) / self.scale
		Za1    = self._Zafun( Z , 1. )
		ishape = 1. / shape
		Zamsi  = self._Zafun( Z , - ishape ) ## Za of Minus Shape Inverse
		
		## Gradient
		grad = np.array( [] )
		
		ploc = self.params._dparams["loc"]
		if not ploc.is_fix():
			loc_vect   = ploc.gradient()   * ( Zamsi - 1 - shape ) / ( self.scale * Za1 )
			grad_loc   = np.dot( ploc.design_.T   , loc_vect   )
			grad = np.hstack( (grad,grad_loc.squeeze()) )
		
		pscale = self.params._dparams["scale"]
		if not pscale.is_fix():
			scale_vect = pscale.gradient() * ( 1. + Z * ( Zamsi - 1 - shape ) / Za1 ) / self.scale
			grad_scale = np.dot( pscale.design_.T , scale_vect )
			grad = np.hstack( (grad,grad_scale.squeeze()) )
		
		pshape = self.params._dparams["shape"]
		if not pshape.is_fix():
			shape_vect = pshape.gradient() * ( ( Zamsi - 1. ) * np.log(Za1) * ishape**2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1 )
			grad_shape = np.dot( pshape.design_.T , shape_vect )
			grad = np.hstack( (grad,grad_shape.squeeze()) )
		return grad



	##}}}
