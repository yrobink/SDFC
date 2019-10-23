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

import numpy         as np
import scipy.special as scs

from SDFC.__AbstractLaw            import AbstractLaw
from SDFC.NonParametric.__mean     import mean
from SDFC.NonParametric.__std      import std
from SDFC.NonParametric.__lmoments import lmoments


#############
## Classes ##
#############

class GPD(AbstractLaw):
	"""
	SDFC.GPD
	========
	
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
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		idx = (self._Y > self.loc)
		Y   = (self._Y[idx] - self.loc[idx]).reshape(-1,1)
		if not pscale.is_fix():
			c_scale = pscale.design_wo1()
			if c_scale is not None:
				c_scale = c_scale.reshape(-1,pscale.n_features-1)[idx]
			self.params.update_coef( std( Y , c_scale , m_Y = 0 , value = False , link = pscale.link ) , "scale" )
		
		if not pshape.is_fix():
			self.params.set_intercept( -1e-8 , "shape" )
	##}}}
	
	def _fit_lmoments( self ): ##{{{
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		idx = (self._Y > self.loc)
		Y   = (self._Y[idx] - self.loc[idx]).reshape(-1,1)
		
		## L-moments
		lmom = lmoments( Y )
		if not pscale.is_fix() and not pshape.is_fix():
			itau     = lmom[0] / lmom[1]
			scale_lm = lmom[0] * ( itau - 1 )
			scale_lm = scale_lm if scale_lm > 0 else 1e-8
			shape_lm = 2 - itau
			self.params.set_intercept( scale_lm , "scale" )
			self.params.set_intercept( shape_lm , "shape" )
		elif not pscale.is_fix():
			scale = lmom[0] * ( 1 - self.shape )
			scale[ np.logical_not(scale > 0) ] = 1e-8
			self.params.update_coef( mean( scale , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
		elif not pshape.is_fix():
			Y /= self.scale[idx].reshape(-1,1)
			lmom = lmoments(Y)
			itau     = lmom[0] / lmom[1]
			self.params.set_intercept( 2 - itau , "shape" )
	##}}}
	
	def _fit_lmoments_experimental( self ): ##{{{
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		idx = (self._Y > self.loc)
		Y   = (self._Y[idx] - self.loc[idx]).reshape(-1,1)
		
		## First step, find lmoments
		c_Y = self.params.merge_covariate()
		if c_Y is None:
			self._fit_lmoments()
			return
		c_Y = c_Y[idx.squeeze(),:]
		lmom = lmoments( Y , c_Y )
		
		if not pscale.is_fix() and not pshape.is_fix():
			itau  = lmom[:,0] / lmom[:,1]
			scale = lmom[:,0] * ( itau - 1 )
			shape = 2 - itau
			self.params.update_coef( mean( scale , pscale.design_wo1() , link = pscale.link , value = False ) , "scale" )
			self.params.update_coef( mean( shape , pshape.design_wo1() , link = pshape.link , value = False ) , "shape" )
		elif not pscale.is_fix():
			scale = lmom[:,0] * ( 1 - self.shape )
			self.params.update_coef( mean( scale , pscale.design_wo1() , link = pscale.link , value = False ) , "scale" )
		elif not pshape.is_fix():
			Y    /= self.scale[idx].reshape(-1,1)
			lmom  = lmoments( Y , pshape.design_wo1() )
			shape = 2 - lmom[:,0] / lmom[:,1]
			self.params.update_coef( mean( shape , pshape.design_wo1() , link = pshape.link , value = False ) , "shape" )
	##}}}
	
	def _fit_mle_initialization(self):##{{{
		self._fit_lmoments_experimental()
#		self._fit_moments()
#		nlll_mom = self._negloglikelihood(self.coef_)
#		grad_mom = np.any(np.isnan(self._gradient_nlll(self.coef_)))
#		
#		self._fit_lmoments()
#		nlll_lmo = self._negloglikelihood(self.coef_)
#		grad_lmo = np.any(np.isnan(self._gradient_nlll(self.coef_)))
#		
#		if grad_mom and grad_lmo and ( nlll_mom < np.inf or nlll_lmo < np.inf ):
#			if nlll_mom < nlll_lmo:
#				self._fit_moments()
#			else:
#				self._fit_lmoments()
#		elif grad_mom and nlll_mom < np.inf:
#			self._fit_moments()
#		elif grad_lmo and nlll_lmo < np.inf:
#			self._fit_lmoments()
#		else:
#			pscale = self.params._dparams["scale"]
#			pshape = self.params._dparams["shape"]
#			if not pscale.is_fix():
#				self.params.set_intercept( 1.   , "scale" )
#			if not pshape.is_fix():
#				self.params.set_intercept( 0.01 , "shape" )
		
	##}}}
	
	def _fit_mle(self):##{{{
		self._fit_mle_initialization()
		AbstractLaw._fit_mle(self)
	##}}}
	
	def _fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self._fit_moments()
		elif self.method == "lmoments":
			self._fit_lmoments()
		elif self.method == "lmoments-experimental":
			self._fit_lmoments_experimental()
		else:
			self._fit_mle()
	##}}}
	
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		## Impossible scale
		if not np.all( self.scale > 0 ):
			return np.inf
		
		## Remove exponential case
		shape = self.shape
		zero_shape = ( np.abs(shape) < 1e-10 )
		if np.any(zero_shape):
			shape[zero_shape] = -1e-10
		
		
		##
		idx   = (self._Y > self.loc).squeeze()
		loc   = self.loc[idx,:]
		scale = self.scale[idx,:]
		shape = shape[idx,:]
		Z = 1. + shape * ( self._Y[idx,:] - loc ) / scale
		
		if not np.all(Z > 0):
			return np.inf
		res = np.sum( np.log( scale ) + np.log(Z) * ( 1 + 1. / shape ) )
		return res
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		
		## Remove exponential case
		shape = self.shape
		zero_shape = ( np.abs(shape) < 1e-10 )
		if np.any(zero_shape):
			shape[zero_shape] = -1e-10
		
		
		##
		idx   = (self._Y > self.loc).squeeze()
		Y      = self._Y[idx,:]
		loc   = self.loc[idx,:]
		scale = self.scale[idx,:]
		shape = shape[idx,:]
		
		Z        = ( Y - loc ) / scale
		ZZ       = 1. + shape * Z
		exponent = 1. + 1. / shape
		
		grad = np.array([])
		
		pscale = self.params._dparams["scale"]
		if not pscale.is_fix():
			gr_scale   = pscale.gradient()[idx,:]
			A = gr_scale * ( - exponent * shape * Z / ZZ / scale + 1. / scale )
			B = pscale.design_[idx,:].T
			grad_scale =  B @ A
			grad       = np.hstack( (grad,grad_scale.squeeze()) )
		
		pshape = self.params._dparams["shape"]
		if not pshape.is_fix():
			gr_shape   = pshape.gradient()[idx,:].reshape(-1,1)
			grad_shape = pshape.design_[idx,:].T @ ( gr_shape * ( - np.log(ZZ) / shape**2 + exponent * Z / ZZ ) ) if np.all( ZZ > 0 ) else np.repeat(np.nan,pshape.n_features)
			grad       = np.hstack( (grad,grad_shape.squeeze()) )
		return grad
	##}}}

