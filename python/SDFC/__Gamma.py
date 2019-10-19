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
import scipy.special  as scp

from SDFC.__AbstractLaw        import AbstractLaw
from SDFC.tools.__LawParam     import LawParam
from SDFC.tools.__LinkFct      import IdLinkFct
from SDFC.NonParametric.__mean import mean
from SDFC.NonParametric.__var  import var


#############
## Classes ##
#############

class Gamma(AbstractLaw):
	"""
	SDFC.Gamma
	==========
	
	Fit parameters of a Gamma law, possibly with co-variable
	
	Attributes
	----------
	method : string
		method used to fit
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
		Initialization of Gamma law
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		"""
		AbstractLaw.__init__( self , ["scale","shape"] , method , n_bootstrap , alpha )
	##}}}
	
	def __str__(self):##{{{
		return self._to_str()
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
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
		n_samples = pscale.n_samples
		
		if not pscale.is_fix() and not pshape.is_fix():
			mX = np.ones( (n_samples,1) )
			vX = np.ones( (n_samples,1) )
			for i in range(1,pscale.n_features):
				for j in range(pshape.n_features):
					mX = np.hstack( (mX,np.reshape( pscale.design_[:,i]    * pshape.design_[:,j] , (n_samples,1) ) ) )
					vX = np.hstack( (vX,np.reshape( pscale.design_[:,i]**2 * pshape.design_[:,j] , (n_samples,1) ) ) )
			m = mean( self._Y , mX )
			v = var(  self._Y , vX )
			
			idx  = np.logical_or( np.abs(m) < 1e-8 , v < 1e-8 )
			cidx = np.logical_not(idx)
			scale = np.zeros_like(m)
			shape = np.zeros_like(m)
			scale[cidx] = v[cidx] / m[cidx]
			shape[cidx] = m[cidx]**2 / v[cidx]
			
			if np.any(idx):
				scale[idx] = scale[cidx].min()
				shape[idx] = shape[cidx].min()
			
			self.params.update_coef( mean( scale , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
			self.params.update_coef( mean( shape , pshape.design_wo1() , value = False , link = pshape.link ) , "shape" )
		elif pscale.is_fix():
			
			m = mean( self._Y  , pshape.design_wo1() * self.scale )
			v = var(  self._Y  , pshape.design_wo1() * self.scale**2 )
			
			shape = m**2 / v
			self.params.update_coef( mean( shape , pshape.design_wo1() , value = False , link = pshape.link ) , "shape" )
			
		elif pshape.is_fix():
			m = mean( self._Y  , pscale.design_wo1()    * self.shape )
			v = var(  self._Y  , pscale.design_wo1()**2 * self.shape )
			
			scale = v / m
			self.params.update_coef( mean( scale , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
	##}}}
	
	def _fit_mle(self):##{{{
		self._fit_moments()
		AbstractLaw._fit_mle(self)
	##}}}
	
	def _fit( self ): ##{{{
		if self.method == "moments":
			self._fit_moments()
		else:
			self._fit_mle()
	##}}}
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		if not np.all(self.scale > 0) or not np.all(self.shape > 0) or not np.all(self._Y > 0):
			return np.Inf
		
		return np.sum( self._Y / self.scale + scp.loggamma(self.shape) + self.shape * np.log(self.scale) - (self.shape-1) * np.log(self._Y) )
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		grad = np.array( [] )
		
		pscale = self.params._dparams["scale"]
		pshape = self.params._dparams["shape"]
		
		if np.all(self.scale > 0) and np.all(self.shape > 0) and np.all(self._Y > 0):
			if not pscale.is_fix():
				grad_scale = pscale.design_.T @ ( ( self.shape / self.scale - self._Y / self.scale**2 ) * pscale.gradient() )
				grad = np.hstack( (grad,grad_scale.squeeze()) )
			if not pshape.is_fix():
				grad_shape = pshape.design_.T @ ( ( scp.digamma(self.shape) + np.log(self.scale) - np.log(self._Y) ) * pshape.gradient() )
				grad = np.hstack( (grad,grad_shape.squeeze()) )
		else:
			grad = np.zeros( coef.size ) + np.nan
		return grad
	##}}}





