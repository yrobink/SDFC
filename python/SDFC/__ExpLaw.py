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

from SDFC.__AbstractLaw        import AbstractLaw
from SDFC.tools.__LawParam     import LawParam
from SDFC.tools.__LinkFct      import IdLinkFct
from SDFC.NonParametric.__mean import mean
from SDFC.NonParametric.__std  import std


#############
## Classes ##
#############

class ExpLaw(AbstractLaw):
	"""
	SDFC.ExpLaw
	===========
	
	Fit parameters of an Exponential law, possibly with co-variable
	
	Attributes
	----------
	method : string
		method used to fit
	scale  : numpy.ndarray
		Scale fitted
	coef_  : numpy.ndarray
		Coefficients fitted
	coefs_bootstrap: numpy.ndarray
		coef_ for each bootstrap
	confidence interval: numpy.ndarray[ shape = (2,coef_.size) ]
		Confidence interval, first line is the alpha/2 quantile, and second line the 1 - alpha/2 quantile
	alpha          : float
		Level of confidence interval
	"""
	
	def __init__( self , method = "MLE" , link_fct_scale = IdLinkFct() , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of ExpLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments", "lmoments", "quantiles" and "MLE" (Maximum Likelihood estimation)
		link_fct_scale : a class herited from SDFC.tools.LinkFct
			Link function for scale, default is SDFC.tools.IdLinkFct(). Interesting option is SDFC.tools.ExpLinkFct().
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , method , n_bootstrap , alpha )
		
		self.scale     = None
		
		self._scale = LawParam( linkFct = link_fct_scale , kind = "scale" )
		
		self._lparams = [self._scale]
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.ExpLaw\n"
		val += "-----------\n"
		val += "* fit_method : {}\n".format(self.method)
		val += "* link_scale : {}\n".format(str(self._scale.linkFct))
		val += "* coef_      : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def fit( self , Y , scale_cov = None ): ##{{{
		"""
		Fit function for ExpLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		"""
		self._size = Y.size
		
		
		if self.n_bootstrap > 0:
			self.coefs_bootstrap = []
			if scale_cov is not None and scale_cov.ndim == 1:
				scale_cov = scale_cov.reshape( (scale_cov.size,1) )
			
			
			for i in range(self.n_bootstrap):
				idx = np.random.choice( self._size , self._size )
				Y_bs         = Y[idx]
				scale_cov_bs = None   if scale_cov is None else scale_cov[idx,:]
				
				self._fit( Y_bs , scale_cov_bs )
				self.coefs_bootstrap.append( self.coef_ )
			
			self.coefs_bootstrap = np.array( self.coefs_bootstrap )
			self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
		
		self._fit( Y , scale_cov )
	##}}}
	
	def bootstrap_law( self , i ):##{{{
		"""
		Return an ExpLaw with coef from bootstrap
		
		Arguments
		---------
		i : integer
			Number of bootstrap
		
		Return
		------
		law : SDFC.ExpLaw
			A ExpLaw, None if n_bootstrap = 0
		"""
		if self.n_bootstrap == 0:
			return None
		law = ExpLaw( self.method , alpha = self.alpha )
		law._scale = self._scale.copy()
		law._scale.set_coef( self._scale.coef_ )
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
	
	def _fit( self , Y , scale_cov = None ):##{{{
		self._Y    = np.ravel(Y)
		self._scale.init( X = scale_cov , fix_values = None , size = self._size )
		
		if self.method == "moments":
			self._fit_moments()
		else:
			self._fit_mle()
		
		self.coef_ = self._scale.coef_
		##}}}
	
	def _fit_moments(self):##{{{
		
		sX = self._scale.design_wo1()
		self._scale.set_coef( mean( self._Y , sX , return_coef = True , linkFct = self._scale.linkFct ) )
		self._scale.update()
		self.scale = self._scale.valueLf()
	##}}}
	
	def _fit_mle(self):##{{{
		
		self._fit_moments()
		
		param_init = self._scale.coef_
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = "BFGS" )
		self._update_param(self.optim_result.x)
	##}}}
	
	def _negloglikelihood( self ): ##{{{
		if not np.all(self.scale > 0):
			return np.Inf
		
		return np.sum( np.log(self.scale) - self._Y / self.scale )
	##}}}
	
	def _update_param( self , param ):##{{{
		
		self._scale.set_coef( param )
		self._scale.update()
		
		self.scale = np.ravel( self._scale.valueLf() )
	##}}}
	
	def _optim_function( self , param ):##{{{
		self._update_param(param)
		return self._negloglikelihood()
	##}}}
	
	def _gradient_optim_function( self , param ): ##{{{
		self._update_param(param)
		grad_scale = np.zeros_like(param) + np.nan
		
		if np.all(self.scale > 0):
			grad_scale = self._scale.design_.T @ ( ( 1. / self.scale - self._Y / self.scale**2 ) * self._scale.valueGrLf() )
		
		return grad_scale
	##}}}

