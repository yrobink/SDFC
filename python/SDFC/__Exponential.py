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

from SDFC.__AbstractLaw        import AbstractLaw
from SDFC.NonParametric.__mean import mean


#############
## Classes ##
#############

class Exponential(AbstractLaw):
	"""
	Class to fit an Exponential law with covariates, available methods are:
	
	moments  : use empirical estimator of mean and standard deviation to find loc and scale, possibly with least square
			   regression if covariates are given
	bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of n_mcmc_iteration sample draw from
	           the posterior P(coef_ | Y)
	mle      : Maximum likelihood estimation
	
	Parameters
	==========
	scale : scale parameter
	"""
	__doc__ += AbstractLaw.__doc__
	
	def __init__( self , method = "MLE" , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of Exponential law
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments", "lmoments", "quantiles" and "MLE" (Maximum Likelihood estimation)
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , ["scale"] , method , n_bootstrap , alpha )
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
	
	
	def _fit_moments(self):##{{{
		
		pscale = self.params._dparams["scale"]
		
		## Fit loc
		if not pscale.is_fix():
			self.params.update_coef( mean( self._Y , pscale.design_wo1() , value = False , link = pscale.link ) , "scale" )
	##}}}
	
	def _initialization_mle(self):##{{{
		self._fit_moments()
	##}}}
	
	def _fit( self ):##{{{
		
		## Fit itself
		if self.method == "moments":
			self._fit_moments()
	##}}}
	
	@AbstractLaw._update_coef
	def _negloglikelihood( self , coef ): ##{{{
		if not np.all(self.scale > 0):
			return np.Inf
		
		return np.sum( np.log(self.scale) - self._Y / self.scale )
	##}}}
	
	@AbstractLaw._update_coef
	def _gradient_nlll( self , coef ): ##{{{
		grad_scale = np.zeros_like(coef) + np.nan
		
		pscale = self.params._dparams["scale"]
		if np.all(self.scale > 0):
			grad_scale = pscale.design_.T @ ( ( 1. / self.scale - self._Y / self.scale**2 ) * pscale.gradient() )
		
		return grad_scale.squeeze()
	##}}}


