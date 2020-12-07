# -*- coding: utf-8 -*-

## Copyright(c) 2020 Yoann Robin
## 
## This file is part of SDFC.
## 
## SDFC is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## SDFC is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with SDFC.  If not, see <https://www.gnu.org/licenses/>.


##############
## Packages ##
##############

import numpy as np
import scipy.optimize as sco
import scipy.stats as sc

from .core.__LHS import LHS
from .core.__RHS import RHS


###############
## Class(es) ##
###############

class AbstractLaw:
	##{{{
	"""
	
	Attributes
	==========
	<param> : np.array
		Value of param fitted, can be loc, scale, name of param of law, etc.
	method : string
		method used to fit
	coef_  : numpy.ndarray
		Coefficients fitted
	info_  : AbstractLaw.Info
		Class containing info of the fit
	
	Fit method
	==========
	
	The method <law>.fit is generic, and takes arguments of the form
	<type for param>_<name of param>, see below.
	In case of Bayesian fit, some others optional parameters are available.
	
	Arguments
	---------
	Y         : numpy.ndarray
		Data to fit
	c_<param> : numpy.ndarray or None
		Covariate of a param to fit
	f_<param> : numpy.ndarray or None
		Fix value of a param
	l_<param> : SDFC.tools.LinkFct (optional)
		Link function of a param
	c_global  : list
		List of covariates, sorted by parameters. *_<param> are ignored if set
	l_global  : Inherit of SDFC.link.MultivariateLink
		Global link function. *_<param> are ignored if set
	
	Optional arguments for MLE fit
	------------------------------
	
	Optional arguments for Bayesian fit
	-----------------------------------
	prior : None or law or prior
		Prior for Bayesian fit, if None a Multivariate Normal law assuming
		independence between parameters is used, if you set it, this must be a
		class which implement the method logpdf(coef), returning the log of
		probability density function
	mcmc_init: None or vector of initial parameters
		Starting point of the MCMC algorithm. If None, prior.rvs() is called.
	transition: None or function
		Transition function for MCMC algorithm, if None is given a normal law
		N(0,0.1) is used.
	n_mcmc_drawn : None or integer
		Number of drawn for MCMC algorithm, if None, the value 10000 is used.
	
	Example
	=======
	Example with a Normal law:
	>> _,X_loc,X_scale,_ = SDFC.tools.Dataset.covariates(2500)
	>> loc   = 1. + 0.8 * X_loc
	>> scale = 0.08 * X_scale
	>> 
	>> Y = numpy.random.normal( loc = loc , scale = scale )
	>> 
	>> ## Define the Normal law estimator, with the MLE method and 10 bootstrap for confidence interval:
	>> law = SDFC.Normal( method = "MLE" , n_bootstrap = 10 )
	>>
	>> ## Now perform the fit, c_loc is the covariate of loc, and c_scale the covariate of scale, and we pass a link function to scale:
	>> law.fit( Y , c_loc = X_loc , c_scale = X_scale , l_scale = SDFC.tools.ExpLink() )
	>> print(law) ## That print a summary of fit
	>>
	>> ## But we can assume that scale is stationary, so no covariates are given:
	>> law.fit( Y , c_loc = X_loc , l_scale = SDFC.tools.ExpLink() )
	>> print(law)
	>>
	>> ## Or the loc can be given, so we need to fit only the scale:
	>> law.fit( Y , f_loc = loc , c_scale = X_scale )
	>>
	>> ## And that works for all laws defined in SDFC, you can call
	>> print(law.kinds_params)
	>> ## to print the name of parameters of each law
	"""
	##}}}
	#####################################################################
	
	class _Info(object):##{{{
		def __init__(self):
			self._cov = None
		
		@property
		def cov_(self):
			return self._cov
	##}}}
	
	
	## Init function
	##==============
	
	def __init__( self , names : list , method : str ): ##{{{
		"""
		Initialization of AbstractLaw
		
		Parameters
		----------
		names  : list
			List of names of the law, e.g. ["loc","scale"] for Normal law.
		method : string
			Method called to fit parameters
		
		"""
		self._method = method.lower()
		self._lhs    = LHS(names,0)
		self._rhs    = RHS(self._lhs)
		self.info_   = AbstractLaw._Info()
	##}}}
	
	
	## Properties
	##===========
	
	@property
	def method(self):##{{{
		return self._method
	##}}}
	
	@property
	def coef_(self):##{{{
		return self._rhs.coef_
	##}}}
	
	@coef_.setter
	def coef_( self , coef_ ): ##{{{
		self._rhs.coef_ = coef_
	##}}}
	
	@property
	def cov_(self):##{{{
		return self.info_.cov_
	##}}}
	
	
	## Fit functions
	##==============
	
	def _random_valid_point(self):##{{{
		"""
		Try to find a valid point in the neighborhood of self.coef_
		"""
		coef_ = self.coef_.copy()
		cov_  = 0.1 * np.identity(coef_.size)
		
		p_coef = coef_.copy()
		n_it   = 1
		while not np.isfinite(self._negloglikelihood(p_coef)) or not np.all(np.isfinite(self._gradient_nlll(p_coef))):
			if n_it % 100 == 0: cov_ *= 2
			p_coef = np.random.multivariate_normal( coef_ , cov_ )
			n_it += 1
		self.coef_ = p_coef
	##}}}
	
	def _fit_MLE(self): ##{{{
		
		self._init_MLE()
		self._random_valid_point()
		
		self.info_.mle_optim_result = sco.minimize( self._negloglikelihood , self.coef_ , jac = self._gradient_nlll , method = "BFGS" )
		self.info_._cov = self.info_.mle_optim_result.hess_inv
		self.coef_ = self.info_.mle_optim_result.x
		
	##}}}
	
	def _fit_Bayesian( self , **kwargs ):##{{{
		## Find numbers of features
		##=========================
		n_features = self._rhs.n_features
		
		## Define prior
		##=============
		prior = kwargs.get("prior")
		if prior is None:
			prior = sc.multivariate_normal( mean = np.zeros(n_features) , cov = 10 * np.identity(n_features) )
		
		## Define transition
		##==================
		transition = kwargs.get("transition")
		if transition is None:
			transition = lambda x : x + np.random.normal( size = n_features , scale = 0.1 )
		
		## Define numbers of iterations of MCMC algorithm
		##===============================================
		n_mcmc_drawn = kwargs.get("n_mcmc_drawn")
		if n_mcmc_drawn is None:
			n_mcmc_drawn = 10000
		
		## MCMC algorithm
		##===============
		draw   = np.zeros( (n_mcmc_drawn,n_features) )
		accept = np.zeros( n_mcmc_drawn , dtype = np.bool )
		
		## Init values
		##============
		init = kwargs.get("mcmc_init")
		if init is None:
			init = prior.rvs()
		
		draw[0,:]     = init
		lll_current   = -self._negloglikelihood(draw[0,:])
		prior_current = prior.logpdf(draw[0,:]).sum()
		p_current     = prior_current + lll_current
		
		for i in range(1,n_mcmc_drawn):
			draw[i,:] = transition(draw[i-1,:])
			
			## Likelihood and probability of new points
			lll_next   = - self._negloglikelihood(draw[i,:])
			prior_next = prior.logpdf(draw[i,:]).sum()
			p_next     = prior_next + lll_next
			
			## Accept or not ?
			p_accept = np.exp( p_next - p_current )
			if np.random.uniform() < p_accept:
				lll_current   = lll_next
				prior_current = prior_next
				p_current     = p_next
				accept[i] = True
			else:
				draw[i,:] = draw[i-1,:]
				accept[i] = False
		
		self.coef_ = np.mean( draw[int(n_mcmc_drawn/2):,:] , axis = 0 )
		
		## Update information
		self.info_.draw         = draw
		self.info_.accept       = accept
		self.info_.n_mcmc_drawn = n_mcmc_drawn
		self.info_.rate_accept  = np.sum(accept) / n_mcmc_drawn
		self.info_._cov         = np.cov(draw.T)
	##}}}
	
	def fit( self , Y , **kwargs ): ##{{{
		
		## Add Y
		self._Y = Y.reshape(-1,1)
		
		## Init LHS/RHS
		self._lhs.n_samples = Y.size
		self._rhs.build(**kwargs)
		
		if self._rhs.n_features == 0:
			raise NameError("All parameters are fixed (n_features == 0), no fit")
		
		self.coef_ = np.zeros(self._rhs.n_features)
		## Now fit
		if self._method not in ["mle","bayesian"] and self._rhs.l_global._special_fit_allowed:
			self._special_fit()
		elif self._method == "mle" :
			self._fit_MLE()
		else:
			self._fit_Bayesian(**kwargs)
		
	##}}}


