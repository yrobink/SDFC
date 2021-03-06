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
import texttable      as tt
from SDFC.tools.__LawParams import LawParams


###########
## Class ##
###########

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
	coefs_bootstrap: numpy.ndarray
		coef_ for each bootstrap
	confidence interval: numpy.ndarray[ shape = (2,coef_.size) ]
		Confidence interval, first line is the alpha/2 quantile, and second line the 1 - alpha/2 quantile
	alpha          : float
		Level of confidence interval
	
	
	Fit method
	==========
	
	The method <law>.fit is generic, and takes arguments of the form <type for param>_<name of param>, see below.
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
	
	Optional arguments for MLE fit
	------------------------------
	
	Optional arguments for Bayesian fit
	-----------------------------------
	prior : None or law or prior
		Prior for Bayesian fit, if None a Multivariate Normal law assuming independence between parameters is used,
		if you set it, this must be a class which implement the method logpdf(coef), returning the log of probability
		density function
	mcmc_init: None or vector of initial parameters
		Starting point of the MCMC algorithm. If None, prior.rvs() is called.
	transition: None or function
		Transition function for MCMC algorithm, if None is given a normal law N(0,0.1) is used.
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
	
	class _Bootstrap:##{{{
		def __init__( self , n_bootstrap , alpha ):##{{{
			self.n_bootstrap         = n_bootstrap
			self.coefs_bootstrap     = None
			self.confidence_interval = None
			self.alpha               = alpha
		##}}}
		
		def _bootstrap_method(func):##{{{
			def wrapper(*args,**kwargs):
				self,Y = args
				if self.n_bootstrap > 0:
					self.coefs_bootstrap = []
					for _ in range(self.n_bootstrap):
						idx = np.random.choice( Y.size , Y.size , replace = True )
						self.params = LawParams( kinds = self.kinds_params )
						self.params.add_params( n_samples = Y.size , resample = idx , **kwargs )
						self._Y = Y.reshape(-1,1)[idx,:]
						if self.method not in [ "mle" , "bayesian" ]:
							self._fit()
						elif self.method == "bayesian":
							self._fit_bayesian()
						else:
							self._fit_mle()
						self.coefs_bootstrap.append( self.coef_ )
					self.coefs_bootstrap = np.array( self.coefs_bootstrap )
					self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
				return func(*args,**kwargs)
			return wrapper
		##}}}
	##}}}
	
	@property
	def n_bootstrap(self):##{{{
		return self._bootstrap.n_bootstrap
	##}}}
	
	@n_bootstrap.setter
	def n_bootstrap( self , n_bootsrap ):##{{{
		self._bootstrap.n_bootstrap = n_bootstrap
	##}}}
	
	@property
	def coefs_bootstrap(self):##{{{
		return self._bootstrap.coefs_bootstrap
	##}}}
	
	@coefs_bootstrap.setter
	def coefs_bootstrap( self , coefs_bootstrap ):##{{{
		self._bootstrap.coefs_bootstrap = coefs_bootstrap
	##}}}
	
	@property
	def confidence_interval(self):##{{{
		return self._bootstrap.confidence_interval
	##}}}
	
	@confidence_interval.setter
	def confidence_interval( self , confidence_interval ):##{{{
		self._bootstrap.confidence_interval = confidence_interval
	##}}}
	
	@property
	def alpha(self):##{{{
		return self._bootstrap.alpha
	##}}}
	
	@alpha.setter
	def alpha( self , alpha ):##{{{
		self._bootstrap.alpha = alpha
	##}}}
	
	
	#####################################################################
	
	class _Info:##{{{
		def __init__( self ):
			self.cov = None
	##}}}
	
	@property
	def cov(self):##{{{
		return self._info.cov
	##}}}
	
	@property
	def info(self):##{{{
		return self._info
	##}}}
	
	#####################################################################
	
	def __init__( self , kinds_params , method , n_bootstrap , alpha ):##{{{
		"""
		Initialization of AbstractLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		
		"""
		self.method    = method.lower()
		self.params    = {}
		self._kinds_params = kinds_params
		
		self._bootstrap = AbstractLaw._Bootstrap( n_bootstrap , alpha )
		self._info      = AbstractLaw._Info()
		
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.AbstractLaw\n"
		val += "----------------\n"
		val += "Base class, if you read this message, you have done a mistake"
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def _to_str( self ):##{{{
		tab = tt.Texttable( max_width = 0 )
		
		## Header
		header = [ str(type(self)).split(".")[-1][:-2] + " ({})".format(self.method) , "Link" , "Type" , "coef" ]
		if self.confidence_interval is not None:
			header += [ "Quantile {}".format(self.alpha/2) , "Quantile {}".format( 1 - self.alpha / 2 ) ]
		tab.header( header )
		
		## Loop on params
		a = 0
		for k in self.params._dparams:
			p = self.params._dparams[k]
			coef = None if p.coef_ is None else str(p.coef_.squeeze().round(3).tolist())
			row = [ p.kind , str(p.link).split(".")[-1] , str(type(p)).split(".")[-1].split("Param")[0] , coef ]
			if self.confidence_interval is not None:
				if not p.is_fix():
					b = a + p.n_features
					row += [ self.confidence_interval[i,a:b].squeeze().round(3).tolist() for i in range(2) ]
					a = b
				else:
					row += [ str(None) , str(None) ]
			tab.add_row( row )
		return tab.draw() + "\n"
	##}}}
	
	@property
	def kinds_params(self):##{{{
		return self._kinds_params
	##}}}
	
	@kinds_params.setter
	def kinds_params( self , kp ):##{{{
		pass
	##}}}
	
	@property
	def coef_(self):##{{{
		return self.params.coef_
	##}}}
	
	@coef_.setter
	def coef_( self , coef_  ):##{{{
		self.params.update_coef(coef_)
	##}}}
	
	def _predict_covariate( self , kind , c_p ): ##{{{
		p = self.params._dparams[kind]
		
		if isinstance(p,CovariateParam) and c_p is not None:
			return p.coef_[0] + c_p.T @ p.coef_[1:]
		return p.value
	##}}}
	
	def _update_coef(func):##{{{
		def wrapper(*args):
			args[0].params.update_coef(args[1])
			return func(*args)
		return wrapper
	##}}}
	
	@_Bootstrap._bootstrap_method
	def fit( self , Y , **kwargs ): ##{{{
		
		## Fit part
		self.params = LawParams( kinds = self.kinds_params )
		self.params.add_params( n_samples = Y.size , resample = None , **kwargs )
		self._Y = Y.reshape(-1,1)
		if self.method not in ["mle","bayesian"]:
			self._fit()
		elif self.method == "bayesian":
			self._fit_bayesian(**kwargs)
		else:
			self._fit_mle()
		del self._Y
	##}}}
	
	def _fit_bayesian( self , **kwargs ):##{{{
		
		## Find numbers of features
		##=========================
		n_features = 0
		for k in self.params._dparams:
			n_features += self.params._dparams[k].n_features
		
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
		draw = np.zeros( (n_mcmc_drawn,n_features) )
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
		
		self.params.update_coef( np.mean( draw , axis = 0 ) )
		
		## Update information
		self._info.draw         = draw
		self._info.accept       = accept
		self._info.n_mcmc_drawn = n_mcmc_drawn
		self._info.rate_accept  = np.sum(accept) / n_mcmc_drawn
		self._info.cov          = np.cov(draw.T)
	##}}}
	
	def _fit_mle( self ):##{{{
		self._initialization_mle()
		self._info.optim_result = sco.minimize( self._negloglikelihood , self.coef_ , jac = self._gradient_nlll , method = "BFGS" )
		self.coef_ = self._info.optim_result.x
		self._info.cov = self._info.optim_result.hess_inv
	##}}}



