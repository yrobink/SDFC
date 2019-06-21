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

class NormalLaw(AbstractLaw):
	"""
	SDFC.NormalLaw
	==============
	
	Fit parameters of a Normal law, possibly with co-variable
	
	Attributes
	----------
	method : string
		method used to fit
	loc    : numpy.ndarray
		Location fitted
	scale  : numpy.ndarray
		Scale fitted
	coef_  : numpy.ndarray
		Coefficients fitted
	n_bootstrap: integer
		Numbers of bootstrap for confidence interval
	coefs_bootstrap: numpy.ndarray
		coef_ for each bootstrap
	confidence interval: numpy.ndarray[ shape = (2,coef_.size) ]
		Confidence interval, first line is the alpha/2 quantile, and second line the 1 - alpha/2 quantile
	alpha          : float
		Level of confidence interval
	"""
	
	def __init__( self , method = "MLE" , link_fct_loc = IdLinkFct() , link_fct_scale = IdLinkFct() , n_bootstrap = 0 , alpha = 0.05 ): ##{{{
		"""
		Initialization of NormalLaw
		
		Parameters
		----------
		method         : string
			Method called to fit parameters, options are "moments" and "MLE" (Maximum Likelihood estimation)
		link_fct_loc   : a class herited from SDFC.tools.LinkFct
			Link function for loc, default is IdLinkFct()
		link_fct_scale : a class herited from SDFC.tools.LinkFct
			Link function for scale, default is SDFC.tools.IdLinkFct(). Interesting option is SDFC.tools.ExpLinkFct().
		n_bootstrap    : integer
			Numbers of bootstrap for confidence interval, default = 0 (no bootstrap)
		alpha          : float
			Level of confidence interval, default = 0.05
		
		"""
		AbstractLaw.__init__( self , method , n_bootstrap , alpha )
		
		self.loc       = None
		self.scale     = None
		
		self._loc   = LawParam( linkFct = link_fct_loc   , kind = "loc"   )
		self._scale = LawParam( linkFct = link_fct_scale , kind = "scale" )
		self._lparams = [self._loc,self._scale]
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.NormalLaw\n"
		val += "--------------\n"
		val += "* fit_method : {}\n".format(self.method)
		val += "* link_loc   : {}\n".format(str(self._loc.linkFct))
		val += "* link_scale : {}\n".format(str(self._scale.linkFct))
		val += "* coef_      : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def fit( self , Y , loc_cov = None , scale_cov = None , floc = None , fscale = None ): ##{{{
		"""
		Fit function for NormalLaw
		
		Arguments
		---------
		
		Y         : numpy.ndarray
			Data to fit
		loc_cov   : None or numpy.ndarray
			Co-variates of loc in columns.
		scale_cov : None or numpy.ndarray
			Co-variates of scale in columns.
		floc      : None or numpy.ndarray
			If not None, fix the value of loc parameter (so not fitted)
		fscale    : None or numpy.ndarray
			If not None, fix the value of scale parameter (so not fitted)
		"""
		self._size = Y.size
		
		
		if self.n_bootstrap > 0:
			self.coefs_bootstrap = []
			if loc_cov is not None and loc_cov.ndim == 1:
				loc_cov = loc_cov.reshape( (loc_cov.size,1) )
			if scale_cov is not None and scale_cov.ndim == 1:
				scale_cov = scale_cov.reshape( (scale_cov.size,1) )
			if np.isscalar(floc):
				floc = np.array([floc]).ravel()
			if np.isscalar(fscale):
				fscale = np.array([fscale]).ravel()
			
			
			for i in range(self.n_bootstrap):
				idx = np.random.choice( self._size , self._size )
				Y_bs         = Y[idx]
				loc_cov_bs   = None   if loc_cov   is None else loc_cov[idx,:]
				scale_cov_bs = None   if scale_cov is None else scale_cov[idx,:]
				floc_bs      = floc   if floc      is None or floc.size == 1   else floc[idx]
				fscale_bs    = fscale if fscale    is None or fscale.size == 1 else fscale[idx]
				
				self._fit( Y_bs , loc_cov_bs , scale_cov_bs , floc_bs , fscale_bs )
				self.coefs_bootstrap.append( self.coef_ )
			
			self.coefs_bootstrap = np.array( self.coefs_bootstrap )
			self.confidence_interval = np.quantile( self.coefs_bootstrap , [ self.alpha / 2. , 1 - self.alpha / 2.] , axis = 0 )
		
		self._fit( Y , loc_cov , scale_cov , floc , fscale )
		
	##}}}
	
	def _fit( self , Y , loc_cov = None , scale_cov = None , floc = None , fscale = None ): ##{{{
		self._Y    = np.ravel(Y)
		self._size = Y.size
		
		self._loc.init(   X = loc_cov   , fix_values = floc   , size = self._size )
		self._scale.init( X = scale_cov , fix_values = fscale , size = self._size )
		
		if self.method == "moments":
			self._fit_moments()
		else:
			self._fit_mle()
		
		self.coef_ = self._concat_param()
		
	##}}}
	
	def _fit_moments(self):##{{{
		
		if self._loc.not_fixed():
			lX = self._loc.design_wo1()
			self._loc.set_coef( mean( self._Y , lX , return_coef = True , linkFct = self._loc.linkFct ) )
			self._loc.update()
		
		if self._scale.not_fixed():
			sX = self._scale.design_wo1()
			self._scale.set_coef( std(  self._Y , sX , m = self._loc.valueLf() , return_coef = True , linkFct = self._scale.linkFct ) )
			self._scale.update()
		
		self.loc   = self._loc.valueLf()
		self.scale = self._scale.valueLf()
	##}}}
	
	def _fit_mle(self):##{{{
		
		self._fit_moments()
		param_init = self._concat_param()
		self.optim_result = sco.minimize( self._optim_function , param_init , jac = self._gradient_optim_function , method = "BFGS" )
		self._update_param( self.optim_result.x )
		
	##}}}
	
	def _split_param( self , param ):##{{{
		param_loc   = None
		param_scale = None
		
		if self._loc.not_fixed():
			param_loc   = param[:self._loc.size]
			if self._scale.not_fixed():
				param_scale = param[self._loc.size:]
		elif self._scale.not_fixed():
			param_scale = param[:self._scale.size]
		
		return param_loc,param_scale
	##}}}
	
	def _concat_param( self ):##{{{
		param = None
		if self._loc.not_fixed() and self._scale.not_fixed():
			param = np.hstack( (self._loc.coef_,self._scale.coef_) )
		elif self._loc.not_fixed():
			param = self._loc.coef_
		elif self._scale.not_fixed():
			param = self._scale.coef_
		return param
	##}}}
	
	def _negloglikelihood( self ): ##{{{
		scale2 = np.power( self.scale , 2 )
		return np.Inf if not np.all( self.scale > 0 ) else np.sum( np.log( scale2 ) ) / 2. + np.sum( np.power( self._Y - self.loc , 2 ) / scale2 ) / 2.
	##}}}
	
	def _update_param( self , param ):##{{{
		param_loc,param_scale = self._split_param(param)
		self._loc.set_coef(   param_loc   )
		self._scale.set_coef( param_scale )
		self._loc.update()
		self._scale.update()
		
		self.loc   = self._loc.valueLf()
		self.scale = self._scale.valueLf()
	##}}}
	
	def _optim_function( self , param ):##{{{
		self._update_param(param)
		return self._negloglikelihood()
	##}}}
	
	def _gradient_optim_function( self , param ): ##{{{
		self._update_param(param)
		
		grad = np.array( [] )
		
		Yc = self._Y - self.loc
		if self._loc.not_fixed():
			grad_loc   = - self._loc.design_.T @ (Yc / self.scale**2 * self._loc.valueGrLf() )
			grad = np.hstack( (grad,grad_loc) )
		if self._scale.not_fixed():
			grad_scale = self._scale.design_.T @ ( ( 1. / self.scale - Yc**2 / self.scale**3 ) * self._scale.valueGrLf() )
			grad = np.hstack( (grad,grad_scale) )
		
		return grad
	##}}}




