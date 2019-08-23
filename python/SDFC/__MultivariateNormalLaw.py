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

class MultivariateNormalLaw:
	
	def __init__( self , method = "MLE" ):##{{{
		self.d      = -1
		self._Y     = None
		self.method = method
		
		self.mean  = None
		self.cov   = None
		self._mean = LawParam( kind = "mean" )
		self._cov  = LawParam( kind = "cov" )
		
		self.coef_ = None
		
		self.optim_result = None
	##}}}
	
	def __str__(self):##{{{
		val  = "SDFC.MultivariateNormalLaw\n"
		val += "--------------\n"
		val += "* fit_method : {}\n".format(self.method)
		val += "* link_mean  : {}\n".format(str(self._mean.linkFct))
		val += "* link_cov   : {}\n".format(str(self._cov.linkFct))
		val += "* coef_      : {}\n".format(self.coef_)
		
		return val
	##}}}
	
	def __repr__(self):##{{{
		return self.__str__()
	##}}}
	
	def fit( self , Y , mean_cov = None , cov_cov = None ):##{{{
		self._Y   = Y
		self.d    = Y.shape[1]
		self._idx = np.triu_indices(self.d)
		
		self._mean.init( X = mean_cov , dim = self.d , size = self._Y.shape[0] )
		self._cov.init(  X = cov_cov  , dim = int(self.d + self.d * (self.d - 1) / 2) , size = self._Y.shape[0] , transform = self._cov_transform )
		
		if self.method == "moments":
			self._fit_moments()
		else:
			self._fit_mle()
		
		self.coef_ = self._concat_param()
	##}}}
	
	def predict_mean( self , mean_cov = None ):##{{{
		return self._predict_param( self._mean , mean_cov )
	##}}}
	
	def predict_cov( self , cov_cov = None ):##{{{
		return self._predict_param( self._cov , cov_cov )
	##}}}
	
	def _predict_param( self , param , cov ):##{{{
		if cov is None or param.shape[0] == 1:
			return param.valueLf()
		else:
			if cov.ndim == 1: cov = cov.reshape( (cov.size,1) )
			design = np.hstack( (np.ones((cov.shape[0],1)),cov) )
			return param._transform( design @ param.coef_[:design.shape[1]] )
	##}}}
	
	def _cov_transform( self , M ):##{{{
		cov = np.zeros( (M.shape[0],self.d,self.d) )
		cov_tmp = np.zeros((self.d,self.d))
		for i in range(M.shape[0]):
			cov_tmp[self._idx]   = M[i,:]
			cov_tmp.T[self._idx] = M[i,:]
			cov[i,:,:] = cov_tmp
		return cov
	##}}}
	
	def _split_param( self , param ):##{{{
		param_mean = None
		param_cov  = None
		
		if self._mean.not_fixed():
			param_mean   = param[:self._mean.size]
			if self._cov.not_fixed():
				param_cov = param[self._mean.size:]
		elif self._cov.not_fixed():
			param_cov = param
		
		return param_mean,param_cov
	##}}}
	
	def _concat_param( self ):##{{{
		param = np.array([])
		for p in [self._mean,self._cov]:
			if p.not_fixed():
				param = np.hstack( (param,p.coef_.ravel()) )
		return param
	##}}}
	
	def _update_param( self , param ):##{{{
		param_mean,param_cov = self._split_param(param)
		self._mean.set_coef( param_mean )
		self._cov.set_coef(  param_cov  )
		self._mean.update()
		self._cov.update()
		
		self.mean = self._mean.valueLf()
		self.cov  = self._cov.valueLf()
	##}}}
	
	def _fit_moments(self):##{{{
		
		if self._mean.not_fixed():
			lX = self._mean.design_wo1()
			coef_mean = np.zeros_like( self._mean.coef_ )
			for i in range(self.d):
				coef_mean[:,i] = sdnp.mean( self._Y[:,i] , lX , return_coef = True )
			self._mean.set_coef( coef_mean )
			self._mean.update()
			self.mean = self._mean.valueLf()
		
		if self._cov.not_fixed():
			lX = self._cov.design_wo1()
			coef_cov = np.zeros( (self.d,self.d,self._cov.shape[0]) )
			for i in range(self.d):
				for j in range(i,self.d):
					coef_cov[i,j,:] = sdnp_cov( self._Y[:,i] , self._Y[:,j] , X = lX , m0 = self.mean[:,i] , m1 = self.mean[:,j] , return_coef = True )
			for k,ij in enumerate(zip(*self._idx)):
				i,j = ij
				self._cov.coef_[:,k] = coef_cov[i,j,:]
		
		self._update_param( self._concat_param() )
	##}}}
	
	def _fit_mle(self):##{{{
		self._fit_moments()
		param_init = self._concat_param()
		self.optim_result = sco.minimize( self._optim_function , x0 = param_init )
		self._update_param( self.optim_result.x )
	##}}}
	
	def _negloglikelihood( self ):##{{{
		
		detCov = np.array( [np.linalg.det(self.cov[i,:,:]) for i in range(self.cov.shape[0])] )
		
		if not np.all(detCov > 0):
			return np.Inf
		
		icov = np.linalg.inv(self.cov)
		
		res = np.sum( np.log(detCov) )
		for i in range(self._Y.shape[0]):
			z = self._Y[i,:] - self.mean[i,:]
			res += z @ np.linalg.solve( self.cov[i,:,:] , z )
		
		return res / 2
	##}}}
	
	def _optim_function( self , param ):##{{{
		self._update_param(param)
		return self._negloglikelihood()
	##}}}
	
	def _gradient_optim_function( self , param ):##{{{
		self._update_param(param)
		
		grad = np.array([])
		
		if self._mean.not_fixed():
			pass
		
		if self._cov.not_fixed():
			pass
		
		return grad
	##}}}


