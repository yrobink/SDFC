
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


#' GPD
#'
#' @description
#' GPD distribution
#'
#' @details
#' Class to fit a GPD law with covariates, available methods are:
#' 
#' moments  : use empirical estimator of mean and standard deviation to find
#'            loc and scale, possibly with least square regression if
#'            covariates are given
#'
#' lmoments : Use L-Moments estimation, only in stationary context
#'
#' lmoments_experimental: Use non-stationary L-Moments with Quantile
#'            Regression, experimental and not published, only
#'            used to find an initialization of MLE
#'
#' bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of
#'            n_mcmc_iteration sample draw from the posterior P(coef_ | Y)
#'
#' mle      : Maximum likelihood estimation
#'
#' Warnings: For this class, f_loc must be always given, because we fit a pareto
#'           beyond the loc parameter!
#' 
#' Parameters:
#'
#' loc   : location parameter
#'
#' scale : scale parameter
#'
#' shape : shape parameter
#'
#' See AbstractLaw for details and examples
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
GPD = R6::R6Class( "GPD" ,
	
	inherit = AbstractLaw,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	fit_moments = function() ##{{{
	{
		coef_ = numeric(length(self$coef_))
		
		
		self$coef_ = coef_
	},
	##}}}
	
	fit_lmoments = function() ##{{{
	{
		coef_ = numeric(length(self$coef_))
		
		il_b  = 1
		il_e  = il_b + private$.rhs$l_global$sizec("loc") - 1
		isc_b = il_e + 1
		isc_e = isc_b + private$.rhs$l_global$sizec("scale") - 1
		ish_b = isc_e + 1
		ish_e = ish_b + private$.rhs$l_global$sizec("shape") - 1
		
		Z = private$.Y - self$loc
		
		lmom = lmoments( Z )
		if( !private$.rhs$lhs$fixed[["scale"]] && !private$.rhs$lhs$fixed[["shape"]] )
		{
			itau     = lmom[1] / lmom[2]
			scale_lm = lmom[1] * ( itau - 1 )
			scale_lm = if( scale_lm > 0) scale_lm else 1e-8
			shape_lm = 2 - itau
			coef[isc_b] = scale_lm
			coef[ish_b] = shape_lm
		}
		else if( !private$.rhs$lhs$fixed[["scale"]] )
		{
			scale = lmom[1] * ( 1 - self.shape )
			if( base::any(!(scale > 0)) )
				scale[ !(scale > 0) ] = 1e-8
			coef[isc_b:isc_e] = SDFC::mean( scale , self$rhs$c_global[["scale"]] , value = FALSE , link =  private$.rhs$l_global$linkc("scale")$l )
		}
		else if( !private$.rhs$lhs$fixed[["shape"]] )
		{
			Y = Z / self$scale
			lmom = lmoments(Y)
			itau     = lmom[1] / lmom[2]
			coef[ish_b] = 2 - itau
		}
		
		self$coef_ = coef_
	},
	##}}}
	
	fit_lmoments_experimental = function() ##{{{
	{
		self$coef_ = coefs
	},
	##}}}
	
	special_fit = function() ##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
		else if( self$method == "lmoments" )
			private$fit_lmoments()
		else if( self$method == "lmoments-experimental" )
			private$fit_lmoments_experimental()
	},
	##}}}
	
	init_MLE = function()##{{{
	{
		if( private$.rhs$l_global$special_fit_allowed )
		{
			private$fit_lmoments_experimental()
		}
		else
		{
			self$coef_ = private$.rhs$l_global$valid_point(self)
		}
	},
	##}}}
	
	negloglikelihood = function(coef) ##{{{
	{
		self$coef_ = coef
		
		if( !base::all(self$scale > 0) )
			return(Inf)
		
		## Remove exponential case
		zero_shape = ( base::abs(self$shape) < 1e-10 )
		shape = self$shape
		if( base::any(zero_shape) )
			shape[zero_shape] = 1e-10
		
		loc   = self$loc
		scale = self$scale
		
		##
		Z = 1 + shape * ( private$.Y - loc ) / scale
		
		if( !base::all(Z > 0) )
			return(Inf)
		
		
		res = base::sum( base::log( scale ) + base::log(Z) * ( 1 + 1. / shape ) )
		
		if( is.finite(res) )
			return(res)
		
		return(Inf)
	},
	##}}}
	
	gradient_nlll = function(coef) ##{{{
	{
		self$coef_ = coef
		
		loc   = self$loc
		scale = self$scale
		shape = self$shape
		
		Z  = ( private$.Y - loc ) / scale
		ZZ = 1. + shape * Z
		C  = ( 1. + 1. / shape ) * Z / ZZ
		
		T0 = 1 / scale - C * shape / scale
		T1 = - base::log(ZZ) / shape^2 + C
		
		
		for( T in list(T0,T1) )
		{
			if( !base::all(is.finite(T)) )
			{
				return( numeric(length(self$coef_)) + NA )
			}
		}
		
		
		## Compute gradient
		jac = private$.rhs$lhs$jacobian
		p   = 1
		if( !private$.rhs$lhs$fixed[["scale"]] )
		{
			jac[p,,] = jac[p,,] * as.vector(T0)
		}
		if( !private$.rhs$lhs$fixed[["shape"]] )
		{
			jac[p,,] = jac[p,,] * as.vector(T1)
		}
		
		return( base::apply( jac , 3 , base::sum ) )
	}
	##}}}
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	loc = function(value)
	{
		return(private$.rhs$lhs$values$loc)
	},
	
	scale = function(value)
	{
		return(private$.rhs$lhs$values$scale)
	},
	
	shape = function(value)
	{
		return(private$.rhs$lhs$values$shape)
	}
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	
	## Public fields ##{{{
	
	#' @field loc [numeric] Location parameter
	#' @field scale [numeric] Scale parameter
	#' @field shape [numeric] Scale parameter
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new GPD object.
    #' @param method [string] method used to fit
	#' @return A new `GPD` object.
	initialize = function( method = "mle" )
	{
		super$initialize( base::c("loc","scale","shape") , method )
	}
	
	##}}}
	
	)
	
	##}}}
	
)




