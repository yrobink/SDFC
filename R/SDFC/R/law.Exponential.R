
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


#' Exponential
#'
#' @description
#' Exponential distribution
#'
#' @details
#' Class to fit a Exponential law with covariates, available methods are:
#' 
#' moments  : use empirical estimator of mean and standard deviation to find
#'            loc and scale, possibly with least square regression if
#'            covariates are given
#'
#' bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of
#'            n_mcmc_iteration sample draw from the posterior P(coef_ | Y)
#'
#' mle      : Maximum likelihood estimation
#'
#' 
#' Parameters:
#'
#' scale : scale parameter
#'
#' See AbstractLaw for details and examples
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
Exponential = R6::R6Class( "Exponential" ,
	
	inherit = AbstractLaw,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	fit_moments = function() ##{{{
	{
		X = self$rhs$c_global[["scale"]]
		self$coef_ = SDFC::mean( private$.Y , X , private$.rhs$l_global$linkc("scale")$l , value = FALSE )
	},
	##}}}
	
	special_fit = function() ##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
	},
	##}}}
	
	init_MLE = function()##{{{
	{
		if( private$.rhs$l_global$special_fit_allowed )
		{
			private$fit_moments()
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
		
		return( base::sum( base::log(self$scale) - private$.Y / self$scale ) )
	},
	##}}}
	
	gradient_nlll = function(coef) ##{{{
	{
		self$coef_ = coef
		
		
		if( !base::all(self$scale > 0) )
		{
			return( numeric(length(self$coef_)) + NA )
		}
		
		jac = private$.rhs$lhs$jacobian
		T0 = 1. / self$scale - private$.Y / self$scale^2
		jac[1,,] = jac[1,,] * as.vector(T0)
		
		
		return( base::apply( jac , 3 , base::sum ) )
	}
	##}}}
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	scale = function(value)
	{
		return(private$.rhs$lhs$values$scale)
	}
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	
	## Public fields ##{{{
	
	#' @field scale [numeric] Scale parameter
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new Exponential object.
    #' @param method [string] method used to fit
	#' @return A new `Exponential` object.
	initialize = function( method = "mle" )
	{
		super$initialize( base::c("scale") , method )
	}
	
	##}}}
	
	)
	
	##}}}
	
)




