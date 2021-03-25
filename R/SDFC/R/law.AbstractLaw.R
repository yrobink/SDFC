
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



#' AbstractLaw 
#'
#' @description
#' Base class of parametric laws
#'
#' @details
#' Just to be inherit
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
AbstractLaw = R6::R6Class( "AbstractLaw" ,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	## Private fields ##{{{
	
	#' @field method [string] Method used to fit the distribution
	.method = NULL,
	#' @field rhs [SDFC::RHS] Right Hand Side
	.rhs = NULL,
	#' @field coef_ [numeric] Coefficients fitted
	.Y = NULL,
	#' @field lhs [SDFC::RHS] Left Hand Side
	#' @field cov [matrix] Covariance matrix of coefs
	
	##}}}
	
	random_valid_point = function() ##{{{
	{
		coef_ = self$coef_
		scale = 0.1
		
		p_coef = coef_
		n_it = 1
		
		while( !(is.finite(private$negloglikelihood(p_coef)) ) )# && base::all(is.finite(private$gradient_nlll(p_coef)))) )
		{
			if( n_it %% 100 == 0 ) scale = scale * 2
			
			p_coef = coef_ + stats::rnorm( n = length(coef_) , mean = 0 , sd = scale )
			n_it = n_it + 1
		}
		
		self$coef_ = p_coef
	},
	##}}}
	
	fit_MLE = function(...) ##{{{
	{
		kwargs = list()
		private$init_MLE()
		coef_ = self$coef_
		
		is.success = FALSE
		n_test     = 0
		max_test   = kwargs[["mle_n_restart"]]
		if( is.null(max_test) ) max_test = 2
		
		while( !is.success && n_test < max_test )
		{
			private$random_valid_point()
			self$info_[["mle_optim_result"]] = stats::optim( self$coef_ , private$negloglikelihood , method = "BFGS" , hessian = TRUE )
			self$coef_ = self$info_[["mle_optim_result"]]$par
			is.success = self$info_[["mle_optim_result"]]$convergence == 0
			n_test = n_test + 1
		}
		self$info_[["cov"]] = base::solve(self$info_[["mle_optim_result"]]$hessian)
		
		if( !is.success )
			self$coef_ = coef_
		
	},
	##}}}
	
	fit_Bayesian = function(...)##{{{
	{
		kwargs = list(...)
	}
	##}}}
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	method = function(value)
	{
		if(missing(value))
		{
			return(private$.method)
		}
		else
		{
			private$.method = base::tolower(value)
		}
	},
	
	lhs = function(value)
	{
		return(private$.rhs$lhs)
	},
	
	rhs = function(value)
	{
		return(private$.rhs)
	},
	
	coef_ = function(value)
	{
		if(missing(value))
		{
			return(private$.rhs$coef_)
		}
		else
		{
			private$.rhs$coef_ = value
		}
	},
	
	cov = function(value)
	{
		if(missing(value))
			return(self$info_$cov)
	}
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	
	## Public fields ##{{{
	
	#' @field info_ [list] List of information about the fit
	info_ = NULL,
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new AbstractLaw object.
    #' @param names [vector(string)] names of parameters
    #' @param method [string] method used to fit
	#' @return A new `AbstractLaw` object.
	initialize = function( names , method )
	{
		self$method  = method
		lhs = LHS$new(names,0)
		private$.rhs = RHS$new(lhs)
		self$info_   = list()
	},
	
	##}}}
	
	## Fit ##{{{
	
	#' @description
    #' Fit the law.
    #' @param Y [numeric] Data to fit
    #' @param ... Many arguments
	fit = function( Y , ... ) 
	{
		kwargs = list(...)
		
		## Add Y
		private$.Y = matrix( Y , ncol = 1 )
		
		## Init LHS/RHS
		private$.rhs$lhs$n_samples = length(Y)
		base::do.call( private$.rhs$build , kwargs )
		
		self$coef_ = numeric(private$.rhs$n_features)
		
		## Now fit
		if( !(self$method %in% base::c("mle","bayesian")) && private$.rhs$l_global$special_fit_allowed )
		{
			private$special_fit()
		}
		else if( self$method == "mle" )
		{
			base::do.call( private$fit_MLE , kwargs )
		}
		else
		{
			base::do.call( private$fit_Bayesian , kwargs )
		}
		
	}
	##}}}
	
	## Fit bootstrap ##{{{
	
	##}}}
	
	)
	
	##}}}
	
)


