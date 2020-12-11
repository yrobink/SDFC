
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


#' QuantileRegression
#'
#' Class to perform a Quantile Regression with covariates
#'
#' @docType class
#'
#' @description
#' Class for quantile regression
#'
#' @details
#' Wrapper which call c++ quantile regression
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ## Data
#' size = 2000
#' c_data = SDFC::dataset(size)
#' 
#' loc   = 0.5 + 2 * c_data$X_loc
#' scale = 1 + 2 * c_data$X_scale
#' Y = stats::rnorm( size , mean = loc , sd = scale )
#' 
#' ## Quantile regression
#' ltau  = base::seq( 0.01 , 0.99 , length = 100 )
#' qr = SDFC::QuantileRegression$new(ltau)
#' qr$fit( Y , base::cbind( c_data$X_loc , c_data$X_scale ) )
#' Yq = if( qr$is_success ) qr$predict() else NULL
#' @export
QuantileRegression = R6::R6Class( "QuantileRegression" ,
	
	## Private list
	##=============
	## {{{
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	qrmethod      = NULL,
	#' @field is_fitted [bool] If QR is fitted
	.is_fitted    = NULL,
	#' @field is_success [bool] If QR fit is a success
	.is_success   = NULL,
	#' @field is_unfeasible [bool] If QR fit is unfeasible
	.is_unfeasible = NULL,
	#' @field coef_ [vector] Value of coef fitted
	.coef_ = NULL
	
	#############
	## Methods ##
	#############
	
	
	),
	##}}}
	
	## Active list
	##=============
	## {{{
	
	active = list(
	
	coef_ = function(value) ##{{{
	{
		if(missing(value))
		{
			return( private$qrmethod$coef() )
		}
	},
	##}}}
	
	is_fitted = function(value) ##{{{
	{
		if(missing(value))
		{
			return( private$qrmethod$is_fitted() )
		}
	},
	##}}}
	
	is_success = function(value) ##{{{
	{
		if(missing(value))
		{
			return( private$qrmethod$is_success() )
		}
	},
	##}}}
	
	is_unfeasible = function(value) ##{{{
	{
		if(missing(value))
		{
			return( private$qrmethod$is_unfeasible() )
		}
	}
	##}}}
	
	),
	##}}}
	
	## Public list
	##============
	## {{{
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#' @field ltau [vector] vector of probabilities
	ltau    = NULL,
	#' @field maxit [integer] Max number of iterations for fit
	maxit   = NULL,
	#' @field tol [float] Numerical tolerance
	tol     = NULL,
	#' @field beta [float] Beta parameters
	beta    = NULL,
	#' @field method [string] Method used for fit
	method  = NULL,
	#' @field verbose [bool] Print or not message for fit
	verbose = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	## initialize ##{{{
	#' @description
    #' Create a new QuantileRegression object.
	#' @param ltau    [vector] vector of probabilities
	#' @param method  [string] Method used for fit
	#' @param verbose [bool]   Print or not message for fit
	#' @return A new `QuantileRegression` object.
	initialize = function( ltau , method = "Frish-Newton" , verbose = FALSE )
	{
		self$ltau        = ltau
		self$maxit       = 50
		self$tol         = 1e-6
		self$beta        = 0.99995
		self$method      = method
		self$verbose     = verbose
		private$qrmethod = FrishNewton$new( ltau , self$maxit , self$tol , self$beta )
	},
	##}}}
	
	
	#############
	## Methods ##
	#############
	
	## fit ##{{{
	#' @description
    #' Fit the quantile regression
    #' @param Y [vector] Dataset to fit
    #' @param X [matrix] Covariates
    #' @return NULL
	fit = function( Y , X )
	{
		if( !is.matrix(X) )
		{
			X = matrix( X , nrow = length(X) , ncol = 1 ) 
		}
		private$qrmethod$fit( Y , X )
		if( private$qrmethod$is_unfeasible() && self$verbose )
		{
			print( "SDFC::QuantileRegression : Unfeasible problem" )
		}
	},
	##}}}
	
	## predict ##{{{
	#' @description
    #' Return the quantile fitted
    #' @return [vector] quantile fitted
	predict = function()
	{
		return( private$qrmethod$predict() )
	}
	##}}}
	
	)
	
	##}}}
	
	
	
)

