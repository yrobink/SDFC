
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
#' @importFrom R6 R6Class
#'
#' @param ltau [vector]
#'        Vector of quantiles where we want to perform the fit
#' @param method [string]
#'        Method used to fit, currently, only "Frish-Newton" is available
#' @param verbose [bool]
#'        Print warning and error message
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(ltau,method,verbose)}}{Initialize Quantile Regression with code{QuantileRegression}}
#'   \item{\code{fit(Y,X)}}{Fit the quantile regression}.
#'   \item{\code{predict()}}{Return the quantile fitted}.
#'   \item{\code{coef()}}{Return coefficients fitted}.
#'   \item{\code{is_fitted()}}{TRUE if fit is already called}.
#'   \item{\code{is_success()}}{TRUE if fit is a success}.
#'   \item{\code{is_unfeasible()}}{TRUE if fit is not feasible}.
#' }
#' @examples
#' ## Data
#' size = 2000
#' c_data = dataset.covariates(size)
#' 
#' loc   = 0.5 + 2 * c_data$X_loc
#' scale = 1 + 2 * c_data$X_scale
#' Y = stats::rnorm( size , mean = loc , sd = scale )
#' 
#' ## Quantile regression
#' ltau  = base::seq( 0.01 , 0.99 , length = 100 )
#' qr = SDFC::QuantileRegression$new(ltau)
#' qr$fit( Y , base::cbind( c_data$X_loc , c_data$X_scale ) )
#' Yq = if( qr$is_success() ) qr$predict() else NULL
#' @export
QuantileRegression = R6::R6Class( "QuantileRegression" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	ltau    = NULL,
	maxit   = NULL,
	tol     = NULL,
	beta    = NULL,
	method  = NULL,
	verbose = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( ltau , method = "Frish-Newton" , verbose = FALSE ) ##{{{
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
	
	
	###############
	## Accessors ##
	###############
	
	coef = function() ##{{{
	{
		return( private$qrmethod$coef() )
	},
	##}}}
	
	
	###########
	## State ##
	###########
	
	is_fitted = function() ##{{{
	{
		return( private$qrmethod$is_fitted() )
	},
	##}}}
	
	is_success = function() ##{{{
	{
		return( private$qrmethod$is_success() )
	},
	##}}}
	
	is_unfeasible = function() ##{{{
	{
		return( private$qrmethod$is_unfeasible() )
	},
	##}}}
	
	#############
	## Methods ##
	#############
	
	fit = function( Y , X ) ##{{{
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
	
	predict = function() ##{{{
	{
		return( private$qrmethod$predict() )
	}
	##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	qrmethod = NULL
	
	
	#############
	## Methods ##
	#############
	
	
	)
)

