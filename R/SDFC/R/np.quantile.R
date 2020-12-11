 
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


#' quantile
#'
#' Compute the quantile with covariates
#'
#' @param Y  [vector] Dataset to fit
#' @param ltau [vector] vector of quantiles to fit
#' @param c_Y  [vector or NULL] Covariate
#' @param value  [bool] if TRUE return variance, else return coefficients of the fit
#' @param ... Arguments of stats::quantile used only if c_Y is NULL
#'
#' @return [vector] Quantile or coefficients of regression
#'
#' @examples
#' ## Data
#' size = 2500
#' t    = base::seq( 0 , 1 , length = size )
#' X0    = t^2
#' loc   = 1. + 2 * X0
#' Y    = stats::rnorm( n = size , mean = loc , sd = 0.1 )
#'
#' q = SDFC::quantile( Y , ltau = base::c(0.25,0.5,0.75) , c_Y = X0 )
#' 
#' @export
quantile = function( Y , ltau , c_Y = NULL , value = TRUE , ... )
{
	out  = NULL
	coef = NULL
	if( is.null(c_Y) )
	{
		kwargs = list(...)
		kwargs[["x"]] = Y
		kwargs[["probs"]] = ltau
		out  = as.vector( base::do.call( stats::quantile , kwargs ) )
		coef = out
	}
	else
	{
		reg = QuantileRegression$new( ltau )
		reg$fit( Y , c_Y )
		out = reg$predict()
		coef = reg$coef_
	}
	
	if( value )
		return(out)
	else
		return(coef)
}
