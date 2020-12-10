
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


#' np_mean
#'
#' Compute the mean with covariates and link function
#'
#' @param Y   [vector] Dataset to fit
#' @param c_Y [vector or NULL] Covariate
#' @param link  [SDFC::LinkFct] link function, default is identity
#' @param value  [bool] if TRUE return mean, else return coefficients of the fit
#'
#' @return [vector] Mean or coefficients of regression
#'
#' @examples
#' ## Data
#' size = 2500
#' t    = base::seq( 0 , 1 , length = size )
#' X0    = t^2
#' loc   = 1. + 2 * X0
#' Y    = stats::rnorm( n = size , mean = loc , sd = 0.1 )
#'
#' m = np_mean( Y , c_Y = X0 )
#' 
#' @export
np_mean = function( Y , c_Y = NULL , link = SDFC::IdLink$new() , value = TRUE )
{
	out  = NULL
	coef = NULL
	if( is.null(c_Y) )
	{
		out  = base::mean(Y)
		coef = link$inverse(out)
	}
	else
	{
		YY   = link$inverse(Y)
		coef = as.vector(stats::lm( YY ~ c_Y )$coefficients)
		out  = link$eval( base::cbind( 1 , c_Y )  %*% coef )
	}
	
	if( value )
		return(out)
	else
		return(coef)
}


