 
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


#' median
#'
#' Compute the median with covariates (in fact just call quantile)
#'
#' @param Y  [vector] Dataset to fit
#' @param c_Y  [vector or NA] Covariate
#' @param value  [bool] if TRUE return variance, else return coefficients of the fit
#' @param ... Arguments of stats::median used only if c_Y is NULL
#'
#' @return [vector] Median or coefficients of regression
#'
#' @examples
#' ## Data
#' size = 2500
#' t    = base::seq( 0 , 1 , length = size )
#' X0    = t^2
#' loc   = 1. + 2 * X0
#' Y    = stats::rnorm( n = size , mean = loc , sd = 0.1 )
#'
#' med = SDFC::median( Y , c_Y = X0 )
#' 
#' @export
median = function( Y , c_Y = NA , value = TRUE , ... )
{
	if( base::any(is.na(c_Y)) )
	{
		kwargs = list(...)
		kwargs[["x"]] = Y
		return(base::do.call( stats::median , kwargs ))
	}
	return( SDFC::quantile( Y , base::c(0.5) , c_Y , value ) )
}
