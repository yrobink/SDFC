 
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


#' var
#'
#' Compute the variance with covariates and link function
#'
#' @param Y   [vector] Dataset to fit
#' @param c_Y [vector or NULL] Covariate
#' @param m_Y [vector or NULL] mean already (or not) estimated. If NULL, m = base::mean(Y) is called
#' @param link  [SDFC::UnivariateLink] link function, default is identity
#' @param value  [bool] if TRUE return variance, else return coefficients of the fit
#' @param ... Arguments of stats::var used only if c_Y is NULL
#'
#' @return [vector] Variance or coefficients of regression
#'
#' @examples
#' ## Data
#' size = 2500
#' t    = base::seq( 0 , 1 , length = size )
#' X0    = t^2
#' X1    = base::cos( 2 * base::pi * t )
#' loc   = 1. + 2 * X0
#' scale = 0.6 + 0.5 * X1
#' Y    = stats::rnorm( n = size , mean = loc , sd = scale )
#'
#' m = SDFC::mean( Y , c_Y = X0 ) ## First fit mean
#' v = SDFC::var( Y , c_Y = X1 , m_Y = m ) ## Now variance
#' 
#' @export
var = function( Y , c_Y = NULL , m_Y = NULL , link = SDFC::ULIdentity$new() , value = TRUE , ... )
{
	out  = NULL
	coef = NULL
	m_Y  = if( is.null(m_Y) ) base::mean(Y) else as.vector(m_Y)
	if( is.null(c_Y) )
	{
		kwargs = list(...)
		kwargs[["x"]] = Y - m_Y
		out  = base::do.call( stats::var , kwargs )
		coef = link$inverse(out)
	}
	else
	{
		Yres = link$inverse( (Y - m_Y)^2 )
		coef = as.vector(stats::lm( Yres ~ c_Y )$coefficients)
		out  = base::abs(link$transform( base::cbind( 1 , c_Y ) %*% coef ))
	}
	
	
	if( value )
		return(out)
	else
		return(coef)
}




