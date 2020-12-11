 
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


#' std
#'
#' Compute the standard deviation with covariates and link function
#'
#' @param Y   [vector] Dataset to fit
#' @param c_Y [vector or NULL] Covariate
#' @param m_Y  [vector or NULL] mean already (or not) estimated. If NULL, m = base::mean(Y) is called
#' @param link  [SDFC::UnivariateLink] link function, default is identity
#' @param value  [bool] if TRUE return mean, else return coefficients of the fit
#'
#' @return [vector] Standard deviation or coefficients of regression
#'
#' @examples
#' ## Data
#' size  = 2500
#' t     = base::seq( 0 , 1 , length = size )
#' X0    = t^2
#' X1    = base::cos( 2 * base::pi * t )
#' loc   = 1. + 2 * X0
#' scale = 0.6 + 0.5 * X1
#' Y     = stats::rnorm( n = size , mean = loc , sd = scale )
#'
#' m = SDFC::mean( Y , c_Y = X0 ) ## First fit mean
#' s = SDFC::std( Y , c_Y = X1 , m_Y = m ) ## Now standard deviation
#' 
#' @export
std = function( Y , c_Y = NULL , m_Y = NULL , link = SDFC::ULIdentity$new() , value = TRUE )
{
	var = SDFC::np_var( Y , c_Y , m_Y , link )
	out = base::sqrt( var )
	
	if( !value )
	{
	
		if( is.null(c_Y) )
		{
			return(link$inverse(out))
		}
		else
		{
			YY = link$inverse(out)
			return( as.vector(stats::lm( YY ~ c_Y )$coefficients) )
		}
	}
	
	return(out)
}




