
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



#' lmoments
#'
#' Compute the L-Moments
#'
#' @param Y  [vector] Dataset
#' @param c_Y  [vector] Covariate. If NULL stationary L-moments are computed
#' @param order [int] order of moments. a vector with elements between 1 and 4
#' @param lq [vector] Vector of quantile used for quantile regression. Default is seq(0.05,0.95,0.01)
#'
#' @return [lmom] L-Moments
#'
#' @examples
#' ## Data
#' size = 2000
#' c_data = dataset(size)
#' 
#' t       = c_data$t
#' X_loc   = c_data$X_loc
#' X_scale = c_data$X_scale
#' loc   = 0.5 + 2 * X_loc
#' scale =   1 + 2 * X_scale
#' Y = stats::rnorm( size , mean = loc , sd = scale )
#' 
#' c_Y = base::cbind( X_loc , X_scale )
#' 
#' lmom = SDFC::lmoments( Y , c_Y = c_Y )
#' @export
lmoments = function( Y , c_Y = NULL , order = NULL , lq = base::seq( 0.05 , 0.95 , 0.01 ) )
{
	
	lmoments_stationary = function(Y)##{{{
	{
		Ys = Y[order(Y)]
		size = length(Ys)
		lmom = numeric(4)
		
		## Order 2
		C0 = base::choose( 1:size , 1 )
		C1 = base::choose( base::seq( size - 1 , 0 , -1 ) , 1 )
		
		## Order 3
		C2 = base::choose( 1:size , 2 )
		C3 = base::choose( base::seq( size - 1 , 0 , -1 ) , 2 )
		
		## Order 4
		C4 = base::choose( 1:size , 3 )
		C5 = base::choose( base::seq( size - 1 , 0 , -1 ) , 3 )
		
		
		lmom[1] = base::mean(Ys)
		lmom[2] = base::sum( ( C0 - C1 ) * Ys ) / ( 2 * base::choose( size , 2 ) )
		lmom[3] = base::sum( ( C2 - 2 * C0 * C1 + C3 ) * Ys ) / ( 3 * base::choose( size , 3 ) )
		lmom[4] = base::sum( ( C4 - 3 * C2 * C1 + 3 * C0 * C3 - C5 ) * Ys ) / ( 4 * base::choose( size , 4 ) )
		
		return(lmom)
	}
	##}}}
	
	if( is.null(order) )
		order = 1:4
	
	if( is.null(c_Y) )
	{
		lmom = lmoments_stationary(Y)
		return(lmom[order])
	}
	
	if( !is.matrix(c_Y) )
		c_Y = matrix( c_Y , nrow = length(c_Y) , ncol = 1 )
	
	Yq = SDFC::quantile( Y , lq , c_Y )
	lmom = base::t( base::apply( Yq , 1 , lmoments_stationary ) )
	
	return( lmom[,order] )
}



