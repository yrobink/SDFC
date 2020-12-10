
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


#' dataset
#'
#' Function used to generate examples of covariates
#'
#' @usage dataset.covariates(size)
#'
#' @param size   [int]    Numbers of values generated
#'
#' @return [list] list containing three covariates and a time axis
#'
#' @examples
#' ## Data
#' c_data = dataset(2000)
#' c_data$t       ## Time axis
#' c_data$t       ## Covariate for loc
#' c_data$t       ## Covariate for scale
#' c_data$t       ## Covariate for shape
#' @export
dataset = function( size )
{
	t       = base::seq( 0 , 1 , length = size )
	X_loc   = t^2 + base::cos( 2 * base::pi * t ) * 0.2
	X_scale = 2 * t^2 - 2 * t + 1
	X_shape = 2 / ( 1 + base::exp( - 8 * ( t - 0.5 ) ) ) - 1
	return( list( t = t , X_loc = X_loc , X_scale = X_scale , X_shape = X_shape ) )
}



