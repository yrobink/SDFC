
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


## UnivariateLink ##{{{

#' UnivariateLink 
#'
#' @description
#' Base class of UnivariateLink
#'
#' @details
#' Just to be inherit
#'
#' @examples
#' ##
#'
#' @export
UnivariateLink = R6::R6Class( "UnivariateLink" ,
	
	public = list(
	
	#' @description
    #' Create a new UnivariateLink object.
	#' @return A new `UnivariateLink` object.
	initialize = function()
	{
	}
	
	)
	
)
##}}}

## ULIdentity ##{{{

#' ULIdentity 
#'
#' @description
#' Identity link function
#'
#' @details
#' Identity link function
#'
#' @examples
#'
#' idlink = SDFC::ULIdentity$new()
#'
#' @export
ULIdentity = R6::R6Class( "ULIdentity" ,
	
	inherit = UnivariateLink,
	
	public = list(
	
	#' @description
    #' Create a new ULIdentity object.
	#' @return A new `ULIdentity` object.
	initialize = function(...)
	{
		super$initialize()
	},
	
	#' @description
    #' Transform method
    #' @param x [vector] Value to transform
    #' @return Value transformed
	transform = function(x)
	{
		return(x)
	},
	
	#' @description
    #' Inverse method
    #' @param x [vector] Value to inverse
    #' @return Value inversed
	inverse = function(x)
	{
		return(x)
	},
	
	#' @description
    #' Jacobian method
    #' @param x [vector] Value to evaluate the jacobian
    #' @return Jacobian at x
	jacobian = function(x)
	{
		return(x)
	}
	
	)
)

##}}}


