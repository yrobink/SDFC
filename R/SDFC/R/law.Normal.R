
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


#' Normal
#'
#' @description
#' Normal distribution
#'
#' @details
#' Class of Normal distribution
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
Normal = R6::R6Class( "Normal" ,
	
	inherit = AbstractLaw,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	loc = function(value)
	{
		return(private$.rhs$lhs$values$loc)
	},
	
	scale = function(value)
	{
		return(private$.rhs$lhs$values$scale)
	}
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	
	## Public fields ##{{{
	
	#' @field loc [numeric] Location parameter
	#' @field scale [numeric] Scale parameter
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new Normal object.
    #' @param method [string] method used to fit
	#' @return A new `Normal` object.
	initialize = function( method = "mle" )
	{
		super$initialize( base::c("loc","scale") , method )
	}
	
	##}}}
	
	)
	
	##}}}
	
)




