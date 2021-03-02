
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


#' LHS 
#'
#' @description
#' Left Hand Side
#'
#' @details
#' Left part of fitted equation
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
LHS = R6::R6Class( "LHS" ,
	
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
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	## Properties
	#' @field names [list(string)] List of names
	names      = NULL,
	#' @field n_lhs [integer] Number of LHS values
	n_lhs      = NULL,
	#' @field n_samples [integer] Number of samples
	n_samples  = NULL,
	#' @field jacobian [matrix] Jacobian matrix
	jacobian   = NULL,
	#' @field values [vector] LHS values
	values = NULL,
	#' @field fixed [integer] Numbers of parameters to fit
	fixed  = NULL,
	
	#' @description
    #' Create a new LHS object.
    #' @param names     names of LHS parameters
    #' @param n_samples Number of samples
	#' @return A new `LHS` object.
	initialize = function(names,n_samples)
	{
		self$names      = names
		self$n_lhs      = length(self$names)
		self$n_samples  = n_samples
		self$values = list()
		self$fixed  = list()
	}
	
	)
	
	##}}}
	
)


