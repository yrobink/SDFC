
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
#' @importFrom R6 R6Class
#' @importFrom methods new
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
#' Identity link function:
#' f(x)    = x
#' f^-1(x) = x
#' f'(x)   = 1
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#'
#' idlink = SDFC::ULIdentity$new()
#'
#' @export
ULIdentity = R6::R6Class( "ULIdentity" ,
	
	inherit = UnivariateLink,
	
	public = list(
	
	## initialize ##{{{
	
	#' @description
    #' Create a new ULIdentity object.
	#' @return A new `ULIdentity` object.
	initialize = function()
	{
		super$initialize()
	},
	
	##}}}
	
	## transform ##{{{
	
	#' @description
    #' Transform method
    #' @param x [vector] Value to transform
    #' @return Value transformed
	transform = function(x)
	{
		return(x)
	},
	##}}}
	
	## inverse ##{{{
	
	#' @description
    #' Inverse method
    #' @param x [vector] Value to inverse
    #' @return Value inversed
	inverse = function(x)
	{
		return(x)
	},
	##}}}
	
	## jacobian ##{{{
	
	#' @description
    #' Jacobian method
    #' @param x [vector] Value to evaluate the jacobian
    #' @return Jacobian at x
	jacobian = function(x)
	{
		return(x)
	}
	##}}}
	
	)
)

##}}}

## ULExponential ##{{{

#' ULExponential 
#'
#' @description
#' Exponential link function
#'
#' @details
#' Exponential link function, usefull to constrain a parameter in [b;+inf[ if
#' s > 0 or ]-inf;b] if s < 0.
#' f(x)    = base::exp(s*x) + b
#' f^-1(x) = base::log(x-b) / s
#' f'(x)   = s * base::exp(s*x)
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#'
#' explink = SDFC::ULExponential$new()
#'
#' @export
ULExponential = R6::R6Class( "ULExponential" ,
	
	inherit = UnivariateLink,
	
	public = list(
	
	## Public attributes ## {{{
	#' @field b [float] Lower / Upper bound of exponential link
	b = NULL,
	#' @field s [float] Speed of exponential link
	s = NULL,
	
	##}}}
	
	## initialize ##{{{
	
	#' @description
    #' Create a new ULExponential object.
	#' @param b [float] Lower / Upper bound of exponential link
	#' @param s [float] Speed of exponential link
	#' @return A new `ULExponential` object.
	initialize = function( b = 0 , s = 1 )
	{
		super$initialize()
		self$b = b
		self$s = s
	},
	
	##}}}
	
	## transform ##{{{
	
	#' @description
    #' Transform method
    #' @param x [vector] Value to transform
    #' @return Value transformed
	transform = function(x)
	{
		return( base::exp( self$s * x ) + self$b )
	},
	##}}}
	
	## inverse ##{{{
	
	#' @description
    #' Inverse method
    #' @param x [vector] Value to inverse
    #' @return Value inversed
	inverse = function(x)
	{
		return( base::log( x - self$b ) / self$s )
	},
	##}}}
	
	## jacobian ##{{{
	
	#' @description
    #' Jacobian method
    #' @param x [vector] Value to evaluate the jacobian
    #' @return Jacobian at x
	jacobian = function(x)
	{
		return(self$s * base::exp(self$s * x))
	}
	##}}}
	
	)
)

##}}}

