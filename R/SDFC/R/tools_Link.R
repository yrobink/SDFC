
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


## Base class for Link {{{

#' Base class for Link function
#'
#' Base class used to define generic link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{LinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of LinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of LinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of LinkFct at x}
#' }
#' @examples
#'
#' @export
AbstractLink = R6::R6Class( "AbstractLink" ,
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
	},

	eval = function(x)
	{
	},
	
	inverse = function(x)
	{
	},
	
	gradient = function(x)
	{
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## ChainLink {{{

#' ChainLink
#'
#' Chain link function to chain two link functions
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#' @param link0 [LinkFct] First link function to apply
#' @param link1 [LinkFct] Second link function to apply
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(link1,link0)}}{This method is used to create object of this class with \code{ChainLink}}
#'   \item{\code{eval(x)}}{Evaluation of ChainLink at x}
#'   \item{\code{inverse(x)}}{Inverse of ChainLink at x}
#'   \item{\code{gradient(x)}}{Gradient of ChainLink at x}
#' }
#' @examples
#'
#' @export
ChainLink = R6::R6Class( "ChainLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	link0 = NULL,
	link1 = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( link1 , link0 )
	{
		super$initialize()
		self$link0 = link0
		self$link1 = link1
	},
	
	eval = function(x)
	{
		return( self$link1$eval( self$link0$eval(x) ) )
	},
	
	inverse = function(x)
	{
		return( self$link0$inverse( self$link1$inverse(x) ) )
	},
	
	gradient = function(x)
	{
		return( self$link0$gradient(x) * self$link1$gradient( self$link0$eval(x) ) )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## IdLink {{{

#' IdLink
#'
#' Identity link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{IdLink}}
#'   \item{\code{eval(x)}}{Evaluation of IdLink at x}
#'   \item{\code{inverse(x)}}{Inverse of IdLink at x}
#'   \item{\code{gradient(x)}}{Gradient of IdLink at x}
#' }
#' @examples
#'
#' @export
IdLink = R6::R6Class( "IdLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		super$initialize()
	},
	
	eval = function(x)
	{
		return(x)
	},
	
	inverse = function(x)
	{
		return(x)
	},
	
	gradient = function(x)
	{
		return( base::rep( 1. , length(x) ) )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## ExpLink {{{

#' ExpLink
#'
#' Exponential link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{ExpLink}}
#'   \item{\code{eval(x)}}{Evaluation of ExpLink at x}
#'   \item{\code{inverse(x)}}{Inverse of ExpLink at x}
#'   \item{\code{gradient(x)}}{Gradient of ExpLink at x}
#' }
#' @examples
#'
#' @export
ExpLink = R6::R6Class( "ExpLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		super$initialize()
	},
	
	eval = function(x)
	{
		return( base::exp(x) )
	},
	
	inverse = function(x)
	{
		return( base::log(x) )
	},
	
	gradient = function(x)
	{
		return( base::exp(x) )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## InverseLink {{{

#' InverseLink
#'
#' Inverse link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{InverseLink}}
#'   \item{\code{eval(x)}}{Evaluation of InverseLink at x}
#'   \item{\code{inverse(x)}}{Inverse of InverseLink at x}
#'   \item{\code{gradient(x)}}{Gradient of InverseLink at x}
#' }
#' @examples
#'
#' @export
InverseLink = R6::R6Class( "InverseLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		super$initialize()
	},
	
	eval = function(x)
	{
		return( 1. / x )
	},
	
	inverse = function(x)
	{
		return( 1. / x )
	},
	
	gradient = function(x)
	{
		return( - 1. / x^2 )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## LogitLink {{{

#' LogitLink
#'
#' Logit link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param a [float] Lower bound of logit
#' @param b [float] Upper bound of logit
#' @param s [float] Speed of logit between a and b
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(a,b,s)}}{This method is used to create object of this class with \code{LogitLink}}
#'   \item{\code{eval(x)}}{Evaluation of LogitLink at x}
#'   \item{\code{inverse(x)}}{Inverse of LogitLink at x}
#'   \item{\code{gradient(x)}}{Gradient of LogitLink at x}
#' }
#' @examples
#'
#' @export
LogitLink = R6::R6Class( "LogitLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	a = 0,
	b = 1,
	s = 1,
	
	#################
	## Constructor ##
	#################
	
	initialize = function( a = 0 , b = 1 , s = 1 )
	{
		super$initialize()
		self$a = a
		self$b = b
		self$s = s
	},
	
	eval = function(x)
	{
		return( (self$b - self$a) / ( 1. + base::exp(- self$s * x) ) + self$a )
	},
	
	inverse = function(x)
	{
		idx_lo = x < self$a
		idx_up = x > self$b
		x[idx_lo] = self$a + 1e-3
		x[idx_up] = self$b - 1e-3
		return( - base::log( (self$b - self$a) / (x - self$a) - 1 ) / self$s )
	},
	
	gradient = function(x)
	{
		e = base::exp( - self$s * x )
		return( self$s * (self$b - self$a) * e / ( 1 + e )^2 )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}


