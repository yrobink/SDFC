
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


#' RHS 
#'
#' @description
#' Right Hand Side
#'
#' @details
#' Right part of fitted equation
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
RHS = R6::R6Class( "RHS" ,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	#' @field coef_ [vector] Coefficient fitted
	.coef_ = NULL
	#' @field n_features [integer] Number of features, i.e. length of coef_
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	coef_ = function(value)
	{
		if( missing(value) )
		{
			return(private$.coef_)
		}
		else
		{
			private$.coef_    = value
			values            = self$l_global$transform( private$.coef_ , self$c_global )
			names(values)     = self$lhs$names
			self$lhs$values   = values
			self$lhs$jacobian = self$l_global$jacobian(  private$.coef_ , self$c_global )
		}
	},
	
	n_features = function(value)
	{
		if( missing(value) )
			return(self$l_global$n_features)
	}
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	## Properties
	#' @field lhs [SDFC::LHS] LHS part
	lhs      = NULL,
	#' @field c_global [list] List of covariate
	c_global = NULL,
	#' @field l_global [SDFC::MultivariateLink] Link function
	l_global = NULL,
	#' @field s_global [list] Length of coef per params of lhs
	s_global = NULL,
	
	#' @description
    #' Create a new RHS object.
    #' @param lhs Left Hand Side
	#' @return A new `RHS` object.
	initialize = function(lhs)
	{
		self$lhs = lhs
	},
	
#	"""
#	Here five kinds of arguments can be passed:
#	- c_<param> : covariate of the param,
#	- l_<param> : link function of the param,
#	- f_<param> : fixed values of the LHS
#	- c_global  : list of all covariates, sorted by lhs order
#	- l_global  : global link function generated the LHS
#	If c_global is set, all arguments (except l_global) are ignored
#	"""
	#' @description
	#' Build all parameters of RHS object
	#' @param ... Many parameters
	build = function(...)
	{
		kwargs = list(...)
		
		## If global covariate and link functions are defined, just set it
		if( !is.null(kwargs[["c_global"]]) )
		{
			self$c_global = kwargs[["c_global"]]
			self$l_global = kwargs[["l_global"]]
			return()
		}
		
		## Else loop on lhs to find global parameters
		self$c_global = list()
		self$s_global = list()
		l_global      = list() ## This list will be passed to MLTensor
		
		for( lhs in self$lhs$names )
		{
			## Start with covariate
			if( !is.null(kwargs[[base::paste0("c_",lhs)]]) )
			{
				c = as.matrix(kwargs[[base::paste0("c_",lhs)]])
				if( length(dim(c)) == 1 )
					c = matrix( c , ncol = length(c) )
				self$c_global[[lhs]]  = c
				self$s_global[[lhs]]  = 1 + base::dim(c)[2]
				self$lhs$fixed[[lhs]] = FALSE
			}
			else
			{
				## No covariate, two choices : lhs is 1d or fixed
				self$c_global[[lhs]] = 0
				if( !is.null(kwargs[[base::paste0("f_",lhs)]]) )
				{
					self$s_global[[lhs]]  = 0
					self$lhs$fixed[[lhs]] = TRUE
				}
				else
				{
					self$s_global[[lhs]]  = 1
				}
			}
			
			## Now the link function
			if( !is.null(kwargs[[base::paste0("f_",lhs)]]) )
			{
				l_global[[lhs]] = MLConstant$new( kwargs[[base::paste0("f_",lhs)]] , n_samples = self$lhs$n_samples )
			}
			else
			{
				l = kwargs[[base::paste0("l_",lhs)]]
				if( is.null(l) || "UnivariateLink" %in% class(l) )
				{
					l = MLLinear$new( c_ = self$c_global[[lhs]] , l_ = l , n_samples = self$lhs$n_samples )
				}
				l_global[[lhs]] = l
			}
		}
		self$l_global = MLTensor$new( l_global , self$s_global , n_features = Reduce("+",self$s_global) , n_samples = self$lhs$n_samples )
	}
	
	)
	
	##}}}
	
)

