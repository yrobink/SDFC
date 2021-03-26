
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


## MultivariateLink ##{{{

#' MultivariateLink 
#'
#' @description
#' Base class of MultivariateLink
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
MultivariateLink = R6::R6Class( "MultivariateLink" ,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	#' @field special_fit_allowed [bool] If special fit (as moments) is allowed
	.special_fit_allowed = NULL,
	#' @field n_features [integer] Numbers of parameters to fit
	.n_features = NULL,
	#' @field n_samples [integer] Numbers of samples from the dataset to fit
	.n_samples  = NULL
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	n_samples = function(value)
	{
		if(missing(value))
		{
			return(private$.n_samples)
		}
	},
	
	n_features = function(value)
	{
		if(missing(value))
		{
			return(private$.n_features)
		}
	},
	
	special_fit_allowed = function(value)
	{
		if(missing(value))
		{
			return(private$.special_fit_allowed)
		}
	}
	
	),
	
	##}}}
	
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	#' @description
    #' Create a new MultivariateLink object.
    #' @param ... Some arguments
	#' @return A new `MultivariateLink` object.
	initialize = function(...)
	{
		kwargs = list(...)
		private$.special_fit_allowed = FALSE
		private$.n_features = kwargs[["n_features"]]
		private$.n_samples  = kwargs[["n_samples"]]
	}
	
	)
	
	##}}}
	
)
##}}}

## MLConstant ##{{{

#' MLConstant 
#'
#' @description
#' Constant link function
#'
#' @details
#' Modelized a link function with fixed values
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
MLConstant = R6::R6Class( "MLConstant" ,
	
	inherit = MultivariateLink,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	.value = NULL
	
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
	
	#' @description
    #' Create a new MLConstant object.
    #' @param value [vector] Fixed value of the parameter
    #' @param ... Some arguments
	#' @return A new `MLConstant` object.
	initialize = function( value , ... )
	{
		kwargs = list(...)
		kwargs[["n_features"]] = 0
		base::do.call( super$initialize , kwargs )
		private$.value = as.vector(value)
		if( length(private$.value) == 1 )
		{
			private$.value = base::rep( value , self$n_samples )
		}
	},
	
	#' @description
    #' Transform function
    #' @param ... Any arguments, not used
	#' @return The value of fixed parameter
	transform = function(...)
	{
		return(private$.value)
	}
	
	
	)
	
	##}}}
	
)
##}}}

## MLLinear ##{{{

#' MLLinear 
#'
#' @description
#' Linear link function
#'
#' @details
#' Modelized a linear form link function
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
MLLinear = R6::R6Class( "MLLinear" ,
	
	inherit = MultivariateLink,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	.l = NULL,
	.c = NULL,
	#' @field l [UnivariateLink] Univariate link
	
	linear_transform = function( coef , X )
	{
		return(self$design_ %*% coef)
	}
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	l = function(value) ##{{{
	{
		if( missing(value) )
			return(private$.l)
	}
	##}}}
	
	),
	
	##}}}
	
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	#' @field design_ [matrix] Design matrix
	design_ = NULL,
	
	#' @description
    #' Create a new MLLinear object.
    #' @param ... Some arguments
	#' @return A new `MLLinear` object.
	initialize = function( ... )
	{
		kwargs = list(...)
		private$.l = kwargs[["l_"]]
		private$.c = kwargs[["c_"]]
		base::do.call( super$initialize , kwargs )
		
		self$design_ = matrix( 1 , nrow = self$n_samples , ncol = 1 )
		if( !base::any(is.na(private$.c)) )
			self$design_ = base::cbind( self$design_ , private$.c )
		private$.n_features = base::ncol(self$design_)
		if( is.null(private$.l) )
			private$.l = ULIdentity$new()
	},
	
	#' @description
    #' Transform function
    #' @param coef coefficients to fit
    #' @param X co-variates
	#' @return The value
	transform = function( coef , X )
	{
		out = private$.l$transform( private$linear_transform(coef,X) )
		return(out)
	},
	
	#' @description
    #' Jacobian of transform function
    #' @param coef coefficients to fit
    #' @param X co-variates
	#' @return The value
	jacobian = function( coef , X )
	{
		out = private$.l$jacobian( matrix( private$linear_transform(coef,X) , ncol = 1 ) )
		out = as.vector(out) * self$design_
		return(out)
	}
	
	)
	
	##}}}
	
)
##}}}

## MLTensor ##{{{

#' MLTensor 
#'
#' @description
#' Tensor link function
#'
#' @details
#' Link function used to build the product of univariate link function
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
MLTensor = R6::R6Class( "MLTensor" ,
	
	inherit = MultivariateLink,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	.l_p = NULL,
	.s_p = NULL,
	.special_fit_allowed = FALSE
	
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
	
	## Init ##{{{
	
	#' @description
    #' Create a new MLTensor object.
    #' @param l_p List of Link function
    #' @param s_p List of numbers of features per marginal
    #' @param ... Some arguments
	#' @return A new `MLTensor` object.
	initialize = function( l_p , s_p , ... )
	{
		kwargs = list(...)
		kwargs[["n_features"]] = Reduce("+",s_p)
		base::do.call( super$initialize , kwargs )
		private$.l_p = l_p
		private$.s_p = s_p
		private$.special_fit_allowed = 1:length(private$.l_p)
		sfa = base::rep(FALSE,length(private$.l_p))
		for( i in 1:length(private$.l_p) )
		{
			sfa[i] = ("MLLinear" %in% class(private$.l_p[[i]])) || ("MLConstant" %in% class(private$.l_p[[i]]))
		}
		private$.special_fit_allowed = base::all(sfa)
	},
	##}}}
	
	## Transform ##{{{
	
	#' @description
    #' Transform function
    #' @param coef coefficients to fit
    #' @param X co-variates
	#' @return The value
	transform = function( coef , X )
	{
		out = list()
		ib = 1
		ie = 1
		for( i in 1:length(X) )
		{
			ie = ie + private$.s_p[[i]] - 1
			out[[i]] = private$.l_p[[i]]$transform( coef[ib:ie] , X[[i]] )
			ib = ib + private$.s_p[[i]]
			ie = ie + 1
		}
		return(out)
	},
	##}}}
	
	## Jacobian ##{{{
	
	#' @description
    #' Jacobian of transform function
    #' @param coef coefficients to fit
    #' @param X co-variates
	#' @return The value
	jacobian = function( coef , X )
	{
		jac = array( 0 , dim = base::c( base::sum(as.numeric(private$.s_p) > 0) , self$n_samples , self$n_features ) )
		ib = 1
		ie = 1
		ii = 1
		for( i in 1:length(X) )
		{
			if( private$.s_p[[i]] > 0 )
			{
				ie = ie + private$.s_p[[i]] - 1
				jac[ii,,ib:ie] = private$.l_p[[i]]$jacobian( coef[ib:ie] , X[[i]] )
				ib = ib + private$.s_p[[i]]
				ie = ie + 1
				ii = ii + 1
			}
		}
		
		return(jac)
	},
	##}}}
	
	## linkc {{{
	
	#' @description
    #' Link function component
    #' @param name name of component
	#' @return The link
	linkc = function(name)
	{
		return(private$.l_p[[name]])
	},
	##}}}
	
	## sizec {{{
	
	#' @description
    #' Size of each component
    #' @param name name of component
	#' @return The size
	sizec = function(name)
	{
		return(private$.s_p[[name]])
	}
	##}}}
	
	)
	
	##}}}
	
)
##}}}

