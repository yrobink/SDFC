
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


#' MultivariateNormal
#'
#' @description
#' Multivariate Normal distribution
#'
#' @details
#' Multivariate Normal generator, used to define a prior for Bayesian fit
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
MultivariateNormal = R6::R6Class( "MultivariateNormal" ,
	
	## Private list ##{{{
	
	private = list(
	
	## Fields ##{{{
	
	#' @field mean [vector] mean of multivariate Normal law
	.mean = NULL,
	#' @field cov [matrix] covariance of multivariate Normal law
	.cov  = NULL,
	#' @field std [matrix] square root of covariance matrix
	.std  = NULL
	
	##}}}
	
	),
	
	##}}}
	
	## Active bindings ##{{{
	
	active = list(
	
	mean = function(value)
	{
		if( missing(value) )
		{
			return(private$.mean)
		}
		else
		{
			private$.mean = value
		}
	},
	
	cov = function(value)
	{
		if( missing(value) )
		{
			return(private$.cov)
		}
		else
		{
			private$.cov = value
			svd          = base::svd(value)
			private$.std = svd$u %*% base::diag(base::sqrt(svd$d)) %*% base::t(svd$u)
		}
	},
	
	std = function(value)
	{
		if( missing(value) )
		{
			return(private$.std)
		}
		else
		{
			private$.std = std
			svd          = base::svd(value)
			private$.cov = svd$u %*% base::diag(svd$d^2) %*% base::t(svd$u)
		}
	}
	
	),
	
	##}}}
	
	## Public list ##{{{
	
	public = list(
	
	## init ##{{{
	
	#' @description
    #' Create a new MultivariateNormal object.
    #' @param mean mean
    #' @param cov covariance matrix
	#' @return A new `MultivariateNormal` object.
	initialize = function( mean , cov )
	{
		self$mean = mean
		self$cov  = cov
	},
	
	##}}}
	
	## rvs ##{{{
	
	#' @description
    #' Random values generator
    #' @param n Number of samples to drawn
	#' @return X matrix of samples
	rvs = function(n)
	{
		n_dim = length(self$mean)
		X = base::matrix( stats::rnorm(n * n_dim) , nrow = n, ncol = n_dim )
		X = X %*% self$std
		X = base::t(base::t(X) + self$mean)
		
		return(X)
	}
	##}}}
	
	##}}}
	
	)
)
