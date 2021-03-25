
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

## ginv {{{

#' ginv (Moore-Penrose generalized inverse of a matrix)
#'
#' Compute the Moore-Penrose generalized inverse of a matrix, code from MASS package
#'
#' @usage ginv(X,tol)
#'
#' @param X [matrix] a matrix
#' @param tol [double] Numerical tolerance
#'
#' @return Xinv [matrix] Generalized inverse
#'
#' @examples
#'
#' X = matrix( base::c(0,0,1,1) , nrow = 2 ) ## Not invertible (det==0)
#' Xinv = ginv(X) ## But generalized inverse exist
#' 
#' @export
ginv = function( X , tol = base::sqrt(.Machine$double.eps) )
{
	if( length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)) )
		stop("'X' must be a numeric or complex matrix")
	
	if( !is.matrix(X) )
		X = as.matrix(X)
	
	Xsvd = svd(X)
	
	if( is.complex(X) )
		Xsvd$u = base::Conj(Xsvd$u)
	
	Positive = ( Xsvd$d > base::max( tol * Xsvd$d[1L] , 0 ) )
	
	out = NULL
	if( base::all(Positive) )
	{
		out = Xsvd$v %*% ( 1. / Xsvd$d * base::t(Xsvd$u) )
	}
	else if( !base::any(Positive) )
	{
		out = array( 0 , dim(X)[2L:1L] )
	}
	else
	{
		out = Xsvd$v[, Positive, drop = FALSE] %*% ( ( 1. / Xsvd$d[Positive] ) * base::t( Xsvd$u[, Positive, drop = FALSE] ) )
	}
	return(out)
}
##}}}


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
	#' @field cov [matrix] covariance matrix of multivariate Normal law
	.cov  = NULL,
	#' @field std [matrix] square root of covariance matrix
	.std  = NULL,
	#' @field icov [matrix] inverse of covariance matrix
	.icov = NULL,
	#' @field det [double] determinant of covariance matrix
	.det  = NULL,
	#' @field rank [integer] rank of covariance matrix
	.rank = NULL,
	
	##}}}
	
	update_icov_det = function() ##{{{
	{
		private$.icov = ginv(private$.cov)
		svd           = base::svd(private$.cov)
		d             = svd$d[svd$d > 0]
		private$.det  = base::prod(d)
		private$.rank = length(d)
	}
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
			private$.cov  = value
			svd           = base::svd(value)
			private$.std  = svd$u %*% base::diag(base::sqrt(svd$d)) %*% base::t(svd$u)
			private$update_icov_det()
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
			private$.std = value
			svd          = base::svd(value)
			private$.cov = svd$u %*% base::diag(svd$d^2) %*% base::t(svd$u)
			private$update_icov_det()
		}
	},
	
	icov = function(value)
	{
		if( missing(value) )
		{
			return(private$.icov)
		}
	},
	
	det = function(value)
	{
		if( missing(value) )
		{
			return(private$.det)
		}
	},
	
	rank = function(value)
	{
		if( missing(value) )
		{
			return(private$.rank)
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
	},
	##}}}
	
	## pdf {{{
	
	#' @description
    #' Probability density function, this function works with degenerate distribution
    #' @param X Compute pdf at values of matrix X
	#' @return pdf The PDF
	pdf = function(X)
	{
		X  = matrix( X , ncol = length(self$mean) )
		Xm = base::t( base::t(X) - self$mean )
		
		e   = base::apply( Xm * base::t(self$icov %*% t(Xm)) , 1 , base::sum )
		val = base::exp( - e / 2 ) / base::sqrt( self$det * (2*pi)^self$rank )
		
		return(val)
	},
	##}}}
	
	## logpdf ##{{{
	
	#' @description
    #' Log of probability density function
    #' @param X Compute log of pdf at values of matrix X
	#' @return logpdf Log of the PDF
	logpdf = function(X)
	{
		return(base::log(self$pdf(X)))
	},
	##}}}
	
	## fit ##{{{
	
	#' @description
    #' Fit the parameters from a dataset
    #' @param X Data to fit
	#' @return NULL
	fit = function(X)
	{
		self$mean = base::apply( X , 2 , base::mean )
		self$cov  = stats::cov(X)
		invisible(NULL)
	}
	##}}}
	
	##}}}
	
	)
)
