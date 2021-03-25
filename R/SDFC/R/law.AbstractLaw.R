
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



#' AbstractLaw 
#'
#' @description
#' Base class of parametric laws
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
AbstractLaw = R6::R6Class( "AbstractLaw" ,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	## Private fields ##{{{
	
	#' @field method [string] Method used to fit the distribution
	.method = NULL,
	#' @field rhs [SDFC::RHS] Right Hand Side
	.rhs = NULL,
	#' @field coef_ [numeric] Coefficients fitted
	.Y = NULL,
	#' @field lhs [SDFC::RHS] Left Hand Side
	#' @field cov [matrix] Covariance matrix of coefs
	
	##}}}
	
	random_valid_point = function() ##{{{
	{
		coef_ = self$coef_
		scale = 0.1
		
		p_coef = coef_
		n_it = 1
		
		while( !(is.finite(private$negloglikelihood(p_coef)) ) )# && base::all(is.finite(private$gradient_nlll(p_coef)))) )
		{
			if( n_it %% 100 == 0 ) scale = scale * 2
			
			p_coef = coef_ + stats::rnorm( n = length(coef_) , mean = 0 , sd = scale )
			n_it = n_it + 1
		}
		
		self$coef_ = p_coef
	},
	##}}}
	
	fit_MLE = function(...) ##{{{
	{
		kwargs = list()
		private$init_MLE()
		coef_ = self$coef_
		
		is.success = FALSE
		n_test     = 0
		max_test   = kwargs[["mle_n_restart"]]
		if( is.null(max_test) ) max_test = 2
		
		while( !is.success && n_test < max_test )
		{
			private$random_valid_point()
			self$info_[["mle_optim_result"]] = stats::optim( self$coef_ , private$negloglikelihood , gr = private$gradient_nlll , method = "BFGS" , hessian = TRUE )
			self$coef_ = self$info_[["mle_optim_result"]]$par
			is.success = self$info_[["mle_optim_result"]]$convergence == 0
			n_test = n_test + 1
		}
		self$info_[["cov"]] = base::solve(self$info_[["mle_optim_result"]]$hessian)
		
		if( !is.success )
			self$coef_ = coef_
		
	},
	##}}}
	
	fit_Bayesian = function(...)##{{{
	{
		kwargs = list(...)
		
		## Find numbers of features
		n_features = private$.rhs$n_features
		
		## Define prior
		prior = kwargs[["prior"]]
		if( is.null(prior) )
		{
			m = numeric(n_features)
			S = base::diag(10 + numeric(n_features))
			prior = MultivariateNormal$new( m , S )
		}
		
		## Define transition
		transition = kwargs[["transition"]]
		if( is.null(transition))
		{
			transition = function(x) { return( x + stats::rnorm(n_features,0,0.1) ) }
		}
		
		## Define numbers of iterations of MCMC algorithm
		n_mcmc_drawn = kwargs[["n_mcmc_drawn"]]
		if( is.null(n_mcmc_drawn) )
		{
			n_mcmc_drawn = 10000
		}
		
		## Init values
		init = kwargs[["mcmc_init"]]
		if( is.null(init) )
		{
			init = prior$rvs()
		}
		
		## MCMC algorithm
		##===============
		draw   = matrix( 0 , nrow = n_mcmc_drawn , ncol = n_features )
		accept = numeric( n_mcmc_drawn )
		
		draw[1,]      = init
		lll_current   = -private$negloglikelihood(draw[1,])
		prior_current = base::sum( prior$logpdf(draw[1,]) )
		p_current     = prior_current + lll_current
		
		for( i in 2:n_mcmc_drawn )
		{
			draw[i,] = transition(draw[i-1,])
			
			## Likelihood and probability of new points
			lll_next   = - private$negloglikelihood(draw[i,])
			prior_next = base::sum( prior$logpdf(draw[i,]) )
			p_next     = prior_next + lll_next
			
			## Accept or not ?
			p_accept = base::exp( p_next - p_current )
			if( stats::runif(1) < p_accept )
			{
				lll_current   = lll_next
				prior_current = prior_next
				p_current     = p_next
				accept[i]     = TRUE
			}
			else
			{
				draw[i,]  = draw[i-1,]
				accept[i] = FALSE
			}
		}
		
		self$coef_ = base::apply( draw[as.integer(n_mcmc_drawn/2):n_mcmc_drawn,] , 2 , base::mean )
		
		## Update information
		self$info_$draw         = draw
		self$info_$accept       = accept
		self$info_$n_mcmc_drawn = n_mcmc_drawn
		self$info_$rate_accept  = base::sum(accept) / n_mcmc_drawn
		self$info_$cov          = stats::cov(draw)
		
	}
	##}}}
	
	),
	##}}}
	
	## Active list
	##============
	##{{{
	
	active = list(
	
	method = function(value)
	{
		if(missing(value))
		{
			return(private$.method)
		}
		else
		{
			private$.method = base::tolower(value)
		}
	},
	
	lhs = function(value)
	{
		return(private$.rhs$lhs)
	},
	
	rhs = function(value)
	{
		return(private$.rhs)
	},
	
	coef_ = function(value)
	{
		if(missing(value))
		{
			return(private$.rhs$coef_)
		}
		else
		{
			private$.rhs$coef_ = value
		}
	},
	
	cov = function(value)
	{
		if(missing(value))
			return(self$info_$cov)
	}
	
	),
	
	##}}}
	
	## Public list
	##============
	##{{{
	
	public = list(
	
	
	## Public fields ##{{{
	
	#' @field info_ [list] List of information about the fit
	info_ = NULL,
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new AbstractLaw object.
    #' @param names [vector(string)] names of parameters
    #' @param method [string] method used to fit
	#' @return A new `AbstractLaw` object.
	initialize = function( names , method )
	{
		self$method  = method
		lhs = LHS$new(names,0)
		private$.rhs = RHS$new(lhs)
		self$info_   = list()
	},
	
	##}}}
	
	## Fit ##{{{
	
	#' @description
    #' Fit the law.
    #' @param Y [numeric] Data to fit
    #' @param ... Many arguments
	fit = function( Y , ... ) 
	{
		kwargs = list(...)
		
		## Add Y
		private$.Y = matrix( Y , ncol = 1 )
		
		## Init LHS/RHS
		private$.rhs$lhs$n_samples = length(Y)
		base::do.call( private$.rhs$build , kwargs )
		
		self$coef_ = numeric(private$.rhs$n_features)
		
		## Now fit
		if( !(self$method %in% base::c("mle","bayesian")) && private$.rhs$l_global$special_fit_allowed )
		{
			private$special_fit()
		}
		else if( self$method == "mle" )
		{
			base::do.call( private$fit_MLE , kwargs )
		}
		else
		{
			base::do.call( private$fit_Bayesian , kwargs )
		}
		
	}
	##}}}
	
	## Fit bootstrap ##{{{
	
	##}}}
	
	)
	
	##}}}
	
)


