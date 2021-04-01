
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


#' Gamma
#'
#' @description
#' Gamma distribution
#'
#' @details
#' Class to fit a Gamma law with covariates, available methods are:
#' 
#' moments  : use empirical estimator of mean and standard deviation to find
#'            loc and scale, possibly with least square regression if
#'            covariates are given
#'
#' bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of
#'            n_mcmc_iteration sample draw from the posterior P(coef_ | Y)
#'
#' mle      : Maximum likelihood estimation
#'
#' 
#' Parameters:
#'
#' scale : scale parameter
#'
#' shape : shape parameter
#'
#' See AbstractLaw for details and examples
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
Gamma = R6::R6Class( "Gamma" ,
	
	inherit = AbstractLaw,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	fit_moments = function() ##{{{
	{
		coef_ = numeric(length(self$coef_))
		
		n_samples = private$.rhs$lhs$n_samples
		
		if( !private$.rhs$lhs$fixed[["scale"]] && !private$.rhs$lhs$fixed[["shape"]] )
		{
			mX = matrix( 1 , nrow = n_samples , ncol = 1 )
			vX = matrix( 1 , nrow = n_samples , ncol = 1 )
			design_scale = private$.rhs$l_global$linkc("scale")$design_
			design_shape = private$.rhs$l_global$linkc("shape")$design_
			if( ncol(design_scale) > 1 && ncol(design_shape) > 1 )
			{
				for( i in 2:ncol(design_scale) )
				{
					for( j in 1:ncol(design_shape) )
					{
						mX = base::cbind( mX , design_scale[,i]   * design_shape[,j]  )
						vX = base::cbind( vX , design_scale[,i]^2 * design_shape[,j]  )
					}
				}
				mX = mX[,2:ncol(mX)]
				vX = vX[,2:ncol(vX)]
			}
			else if( ncol(design_scale) > 1 )
			{
				mX = design_scale[,2:ncol(design_scale)]
				vX = NA
			}
			else if (ncol(design_shape) > 1)
			{
				mX = NA
				vX = design_shape[,2:ncol(design_shape)]^2
			}
			else
			{
				mX = NA
				vX = NA
			}
			
			m = SDFC::mean( private$.Y , mX )
			v = SDFC::var(  private$.Y , vX )
			if( length(m) == 1 ) m = as.vector(m) + numeric(n_samples)
			if( length(v) == 1 ) v = as.vector(v) + numeric(n_samples)
			
			idx   = ( base::abs(m) < 1e-8 ) | ( v < 1e-8 )
			cidx  = !idx
			scale = numeric(length(m))
			shape = numeric(length(m))
			scale[cidx] = base::abs( v[cidx] / m[cidx] )
			shape[cidx] = base::abs( m[cidx]^2 / v[cidx] )
			
			if( base::any(idx) )
			{
				scale[idx] = base::min(scale[cidx])
				shape[idx] = base::min(shape[cidx])
			}
			
			dsc = if( ncol(design_scale) > 1 ) design_scale[,2:ncol(design_scale)] else NA
			dsh = if( ncol(design_shape) > 1 ) design_shape[,2:ncol(design_shape)] else NA
			coefsc = SDFC::mean( scale , dsc , value = FALSE , link = private$.rhs$l_global$linkc("scale")$l )
			coefsh = SDFC::mean( shape , dsh , value = FALSE , link = private$.rhs$l_global$linkc("shape")$l )
			coef_  = base::c(coefsc,coefsh)
		}
		else if( private$.rhs$lhs$fixed[["scale"]] )
		{
			design_shape = private$.rhs$l_global$linkc("shape")$design_
			
			dsh = if( ncol(design_shape) > 1 ) design_shape[,2:ncol(design_shape)] else NA
			m = SDFC::mean( private$.Y  , dsh * self$scale   )
			v = SDFC::var(  private$.Y  , dsh * self$scale^2 )
			
			shape = base::abs( m^2 / v )
			coef_ = SDFC::mean( shape , dsh , value = FALSE , link = private$.rhs$l_global$linkc("shape")$l )
		}
		else if( private$.rhs$lhs$fixed[["shape"]] )
		{
			design_scale = private$.rhs$l_global$linkc("scale")$design_
			
			dsc = if( ncol(design_scale) > 1 ) design_scale[,2:ncol(design_scale)] else NA
			m = SDFC::mean( private$.Y  , dsc   * self$shape )
			v = SDFC::var(  private$.Y  , dsc^2 * self$shape )
			
			scale = base::abs( v / m )
			coef_ = SDFC::mean( scale , dsc , value = FALSE , link = private$.rhs$l_global$linkc("scale")$l )
		}
		
		self$coef_ = coef_
	},
	##}}}
	
	special_fit = function() ##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
	},
	##}}}
	
	init_MLE = function()##{{{
	{
		if( private$.rhs$l_global$special_fit_allowed )
		{
			private$fit_moments()
		}
		else
		{
			self$coef_ = private$.rhs$l_global$valid_point(self)
		}
	},
	##}}}
	
	negloglikelihood = function(coef) ##{{{
	{
		self$coef_ = coef
		if( !base::all(self$scale > 0) || !base::all(self$shape > 0) || !base::all(private$.Y > 0) )
			return(Inf)
		
		return( base::sum( private$.Y / self$scale + base::lgamma(self$shape) + self$shape * base::log(self$scale) - (self$shape-1) * base::log(private$.Y) ) )
	},
	##}}}
	
	gradient_nlll = function(coef) ##{{{
	{
		self$coef_ = coef
		
		if( !base::all(self$scale > 0) || !base::all(self$shape > 0) || !base::all(private$.Y > 0) )
		{
			return( numeric(length(self$coef_)) + NA )
		}
		
		T0 = - private$.Y / self$scale^2 + self$shape / self$scale
		T1 = base::digamma(self$shape) + base::log(self$scale) - base::log(private$.Y)
		
		## Compute gradient
		jac = private$.rhs$lhs$jacobian
		p   = 1
		if( !private$.rhs$lhs$fixed[["scale"]] )
		{
			jac[p,,] = jac[p,,] * as.vector(T0)
		}
		if( !private$.rhs$lhs$fixed[["shape"]] )
		{
			jac[p,,] = jac[p,,] * as.vector(T1)
		}
		
		return( base::apply( jac , 3 , base::sum ) )
	}
	##}}}
	
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
	},
	
	shape = function(value)
	{
		return(private$.rhs$lhs$values$shape)
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
	#' @field shape [numeric] Scale parameter
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new Gamma object.
    #' @param method [string] method used to fit
	#' @return A new `Gamma` object.
	initialize = function( method = "mle" )
	{
		super$initialize( base::c("scale","shape") , method )
	}
	
	##}}}
	
	)
	
	##}}}
	
)




