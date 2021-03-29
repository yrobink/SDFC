
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


#' GEV
#'
#' @description
#' GEV distribution
#'
#' @details
#' Class to fit a GEV law with covariates, available methods are:
#' 
#' moments  : use empirical estimator of mean and standard deviation to find
#'            loc and scale, possibly with least square regression if
#'            covariates are given
#'
#' lmoments : Use L-Moments estimation, only in stationary context
#'
#' lmoments_experimental: Use non-stationary L-Moments with Quantile
#'            Regression, experimental and not published, only
#'            used to find an initialization of MLE
#'
#' bayesian : Bayesian estimation, i.e. the coefficient fitted is the mean of
#'            n_mcmc_iteration sample draw from the posterior P(coef_ | Y)
#'
#' mle      : Maximum likelihood estimation
#'
#' 
#' Parameters:
#'
#' loc   : location parameter
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
GEV = R6::R6Class( "GEV" ,
	
	inherit = AbstractLaw,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	fit_moments = function() ##{{{
	{
		coef_ = numeric(length(self$coef_))
		
		m = base::mean(private$.Y)
		s = base::sqrt(6) * stats::sd(private$.Y) / base::pi
		
		iloc   = m - 0.57722 * s
		iscale = base::max( 0.1 , base::log(s) )
		ishape = 1e-8
		
		il_b  = 1
		il_e  = il_b + private$.rhs$l_global$sizec("loc") - 1
		isc_b = il_e + 1
		isc_e = isc_b + private$.rhs$l_global$sizec("scale") - 1
		ish_b = isc_e + 1
		ish_e = ish_b + private$.rhs$l_global$sizec("shape") - 1
		
		## Fit scale
		if( !private$.rhs$lhs$fixed[["scale"]] )
		{
			coef_[isc_b] = private$.rhs$l_global$linkc("scale")$l$inverse(iscale)
		}
		
		## Fit loc
		if( !private$.rhs$lhs$fixed[["loc"]] )
		{
			if( private$.rhs$lhs$fixed[["scale"]] )
			{
				iloc = m - 0.57722 * base::exp(self$scale)
				X_loc = self$rhs$c_global[["loc"]]
				coef_[il_b:il_e] = SDFC::mean( iloc , X_loc , private$.rhs$l_global$linkc("loc")$l , value = FALSE )
			}
			else
			{
				coef_[il_b] = private$.rhs$l_global$linkc("loc")$l$inverse(iloc)
			}
		}
		
		## Fit shape
		if( !private$.rhs$lhs$fixed[["shape"]] )
		{
			coef_[ish_b] = private$.rhs$l_global$linkc("shape")$l$inverse(ishape)
		}
		
		self$coef_ = coef_
	},
	##}}}
	
	fit_lmoments = function() ##{{{
	{
		coef_ = numeric(length(self$coef_))
		
		
		lmom = lmoments( private$.Y )
		
		tau3  = lmom[3] / lmom[2]
		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
		kappa = 7.8590 * co + 2.9554 * co**2
		g     = base::gamma( 1. + kappa )
		
		
		iscale = lmom[2] * kappa / ( (1 - 2^( - kappa ) ) * g )
		iloc   = lmom[1] - iscale * (1 - g) / kappa
		ishape = - kappa
		
		
		il_b  = 1
		il_e  = il_b + private$.rhs$l_global$sizec("loc") - 1
		isc_b = il_e + 1
		isc_e = isc_b + private$.rhs$l_global$sizec("scale") - 1
		ish_b = isc_e + 1
		ish_e = ish_b + private$.rhs$l_global$sizec("shape") - 1
		
		## Fit scale
		if( !private$.rhs$lhs$fixed[["scale"]] )
		{
			coef_[isc_b] = private$.rhs$l_global$linkc("scale")$l$inverse(iscale)
		}
		
		## Fit loc
		if( !private$.rhs$lhs$fixed[["loc"]] )
		{
			if( private$.rhs$lhs$fixed[["scale"]] )
			{
				iloc = lmom[1] - self$scale * (1 - g) / kappa
				X_loc = self$rhs$c_global[["loc"]]
				coef_[il_b:il_e] = SDFC::mean( iloc , X_loc , private$.rhs$l_global$linkc("loc")$l , value = FALSE )
			}
			else
			{
				coef_[il_b] = private$.rhs$l_global$linkc("loc")$l$inverse(iloc)
			}
		}
		
		## Fit shape
		if( !private$.rhs$lhs$fixed[["shape"]] )
		{
			coef_[ish_b] = private$.rhs$l_global$linkc("shape")$l$inverse(ishape)
		}
		
		self$coef_ = coef_
	},
	##}}}
	
	fit_lmoments_experimental = function() ##{{{
	{
		## First step, find lmoments
		c_Y  = matrix( 1 , nrow = length(private$.Y) )
		rank = 1
		for( c in private$.rhs$c_global )
		{
			if( is.matrix(c) )
			{
				for( i in 1:ncol(c) )
				{
					c_Y2  = base::cbind(c_Y,c)
					rank2 = base::qr(c_Y2)$rank
					if( rank2 > rank )
					{
						c_Y   = c_Y2
						rank = rank2
					}
				}
			}
		}
		if( rank == 1 )
		{
			c_Y = NULL
		}
		else
		{
			c_Y = c_Y[,2:rank]
		}
		lmom = SDFC::lmoments( private$.Y , c_Y )
		
		## Now find loc/scale/shape
		tau3  = lmom[,3] / lmom[,2]
		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
		shape = - 7.8590 * co - 2.9554 * co^2
		
		## Find scale
		gshape = base::gamma( 1 - shape )
		scale  = - lmom[,2] * shape / ( gshape * ( 1 - 2^shape ) )
		
		if( !(base::min(scale) > 0 ) )
		{
			idx = !(scale > 0)
			scale[idx] = 1e-3
		}
		
		## Find loc
		loc = lmom[,1] - scale * ( gshape - 1 ) / shape
		
		## Find coefs
		coefs = base::c()
		for( p in base::c("loc","scale","shape") )
		{
			if( !private$.rhs$lhs$fixed[[p]] )
			{
				X = self$rhs$c_global[[p]]
				coefs = base::c(coefs,SDFC::mean( private$.Y , X , private$.rhs$l_global$linkc(p)$l , value = FALSE ))
			}
		}
		
		self$coef_ = coefs
	},
	##}}}
	
	special_fit = function() ##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
		else if( self$method == "lmoments" )
			private$fit_lmoments()
		else if( self$method == "lmoments-experimental" )
			private$fit_lmoments_experimental()
	},
	##}}}
	
	init_MLE = function()##{{{
	{
		if( private$.rhs$l_global$special_fit_allowed )
		{
			private$fit_lmoments_experimental()
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
		
		if( !base::all(self$scale > 0) )
			return(Inf)
		if( !base::all(is.finite(self$scale)) )
			return(Inf)
		
		scale2 = self$scale^2
		
		nlll = base::sum(base::log(scale2)) / 2 + base::sum( (private$.Y - self$loc)^2 / scale2 ) / 2
		return(nlll)
	},
	##}}}
	
	gradient_nlll = function(coef) ##{{{
	{
		self$coef_ = coef
		
		loc   = self$loc
		scale = self$scale
		Z     = (private$.Y - loc) / scale
		
		## Compute gradient
		T0  = - Z / scale
		T1  = - private$.Y * Z / scale^2 + loc * Z / scale^2 + 1 / scale
		jac = private$.rhs$lhs$jacobian
		p   = 1
		if( !private$.rhs$lhs$fixed[["loc"]] )
		{
			jac[p,,] = jac[p,,] * as.vector(T0)
			p = p + 1
		}
		if( !private$.rhs$lhs$fixed[["scale"]] )
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
    #' Create a new GEV object.
    #' @param method [string] method used to fit
	#' @return A new `GEV` object.
	initialize = function( method = "mle" )
	{
		super$initialize( base::c("loc","scale","shape") , method )
	}
	
	##}}}
	
	)
	
	##}}}
	
)




