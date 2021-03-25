
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


#' Normal
#'
#' @description
#' Normal distribution
#'
#' @details
#' Class of Normal distribution
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
#' @examples
#' ##
#'
#' @export
Normal = R6::R6Class( "Normal" ,
	
	inherit = AbstractLaw,
	
	## Private list
	##=============
	##{{{
	
	private = list(
	
	fit_moments = function() ##{{{
	{
		coef_ = numeric(length(self$coef_))
		
		## Find loc
		if( !private$.rhs$lhs$fixed[["loc"]] )
		{
			X_loc = self$rhs$c_global[["loc"]]
			coef_[1:private$.rhs$l_global$sizec("loc")] = SDFC::mean( private$.Y , X_loc , private$.rhs$l_global$linkc("loc")$l , value = FALSE )
			self$coef_ = coef_
		}
		## Find scale
		if( !private$.rhs$lhs$fixed[["scale"]] )
		{
			X_scale = self$rhs$c_global[["scale"]]
			ib = private$.rhs$l_global$sizec("loc")+1
			ie = length(self$coef_)
			coef_[ib:ie] = SDFC::std( private$.Y , X_scale , self$loc , private$.rhs$l_global$linkc("scale")$l , value = FALSE )
			self$coef_ = coef_
		}
		
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
		
		if( !base::all(self$scale > 0) )
			return(Inf)
		if( !base::all(is.finite(self$scale)) )
			return(Inf)
		
		scale2 = self$scale^2
		
		nlll = base::sum(base::log(scale2)) / 2 + base::sum( (private$.Y - self$loc)^2 / scale2 ) / 2
		return(nlll)
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
	
	##}}}
	
	## Init ##{{{
	
	#' @description
    #' Create a new Normal object.
    #' @param method [string] method used to fit
	#' @return A new `Normal` object.
	initialize = function( method = "mle" )
	{
		super$initialize( base::c("loc","scale") , method )
	}
	
	##}}}
	
	)
	
	##}}}
	
)




