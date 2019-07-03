
##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SDFC (Statistical    ##
## Distribution Fit with Covariates) library. This library makes it possible    ##
## to regress the parameters of some statistical law with co-variates.          ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SDFC (Statistical Distribution Fit with Covariates). Cette librairie         ##
## permet de calculer de regresser les parametres de lois statistiques selon    ##
## plusieurs co-variables                                                       ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################


########################################################################################################################
##                                                                                                                    ##
## Generalized Extreme Value function like R                                                                          ##
##                                                                                                                    ##
########################################################################################################################

dgev = function( x , loc = 0 , scale = 1 , shape = 0 , log = FALSE ) ##{{{
{
	size_x = length(x)
	loc   = if( length(loc)   == size_x ) loc   else base::rep( loc[1]   , size_x ) 
	scale = if( length(scale) == size_x ) scale else base::rep( scale[1] , size_x ) 
	shape = if( length(shape) == size_x ) shape else base::rep( shape[1] , size_x ) 
	
	
	Z     = ( x - loc ) / scale
	valid = (1 + shape * Z > 0)
	shape_zero  = ( base::abs(shape) < 1e-10 )
	cshape_zero = !shape_zero
	
	TX = numeric(length(x)) + NA
	
	if( base::any(shape_zero) )
	{
		TX[shape_zero] = base::exp( - Z[shape_zero] )
	}
	if( base::any(cshape_zero) )
	{
		TX[cshape_zero] = ( 1 + shape[cshape_zero] * Z[cshape_zero] )^( - 1. / shape[cshape_zero] )
	}
	
	out = TX^( shape + 1  ) * base::exp( - TX ) / scale
	if( base::any(!valid) )
	{
		out[!valid] = 0
	}
	
	if( log )
		return(base::log(out))
	else
		return(out)

}
##}}}

pgev = function( q , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE ) ##{{{
{
	if( !lower.tail )
	{
		return( 1. - pgev( q , loc , scale , shape , lower.tail = TRUE ) )
	}
	
	size_q = base::length(q)
	loc   = if( length(loc)   == size_q ) loc   else base::rep( loc[1]   , size_q )
	scale = if( length(scale) == size_q ) scale else base::rep( scale[1] , size_q )
	shape = if( length(shape) == size_q ) shape else base::rep( shape[1] , size_q )
	
	shape_zero  = ( base::abs(shape) < 1e-10 )
	cshape_zero = !shape_zero
	
	Z = ( q - loc ) / scale
	out = numeric(size_q) + NA
	
	if( base::any(shape_zero) )
	{
		out[shape_zero] = base::exp( - base::exp( - Z[shape_zero] ) )
	}
	
	if( base::any(cshape_zero) )
	{
		out[cshape_zero] = base::exp( - ( 1. + shape[cshape_zero] * Z[cshape_zero] )^( - 1. / shape[cshape_zero] ) )
	}
	
	valid = (1 + shape * Z > 0)
	if( base::any(!valid) )
	{
		out[(shape > 0) & !valid] = 0
		out[(shape < 0) & !valid] = 1
	}
	return(out)
}
##}}}

qgev = function( p , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE ) ##{{{
{
	if( !lower.tail )
	{
		return( qgev( 1. - p , loc , scale , shape , lower.tail = TRUE ) )
	}
	
	size_p = base::length(p)
	loc   = if( length(loc)   == size_p ) loc   else base::rep( loc[1]   , length(p) )
	scale = if( length(scale) == size_p ) scale else base::rep( scale[1] , length(p) )
	shape = if( length(shape) == size_p ) shape else base::rep( shape[1] , length(p) )
	
	shape_zero  = ( base::abs(shape) < 1e-10 )
	cshape_zero = !shape_zero
	
	out = numeric(length(p)) + NA
	if( base::any(shape_zero) )
	{
		out[shape_zero] = loc[shape_zero] - scale[shape_zero] * base::log( - base::log(p[shape_zero]) )
	}
	
	if( base::any(cshape_zero) )
	{
		out[cshape_zero] = loc[cshape_zero] + scale[cshape_zero] * ( ( - base::log(p[cshape_zero]) )^(- shape[cshape_zero]) - 1. ) / shape[cshape_zero]
	}
	
	return(out)
}
##}}}

rgev = function( n = 1 , loc = 0 , scale = 1 , shape = 0 ) ##{{{
{
	p = stats::runif( n = n )
	return( qgev( p , loc , scale , shape ) )
}
##}}}



########################################################################################################################
##                                                                                                                    ##
## Generalized Extreme value class                                                                                    ##
##                                                                                                                    ##
########################################################################################################################


#' GEVLaw (Generalized Extreme Value)
#'
#' Class to generate, fit, and use a generalized extreme values law.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param loc  [vector]
#'        Location parameters
#' @param scale  [vector]
#'        Scale parameters
#' @param shape  [vector]
#'        Shape parameters
#' @param scale_cov  [matrix]
#'        Scale covariate for fit
#' @param shape_cov  [matrix]
#'        Shape covariate for fit
#' @param use_phi [bool]
#'        Use exponential function as link function for scale
#' @param method [string]
#'        Optimization method, default "BFGS"
#' @param verbose [bool]
#'        Print warning and error message
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(loc,scale,shape,use_phi,method,verbose)}}{Initialize GEV law with code{GEVLaw}}
#'   \item{\code{rvs(size,loc,scale,shape)}}{Random number generator}.
#'   \item{\code{density(x,loc,scale,shape)}}{Density}.
#'   \item{\code{cdf(Y,loc,scale,shape)}}{Cumulative distribution function.}.
#'   \item{\code{icdf(p,loc,scale,shape)}}{Inverse of Cumulative distribution function.}.
#'   \item{\code{sf(Y,loc,scale,shape)}}{Survival function.}.
#'   \item{\code{isf(p,loc,scale,shape)}}{Inverse of survival function.}.
#'   \item{\code{fit(Y,loc,scale_cov,shape_cov)}}{Fit the GEV law}.
#' }
#' @examples
#' ## Data
#' @export
GEVLaw = R6::R6Class( "GEVLaw" , ##{{{
	
	inherit = AbstractLaw,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	loc          = NULL,
	scale        = NULL,
	shape        = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( method = "MLE" , link_fct_loc = IdLinkFct$new() , link_fct_scale = IdLinkFct$new() , link_fct_shape = IdLinkFct$new() , n_bootstrap = 0 , alpha = 0.05 ) ##{{{
	{
		super$initialize( method , n_bootstrap , alpha )
		
		self$loc       = NULL
		self$scale     = NULL
		self$shape     = NULL
		
		private$loc_   = LawParam$new( linkFct = link_fct_loc   , kind = "loc"   )
		private$scale_ = LawParam$new( linkFct = link_fct_scale , kind = "scale" )
		private$shape_ = LawParam$new( linkFct = link_fct_shape , kind = "shape" )
	},
	##}}}
	
	
	###############
	## Functions ##
	###############
	
	fit = function( Y , loc_cov = NULL , scale_cov = NULL , shape_cov = NULL , floc = NULL , fscale = NULL , fshape = NULL ) ##{{{
	{
		Y = as.vector(Y)
		private$size_ = length(Y)
		
		##=> Bootstrap here
		if( self$n_bootstrap > 0 )
		{
			if( !is.null(loc_cov) && !is.matrix(loc_cov) )
				loc_cov = matrix( loc_cov , nrow = private$size_ , ncol = 1 )
			if( !is.null(scale_cov) && !is.matrix(scale_cov) )
				scale_cov = matrix( scale_cov , nrow = private$size_ , ncol = 1 )
			if( !is.null(shape_cov) && !is.matrix(shape_cov) )
				shape_cov = matrix( shape_cov , nrow = private$size_ , ncol = 1 )
			
			self$coefs_bootstrap = base::c()
			
			for( i in 1:self$n_bootstrap )
			{
				idx = base::sample( 1:private$size_ , private$size_ , replace = TRUE )
				loc_cov_bs   = if( is.null(loc_cov) )   loc_cov   else loc_cov[idx,]
				scale_cov_bs = if( is.null(scale_cov) ) scale_cov else scale_cov[idx,]
				shape_cov_bs = if( is.null(shape_cov) ) shape_cov else shape_cov[idx,]
				floc_bs      = if( is.null(floc) || length(floc) == 1 )     floc      else floc[idx]
				fscale_bs    = if( is.null(fscale) || length(fscale) == 1 ) fscale    else fscale[idx]
				fshape_bs    = if( is.null(fshape) || length(fshape) == 1 ) fshape    else fshape[idx]
				
				private$fit_( Y[idx] , loc_cov_bs , scale_cov_bs , shape_cov_bs , floc_bs , fscale_bs , fshape_bs )
				self$coefs_bootstrap = base::rbind( self$coefs_bootstrap , self$coef_ )
			}
			self$confidence_interval = base::apply( self$coefs_bootstrap , 2 , stats::quantile , probs = base::c( self$alpha / 2. , 1. - self$alpha / 2. ) )
		}
		
		private$fit_( Y , loc_cov , scale_cov , shape_cov , floc , fscale , fshape )
	}
	##}}}
	
	),
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	loc_   = NULL,
	scale_ = NULL,
	shape_ = NULL,
	
	
	###############
	## Functions ##
	###############
	
	fit_ = function( Y , loc_cov = NULL , scale_cov = NULL , shape_cov = NULL , floc = NULL , fscale = NULL , fshape = NULL ) ##{{{
	{
		private$Y_    = as.vector(Y)
		
		private$loc_$init(   X = loc_cov   , fix_values = floc   , size = private$size_ )
		private$scale_$init( X = scale_cov , fix_values = fscale , size = private$size_ )
		private$shape_$init( X = shape_cov , fix_values = fshape , size = private$size_ )
		
		if( self$method == "moments" )
		{
			private$fit_moments()
		}
		else if( self$method == "lmoments" )
		{
			private$fit_lmoments()
		}
		else if( self$method == "quantiles" )
		{
			private$fit_quantiles()
		}
		else
		{
			private$fit_mle()
		}
		
		self$coef_ = private$concat_param()
	},
	##}}}
	
	fit_moments = function() ##{{{
	{
		m = base::mean(private$Y_)
		s = base::sqrt(6) * stats::sd(private$Y_) / base::pi
		
		iloc   = m - 0.57722 * s
		iscale = base::log(s)
		ishape = 1e-8
		
		private$loc_$set_intercept(   private$loc_$linkFct$inverse( iloc )     )
		private$scale_$set_intercept( private$scale_$linkFct$inverse( iscale ) )
		private$shape_$set_intercept( private$shape_$linkFct$inverse( ishape ) )
		
		
		private$loc_$update()
		private$scale_$update()
		private$shape_$update()
		self$loc   = private$loc_$valueLf()
		self$scale = private$scale_$valueLf()
		self$shape = private$shape_$valueLf()
	},
	##}}}
	
	fit_lmoments = function() ##{{{
	{
		lmom1 = np_lmoments( private$Y_ , 1 )
		lmom2 = np_lmoments( private$Y_ , 2 )
		lmom3 = np_lmoments( private$Y_ , 3 )
		
		tau3  = lmom3 / lmom2
		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
		kappa = 7.8590 * co + 2.9554 * co**2
		g     = base::gamma( 1. + kappa )
		
		iscale = lmom2 * kappa / ( (1 - 2^(-kappa)) * g )
		iloc   = lmom1 - iscale * (1 - g) / kappa
		ishape = - kappa
		
		private$loc_$set_intercept(   private$loc_$linkFct$inverse( iloc )     )
		private$scale_$set_intercept( private$scale_$linkFct$inverse( iscale ) )
		private$shape_$set_intercept( private$shape_$linkFct$inverse( ishape ) )
		
		private$loc_$update()
		private$scale_$update()
		private$shape_$update()
		self$loc   = private$loc_$valueLf()
		self$scale = private$scale_$valueLf()
		self$shape = private$shape_$valueLf()
	},
	##}}}
	
	fit_quantiles = function() ##{{{
	{
		## Fit the loc
		if( private$loc_$not_fixed() )
		{
			if( private$loc_$size_ == 1 )
			{
				private$loc_$set_intercept( private$loc_$linkFct$inverse( as.vector(stats::quantile( private$Y_ , probs = base::exp(-1) )) ) )
			}
			else
			{
				loc = as.vector(np_quantile( private$Y_ , ltau = base::c(base::exp(-1)) , X = private$loc_$design_wo1() ))
				private$loc_$set_coef( np_mean( loc , X = private$loc_$design_wo1() , linkFct = private$loc_$linkFct , return_coef = TRUE ) )
			}
		}
		private$loc_$update()
		self$loc = private$loc_$valueLf()
		
		
		## Fit the scale
		if( private$scale_$not_fixed() )
		{
			probs = base::c( 0.25 , 0.5 , 0.75 )
			coef  = - 1. / base::log( - base::log(probs) )
			if( private$scale_$size_ == 1 )
			{
				private$scale_$set_intercept( private$scale_$linkFct$inverse(base::mean( as.vector( stats::quantile( private$Y_ - self$loc , probs ) ) * coef )) )
			}
			else
			{
				qreg = np_quantile( private$Y_ - self$loc , ltau = probs , X = private$scale_$design_wo1() )
				fscale = ( qreg %*% coef ) / length(probs)
				private$scale_$set_coef( np_mean( fscale , X = private$scale_$design_wo1() , linkFct = private$scale_$linkFct , return_coef = TRUE ) )
			}
		}
		private$scale_$update()
		self$scale = private$scale_$valueLf()
		
		
		## Fit the shape
		if( private$shape_$not_fixed() )
		{
			p0 = 0.1
			p1 = 0.9
			llp0 = base::log( - base::log(p0) )
			llp1 = base::log( - base::log(p1) )
			if( private$shape_$size_ == 1 )
			{
				q = stats::quantile( ( private$Y_ - self$loc ) / self$scale , base::c(p0,p1) )
				kappa = q[1] / q[2]
				sh = 2 * (llp0 - kappa * llp1 ) / ( llp0^2 - kappa * llp1^2 )
				private$shape_$set_intercept( private$shape_$linkFct$inverse(sh) )
			}
			else
			{
				q = np_quantile( ( private$Y_ - self$loc ) / self$scale , ltau = base::c(p0,p1) , X = private$shape_$design_wo1() )
				kappa = q[,1] / q[,2]
				fshape = 2 * (llp0 - kappa * llp1 ) / ( llp0^2 - kappa * llp1^2 )
				private$shape_$set_coef( np_mean( fshape , X = private$shape_$design_wo1() , linkFct = private$shape_$linkFct , return_coef = TRUE ) )
			}
		}
		private$shape_$update()
		self$shape = private$shape_$valueLf()
	},
	##}}}
	
	fit_mle = function() ##{{{
	{
		private$fit_quantiles()
		param_init = private$concat_param()
		
		## Test for initial value
		nll  = private$optim_function(param_init)
		gnll = private$gradient_optim_function(param_init)
		
		if( !is.finite(nll) || !is.finite(gnll) )
		{
			private$shape_$set_coef( numeric( private$shape_$size_ ) )
			param_init = private$concat_param()
		}
		
		optim_result = stats::optim( param_init , fn = private$optim_function , gr = private$gradient_optim_function , method = "BFGS" )
		private$update_param( optim_result$par )
	},
	##}}}
	
	split_param = function( param )##{{{
	{
		param_loc   = NULL
		param_scale = NULL
		param_shape = NULL
		s0 = private$loc_$size_
		s1 = private$scale_$size_
		s2 = private$shape_$size_
		
		if( private$loc_$not_fixed() && private$scale_$not_fixed() && private$shape_$not_fixed() )
		{
			param_loc   = param[1:s0]
			param_scale = param[(s0+1):(s0+s1)]
			param_shape = param[(s0+s1+1):(s0+s1+s2)]
		}
		else if( private$loc_$not_fixed() && private$scale_$not_fixed() )
		{
			s0 = private$loc_$size_
			s1 = private$scale_$size_
			param_loc   = param[1:s0]
			param_scale = param[(s0+1):(s0+s1)]
		}
		else if( private$loc_$not_fixed() && private$shape_$not_fixed() )
		{
			s0 = private$loc_$size_
			s1 = private$shape_$size_
			param_loc   = param[1:s0]
			param_shape = param[(s0+1):(s0+s2)]
		}
		else if( private$scale_$not_fixed() && private$shape_$not_fixed() )
		{
			param_scale = param[1:s1]
			param_shape = param[(s1+1):(s1+s2)]
		}
		else if( private$loc_$not_fixed() )
		{
			param_loc = param
		}
		else if( private$scale_$not_fixed() )
		{
			param_scale = param
		}
		else if( private$shape_$not_fixed() )
		{
			param_shape = param
		}
		
		return( list( loc = param_loc , scale = param_scale , shape = param_shape ) )
	},
	##}}}
	
	concat_param = function()##{{{
	{
		param = NULL
		param_loc   = if( private$loc_$not_fixed() )   private$loc_$coef_   else NULL
		param_scale = if( private$scale_$not_fixed() ) private$scale_$coef_ else NULL
		param_shape = if( private$shape_$not_fixed() ) private$shape_$coef_ else NULL
		
		param = base::c( param_loc , param_scale , param_shape )

#		if( private$loc_$not_fixed() && private$scale_$not_fixed() && private$shape_$not_fixed() )
#		{
#			param = base::c( private$loc_$coef_ , private$scale_$coef_ , private$shape$coef_ )
#		}
#		else if( private$loc_$not_fixed() && private$scale_$not_fixed() )
#		{
#			param = base::c( private$loc_$coef_ , private$scale_$coef_ )
#		}
#		else if( private$loc_$not_fixed() && private$shape_$not_fixed() )
#		{
#			param = base::c( private$loc_$coef_ , private$shape$coef_ )
#		}
#		else if( private$scale_$not_fixed() && private$shape_$not_fixed() )
#		{
#			param = base::c( private$scale_$coef_ , private$shape$coef_ )
#		}
#		else if( private$loc_$not_fixed() )
#		{
#			param = private$loc_$coef_
#		}
#		else if( private$scale_$not_fixed() )
#		{
#			param = private$scale_$coef_
#		}
#		else if( private$shape_$not_fixed() )
#		{
#			param = private$shape_$coef_
#		}
		
		return( param )
	},
	##}}}
	
	update_param = function( param ) ##{{{
	{
		param_sp = private$split_param(param)
		private$loc_$set_coef(   param_sp$loc )
		private$scale_$set_coef( param_sp$scale )
		private$shape_$set_coef( param_sp$shape )
		private$loc_$update()
		private$scale_$update()
		private$shape_$update()
		self$loc   = private$loc_$valueLf()
		self$scale = private$scale_$valueLf()
		self$shape = private$shape_$valueLf()
	},
	##}}}
	
	optim_function = function( param )##{{{
	{
		private$update_param(param)
		return( private$negloglikelihood() )
	},
	##}}}
	
	negloglikelihood = function() ##{{{
	{
		## Impossible scale
		if( base::any( self$scale <= 0 ) )
			return(Inf)
		
		## Fuck exponential case
		zero_shape = ( base::abs(self$shape) < 1e-10 )
		if( base::any(zero_shape) )
			self$shape[zero_shape] = -1e-10
		
		##
		Z = 1 + self$shape * ( private$Y_ - self$loc ) / self$scale
		
		if( base::any(Z <= 0) )
			return( Inf )
		
		res = base::sum( ( 1. + 1. / self$shape ) * base::log(Z) + Z^( - 1. / self$shape ) + base::log(self$scale) )
		
		if( is.finite(res) )
			return(res)
		else
			return(Inf)
	},
	##}}}
	
	logZafun = function( Z , alpha ) ##{{{
	{
		return( alpha * base::log( 1. + self$shape * Z ) )
	},
	##}}}
	
	Zafun = function( Z , alpha ) ##{{{
	{
		return( base::exp( private$logZafun( Z , alpha ) ) )
	},
	##}}}
	
	gradient_optim_function = function( param ) ##{{{
	{
		private$update_param(param)
		
		## Impossible
		if( base::any( 1. + self$shape * ( private$Y_ - self$loc ) / self$scale <= 0 ) )
			return( numeric( length(param) ) + NA )
		
		## Usefull values
		Z      = ( private$Y_ - self$loc ) / self$scale
		ishape = 1. / self$shape
		Za1    = private$Zafun( Z , 1. )
		Zamsi  = private$Zafun( Z , - ishape ) ## Za of Minus Shape Inverse
		
		##
		grad = base::c()
		if( private$loc_$not_fixed() )
		{
			loc_vect = private$loc_$valueGrLf()   * ( Zamsi - 1 - self$shape ) / ( self$scale * Za1 )
			grad_loc = base::t(private$loc_$design_) %*% loc_vect
			grad     = base::c(grad,grad_loc)
		}
		if( private$scale_$not_fixed() )
		{
			scale_vect = private$scale_$valueGrLf() * ( 1. + Z * ( Zamsi - 1 - self$shape ) / Za1 ) / self$scale
			grad_scale = base::t(private$scale_$design_) %*% scale_vect
			grad       = base::c(grad,grad_scale)
		}
		if( private$shape_$not_fixed() )
		{
			shape_vect = private$shape_$valueGrLf() * ( ( Zamsi - 1. ) * base::log(Za1) * ishape^2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1 )
			grad_shape = base::t(private$shape_$design_) %*% shape_vect
			grad       = base::c(grad,grad_shape)
		}
		
		return(grad)
	}
	##}}}
	
	
	)
)
##}}}


