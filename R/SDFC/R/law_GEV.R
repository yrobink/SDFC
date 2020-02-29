
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


## dgev {{{

#' dgev
#'
#' Density function of Generalized Extreme Value distribution
#'
#' @param x      [vector] Vector of values
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param log    [bool]   Return log of density if TRUE, default is FALSE
#'
#' @return [vector] Density of GEV at x
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' x = base::seq( -5 , 5 , length = 1000 )
#' y = dgev( x , loc = loc , scale = scale , shape = shape )
#' @export
dgev = function( x , loc = 0 , scale = 1 , shape = 0 , log = FALSE )
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

## pgev {{{

#' pgev
#'
#' Cumulative distribution function (or survival function) of Generalized Extreme Value distribution
#'
#' @param q      [vector] Vector of quantiles
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param lower.tail [bool] Return CDF if TRUE, else return survival function
#'
#' @return [vector] CDF (or SF) of GEV at x
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' x = base::seq( -5 , 5 , length = 1000 )
#' cdfx = pgev( x , loc = loc , scale = scale , shape = shape )
#' @export
pgev = function( q , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE )
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

## qgev {{{

#' qgev
#'
#' Inverse of CDF (or SF) function of Generalized Extreme Value distribution
#'
#' @param p      [vector] Vector of probabilities
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param lower.tail [bool] Return inverse of CDF if TRUE, else return inverse of survival function
#'
#' @return [vector] Inverse of CDF or SF of GEV for probabilities p
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' p = base::seq( 0.01 , 0.99 , length = 100 )
#' q = qgev( p , loc = loc , scale = scale , shape = shape )
#' @export
qgev = function( p , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE )
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

## rgev {{{

#' rgev
#'
#' Random value generator of Generalized Extreme Value distribution
#'
#' @param n      [int]    Numbers of values generated
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param log    [bool]   Return log of density if TRUE, default is FALSE
#'
#' @return [vector] Random value following a GEV(loc,scale,shape)
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' gev = rgev( 100 , loc = loc , scale = scale , shape = shape )
#' @export
rgev = function( n = 1 , loc = 0 , scale = 1 , shape = 0 )
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


## GEV {{{

#' GEV (Generalized Extreme Value distribution
#'
#' Class to fit a GEV law.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param method [string]
#'        Fit method, "moments", "lmoments", "lmoments_experimental", "MLE" and "bayesian" are available.
#' @param n_bootstrap [int]
#'        Number of bootstrap, default 0
#' @param alpha [float]
#'        Level of confidence interval, default 0.05
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(method,n_bootstrap,alpha)}}{Initialize Normal law with code{NormalLaw}}
#'   \item{\code{fit(Y,...)}}{Fit the Normal law}.
#' }
#' @examples
#' ## Start by generate non-stationary Normal dataset
#' size = 2000
#' c_data = SDFC::Dataset$covariates(size)
#' 
#' t       = c_data$t
#' X_loc   = c_data$X_loc
#' X_scale = c_data$X_scale
#' loc   = 0.5 + 2 * X_loc
#' scale =   1 + 2 * X_scale
#' Y = stats::rnorm( size , mean = loc , sd = scale )
#' 
#' ## Regression with MLE
#' law = SDFC::Normal$new( "mle" )
#' law$fit( Y , c_loc = X_loc , c_scale = X_scale )
#' 
#' ## Assuming scale is known (available for any covariates)
#' law = SDFC::Normal$new( "mle" )
#' law$fit( Y , c_loc = X_loc , f_scale = scale )
#' 
#' ## And if we want a link function
#' law = SDFC::Normal$new( "mle" )
#' law$fit( Y , c_loc = X_loc , c_scale = scale , l_scale = SDFC::ExpLink$new() )
#' 
#' ## If we do not give a parameter, it is assumed constant
#' law = SDFC::Normal$new( "mle" )
#' law$fit( Y , c_scale = X_scale )
#' 
#' @export
GEV = R6::R6Class( "GEV" ,##{{{
	
	inherit = AbstractLaw,
	
	##################
	## Private list ##
	##{{{
	
	private = list(
	
	## Arguments
	##==========
	
	## Methods
	##========
	
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
	
	
	fit_moments = function()##{{{
	{
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		m = base::mean(private$Y)
		s = base::sqrt(6) * stats::sd(private$Y) / base::pi
		
		iloc   = m - 0.57722 * s
		iscale = base::log(s)
		ishape = 1e-8
		
		if( !pscale$is_fix() )
			self$params$set_intercept( iscale , "scale" )
		
		if( !ploc$is_fix() )
		{
			if( pscale$is_fix() )
			{
				iloc = m - 0.57722 * base::exp(pscale$value)
				self$params$update_coef( np_mean( iloc , ploc$design_wo1() , value = FALSE , link = ploc$link ) , "loc" )
			}
			else
			{
				self$params$set_intercept( iloc , "loc" )
			}
		}
		
		if( !pshape$is_fix() )
			self$params$set_intercept( ishape , "shape" )
		
	},
	##}}}
	
	fit_lmoments = function() ##{{{
	{
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		lmom = np_lmoments( private$Y )
		
		tau3  = lmom[3] / lmom[2]
		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
		kappa = 7.8590 * co + 2.9554 * co^2
		g     = base::gamma( 1. + kappa )
		
		
		iscale = lmom[2] * kappa / ( ( 1 - 2^( - kappa ) ) * g )
		iloc   = lmom[1] - iscale * (1 - g) / kappa
		ishape = - kappa
		
		## Fit scale
		if( !pscale$is_fix() )
			self$params$set_intercept( iscale , "scale" )
		
		## Fit loc
		if( !ploc$is_fix() )
		{
			if( pscale$is_fix() )
			{
				iloc = lmom[1] - pscale$value * (1 - g) / kappa
				self$params$update_coef( np_mean( iloc , ploc$design_wo1() , value = FALSE , link = ploc$link ) , "loc" )
			}
			else
			{
				self$params$set_intercept( iloc , "loc" )
			}
		}
		
		## Fit shape
		if( !pshape$is_fix() )
			self$params$set_intercept( ishape , "shape" )
	},
	##}}}
	
	fit_lmoments_experimental = function()##{{{
	{
		c_Y = self$params$merge_covariate()
		if( is.null(c_Y) )
		{
			private$fit_lmoments()
			return()
		}
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		## First step, find lmoments
		lmom = np_lmoments( private$Y , c_Y = c_Y )
		
		## Find shape
		uni_shape_solver = function(tau)
		{
			interv = base::c( -1 , 1 )
			fct = function(x) { 3 / 2 + tau / 2 - ( 1 - 3^x ) / (1 - 2^x ) }
			while( base::prod(fct(interv)) > 0 )
			{
				interv = 2 * interv
			}
			opt = stats::uniroot( fct , interv )
			return( opt$root )
		}
		tau3 = matrix( lmom[,3] / lmom[,2] , ncol = 1 )
		shape = base::apply( tau3 , 1 , uni_shape_solver )
		
		## Find scale
		gshape = base::gamma( 1 - shape )
		scale = - lmom[,2] * shape / ( gshape * ( 1 - 2^shape ) )
		
		## Find loc
		loc = lmom[,1] - scale * ( gshape - 1 ) / shape
		
		if( !ploc$is_fix() )
			self$params$update_coef( np_mean( loc   , ploc$design_wo1()   , value = FALSE , link = ploc$link   ) , "loc"   )
		if( !pscale$is_fix() )
			self$params$update_coef( np_mean( scale , pscale$design_wo1() , value = FALSE , link = pscale$link ) , "scale" )
		if( !pshape$is_fix() )
			self$params$update_coef( np_mean( shape , pshape$design_wo1() , value = FALSE , link = pshape$link ) , "shape" )
	},
	##}}}
	
	initialization_mle = function()##{{{
	{
		private$fit_lmoments()
	},
	##}}}
	
	fit_ = function()##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
		else if( self$method == "lmoments" )
			private$fit_lmoments()
		else if( self$method == "lmoments-experimental" )
			private$fit_lmoments_experimental()
	},
	##}}}
	
	negloglikelihood = function( coef )##{{{
	{
		self$coef_ = coef
		## Impossible scale
		if( !base::all( self$scale > 0 ) )
			return(Inf)
		
		## Fuck exponential case
		shape = self$shape
		zero_shape = ( base::abs(shape) < 1e-10 )
		if( base::any(zero_shape) )
			shape[zero_shape] = -1e-10
		
		##
		Z = 1 + shape * ( private$Y - self$loc ) / self$scale
		
		if( !base::all(Z > 0) )
			return( Inf )
		
		res = base::sum( ( 1. + 1. / shape ) * base::log(Z) + Z^( - 1. / shape ) + base::log(self$scale) )
	
		if( is.finite(res) )
			return(res)
		else
			return(Inf)
	},
	##}}}
	
	gradient_nlll = function( coef ) ##{{{
	{
		self$coef_ = coef
		
		## Impossible
		if( !base::all( 1. + self$shape * ( private$Y - self$loc ) / self$scale > 0 ) )
			return( numeric( length(coef) ) + NA )
		
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		
		## Usefull values
		Z      = ( private$Y - self$loc ) / self$scale
		ishape = 1. / self$shape
		Za1    = private$Zafun( Z , 1. )
		Zamsi  = private$Zafun( Z , - ishape ) ## Za of Minus Shape Inverse
		
		grad = base::c()
		if( !ploc$is_fix() )
		{
			loc_vect = ploc$gradient() * ( Zamsi - 1 - self$shape ) / ( self$scale * Za1 )
			grad_loc = base::t(ploc$design_) %*% loc_vect
			grad     = base::c( grad , grad_loc )
		}
		if( !pscale$is_fix() )
		{
			scale_vect = pscale$gradient() * ( 1. + Z * ( Zamsi - 1 - self$shape ) / Za1 ) / self$scale
			grad_scale = base::t(pscale$design_) %*% scale_vect
			grad       = base::c( grad , grad_scale )
		}
		if( !pshape$is_fix() )
		{
			shape_vect = pshape$gradient() * ( ( Zamsi - 1. ) * base::log(Za1) * ishape^2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1 )
			grad_shape = base::t(pshape$design_) %*% shape_vect
			grad       = base::c( grad , grad_shape )
		}
		
		return( grad )
	}
	##}}}
	
	),
	##}}}
	##################
	
	#################
	## Public list ##
	##{{{
	
	public = list(
	
	## Arguments
	##==========
	
	## Constructor
	##============
	initialize = function( method = "mle" , n_bootstrap = 0 , alpha = 0.05 )
	{
		super$initialize( base::c( "loc" , "scale" , "shape" ) , method , n_bootstrap , alpha )
	}
	
	
	## Methods
	##========
	
	),
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	loc = function( l )##{{{
	{
		if( missing(l) )
			return( self$params$dparams_[["loc"]]$value )
	},
	##}}}
	
	scale = function( s )##{{{
	{
		if( missing(s) )
			return( self$params$dparams_[["scale"]]$value )
	},
	##}}}
	
	shape = function( s )##{{{
	{
		if( missing(s) )
			return( self$params$dparams_[["shape"]]$value )
	}
	##}}}
	
	)
	##}}}
	#################

)
##}}}
##}}}




## ==> ## GEVLaw {{{
## ==> 
## ==> GEVLaw = R6::R6Class( "GEVLaw" , 
## ==> 	
## ==> 	inherit = AbstractLaw,
## ==> 	
## ==> 	public = list(
## ==> 	
## ==> 	###############
## ==> 	## Arguments ##
## ==> 	###############
## ==> 	
## ==> 	loc    = NULL,
## ==> 	scale  = NULL,
## ==> 	shape  = NULL,
## ==> 	loc_   = NULL,
## ==> 	scale_ = NULL,
## ==> 	shape_ = NULL,
## ==> 	
## ==> 	
## ==> 	#################
## ==> 	## Constructor ##
## ==> 	#################
## ==> 	
## ==> 	initialize = function( method = "MLE" , link_fct_loc = IdLinkFct$new() , link_fct_scale = IdLinkFct$new() , link_fct_shape = IdLinkFct$new() , n_bootstrap = 0 , alpha = 0.05 ) ##{{{
## ==> 	{
## ==> 		super$initialize( method , n_bootstrap , alpha )
## ==> 		
## ==> 		self$loc       = NULL
## ==> 		self$scale     = NULL
## ==> 		self$shape     = NULL
## ==> 		
## ==> 		self$loc_   = LawParam$new( linkFct = link_fct_loc   , kind = "loc"   )
## ==> 		self$scale_ = LawParam$new( linkFct = link_fct_scale , kind = "scale" )
## ==> 		self$shape_ = LawParam$new( linkFct = link_fct_shape , kind = "shape" )
## ==> 	},
## ==> 	##}}}
## ==> 	
## ==> 	
## ==> 	###############
## ==> 	## Functions ##
## ==> 	###############
## ==> 	
## ==> 	fit = function( Y , loc_cov = NULL , scale_cov = NULL , shape_cov = NULL , floc = NULL , fscale = NULL , fshape = NULL ) ##{{{
## ==> 	{
## ==> 		Y = as.vector(Y)
## ==> 		private$size_ = length(Y)
## ==> 		
## ==> 		##=> Bootstrap here
## ==> 		if( self$n_bootstrap > 0 )
## ==> 		{
## ==> 			if( !is.null(loc_cov) && !is.matrix(loc_cov) )
## ==> 				loc_cov = matrix( loc_cov , nrow = private$size_ , ncol = 1 )
## ==> 			if( !is.null(scale_cov) && !is.matrix(scale_cov) )
## ==> 				scale_cov = matrix( scale_cov , nrow = private$size_ , ncol = 1 )
## ==> 			if( !is.null(shape_cov) && !is.matrix(shape_cov) )
## ==> 				shape_cov = matrix( shape_cov , nrow = private$size_ , ncol = 1 )
## ==> 			
## ==> 			self$coefs_bootstrap = base::c()
## ==> 			
## ==> 			for( i in 1:self$n_bootstrap )
## ==> 			{
## ==> 				idx = base::sample( 1:private$size_ , private$size_ , replace = TRUE )
## ==> 				loc_cov_bs   = if( is.null(loc_cov) )   loc_cov   else loc_cov[idx,]
## ==> 				scale_cov_bs = if( is.null(scale_cov) ) scale_cov else scale_cov[idx,]
## ==> 				shape_cov_bs = if( is.null(shape_cov) ) shape_cov else shape_cov[idx,]
## ==> 				floc_bs      = if( is.null(floc) || length(floc) == 1 )     floc      else floc[idx]
## ==> 				fscale_bs    = if( is.null(fscale) || length(fscale) == 1 ) fscale    else fscale[idx]
## ==> 				fshape_bs    = if( is.null(fshape) || length(fshape) == 1 ) fshape    else fshape[idx]
## ==> 				
## ==> 				private$fit_( Y[idx] , loc_cov_bs , scale_cov_bs , shape_cov_bs , floc_bs , fscale_bs , fshape_bs )
## ==> 				self$coefs_bootstrap = base::rbind( self$coefs_bootstrap , self$coef_ )
## ==> 			}
## ==> 			self$confidence_interval = base::apply( self$coefs_bootstrap , 2 , stats::quantile , probs = base::c( self$alpha / 2. , 1. - self$alpha / 2. ) )
## ==> 		}
## ==> 		
## ==> 		private$fit_( Y , loc_cov , scale_cov , shape_cov , floc , fscale , fshape )
## ==> 	}
## ==> 	##}}}
## ==> 	
## ==> 	),
## ==> 	
## ==> 	private = list(
## ==> 	
## ==> 	###############
## ==> 	## Arguments ##
## ==> 	###############
## ==> 	
## ==> 	
## ==> 	
## ==> 	###############
## ==> 	## Functions ##
## ==> 	###############
## ==> 	
## ==> 	fit_ = function( Y , loc_cov = NULL , scale_cov = NULL , shape_cov = NULL , floc = NULL , fscale = NULL , fshape = NULL ) ##{{{
## ==> 	{
## ==> 		private$Y_    = as.vector(Y)
## ==> 		
## ==> 		self$loc_$init(   X = loc_cov   , fix_values = floc   , size = private$size_ )
## ==> 		self$scale_$init( X = scale_cov , fix_values = fscale , size = private$size_ )
## ==> 		self$shape_$init( X = shape_cov , fix_values = fshape , size = private$size_ )
## ==> 		
## ==> 		if( self$method == "moments" )
## ==> 		{
## ==> 			private$fit_moments()
## ==> 		}
## ==> 		else if( self$method == "lmoments" )
## ==> 		{
## ==> 			private$fit_lmoments()
## ==> 		}
## ==> 		else if( self$method == "quantiles" )
## ==> 		{
## ==> 			private$fit_quantiles()
## ==> 		}
## ==> 		else
## ==> 		{
## ==> 			private$fit_mle()
## ==> 		}
## ==> 		
## ==> 		self$coef_ = private$concat_param()
## ==> 	},
## ==> 	##}}}
## ==> 	
## ==> 	fit_moments = function() ##{{{
## ==> 	{
## ==> 		m = base::mean(private$Y_)
## ==> 		s = base::sqrt(6) * stats::sd(private$Y_) / base::pi
## ==> 		
## ==> 		iloc   = m - 0.57722 * s
## ==> 		iscale = base::log(s)
## ==> 		ishape = 1e-8
## ==> 		
## ==> 		self$loc_$set_intercept(   self$loc_$linkFct$inverse( iloc )     )
## ==> 		self$scale_$set_intercept( self$scale_$linkFct$inverse( iscale ) )
## ==> 		self$shape_$set_intercept( self$shape_$linkFct$inverse( ishape ) )
## ==> 		
## ==> 		
## ==> 		self$loc_$update()
## ==> 		self$scale_$update()
## ==> 		self$shape_$update()
## ==> 		self$loc   = self$loc_$valueLf()
## ==> 		self$scale = self$scale_$valueLf()
## ==> 		self$shape = self$shape_$valueLf()
## ==> 	},
## ==> 	##}}}
## ==> 	
## ==> 	fit_lmoments = function() ##{{{
## ==> 	{
## ==> 		lmom1 = np_lmoments( private$Y_ , 1 )
## ==> 		lmom2 = np_lmoments( private$Y_ , 2 )
## ==> 		lmom3 = np_lmoments( private$Y_ , 3 )
## ==> 		
## ==> 		tau3  = lmom3 / lmom2
## ==> 		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
## ==> 		kappa = 7.8590 * co + 2.9554 * co**2
## ==> 		g     = base::gamma( 1. + kappa )
## ==> 		
## ==> 		iscale = lmom2 * kappa / ( (1 - 2^(-kappa)) * g )
## ==> 		iloc   = lmom1 - iscale * (1 - g) / kappa
## ==> 		ishape = - kappa
## ==> 		
## ==> 		self$loc_$set_intercept(   self$loc_$linkFct$inverse( iloc )     )
## ==> 		self$scale_$set_intercept( self$scale_$linkFct$inverse( iscale ) )
## ==> 		self$shape_$set_intercept( self$shape_$linkFct$inverse( ishape ) )
## ==> 		
## ==> 		self$loc_$update()
## ==> 		self$scale_$update()
## ==> 		self$shape_$update()
## ==> 		self$loc   = self$loc_$valueLf()
## ==> 		self$scale = self$scale_$valueLf()
## ==> 		self$shape = self$shape_$valueLf()
## ==> 	},
## ==> 	##}}}
## ==> 	
## ==> 	fit_quantiles = function() ##{{{
## ==> 	{
## ==> 		## Fit the loc
## ==> 		if( self$loc_$not_fixed() )
## ==> 		{
## ==> 			if( self$loc_$size_ == 1 )
## ==> 			{
## ==> 				self$loc_$set_intercept( self$loc_$linkFct$inverse( as.vector(stats::quantile( private$Y_ , probs = base::exp(-1) )) ) )
## ==> 			}
## ==> 			else
## ==> 			{
## ==> 				loc = as.vector(np_quantile( private$Y_ , ltau = base::c(base::exp(-1)) , X = self$loc_$design_wo1() ))
## ==> 				self$loc_$set_coef( np_mean( loc , X = self$loc_$design_wo1() , linkFct = self$loc_$linkFct , return_coef = TRUE ) )
## ==> 			}
## ==> 		}
## ==> 		self$loc_$update()
## ==> 		self$loc = self$loc_$valueLf()
## ==> 		
## ==> 		
## ==> 		## Fit the scale
## ==> 		if( self$scale_$not_fixed() )
## ==> 		{
## ==> 			probs = base::c( 0.25 , 0.5 , 0.75 )
## ==> 			coef  = - 1. / base::log( - base::log(probs) )
## ==> 			if( self$scale_$size_ == 1 )
## ==> 			{
## ==> 				self$scale_$set_intercept( self$scale_$linkFct$inverse(base::mean( as.vector( stats::quantile( private$Y_ - self$loc , probs ) ) * coef )) )
## ==> 			}
## ==> 			else
## ==> 			{
## ==> 				qreg = np_quantile( private$Y_ - self$loc , ltau = probs , X = self$scale_$design_wo1() )
## ==> 				fscale = ( qreg %*% coef ) / length(probs)
## ==> 				self$scale_$set_coef( np_mean( fscale , X = self$scale_$design_wo1() , linkFct = self$scale_$linkFct , return_coef = TRUE ) )
## ==> 			}
## ==> 		}
## ==> 		self$scale_$update()
## ==> 		self$scale = self$scale_$valueLf()
## ==> 		
## ==> 		
## ==> 		## Fit the shape
## ==> 		if( self$shape_$not_fixed() )
## ==> 		{
## ==> 			p0 = 0.1
## ==> 			p1 = 0.9
## ==> 			llp0 = base::log( - base::log(p0) )
## ==> 			llp1 = base::log( - base::log(p1) )
## ==> 			if( self$shape_$size_ == 1 )
## ==> 			{
## ==> 				q = stats::quantile( ( private$Y_ - self$loc ) / self$scale , base::c(p0,p1) )
## ==> 				kappa = q[1] / q[2]
## ==> 				sh = 2 * (llp0 - kappa * llp1 ) / ( llp0^2 - kappa * llp1^2 )
## ==> 				self$shape_$set_intercept( self$shape_$linkFct$inverse(sh) )
## ==> 			}
## ==> 			else
## ==> 			{
## ==> 				q = np_quantile( ( private$Y_ - self$loc ) / self$scale , ltau = base::c(p0,p1) , X = self$shape_$design_wo1() )
## ==> 				kappa = q[,1] / q[,2]
## ==> 				fshape = 2 * (llp0 - kappa * llp1 ) / ( llp0^2 - kappa * llp1^2 )
## ==> 				self$shape_$set_coef( np_mean( fshape , X = self$shape_$design_wo1() , linkFct = self$shape_$linkFct , return_coef = TRUE ) )
## ==> 			}
## ==> 		}
## ==> 		self$shape_$update()
## ==> 		self$shape = self$shape_$valueLf()
## ==> 	},
## ==> 	##}}}
## ==> 	
