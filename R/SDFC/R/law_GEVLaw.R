
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
	valid = (Z > 0) & (1 + shape * Z > 0)
	
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
		out[!valid] = -Inf
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
	loc   = if( length(loc)   == size_p ) loc   else base::rep( loc[1]   , length(p) )
	scale = if( length(scale) == size_p ) scale else base::rep( scale[1] , length(p) )
	shape = if( length(shape) == size_p ) shape else base::rep( shape[1] , length(p) )
	
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
	}
	##}}}
	
#	fit = function( Y , loc_cov = NULL , scale_cov = NULL , shape_cov = NULL ) ##{{{
#	{
#		self$Y    = Y
#		self$size = length(Y)
#		
#		## Design matrix
#		size = length(Y)
#		self$loc_design   = base::cbind( base::rep(1,self$size) , loc_cov )
#		self$scale_design = base::cbind( base::rep(1,self$size) , scale_cov )
#		self$shape_design = base::cbind( base::rep(1,self$size) , shape_cov )
#		self$nloc         = base::ncol(self$loc_design)
#		self$nscale       = base::ncol(self$scale_design)
#		self$nshape       = base::ncol(self$shape_design)
#		self$ncov         = self$nloc + self$nscale + self$nshape
#		
#		## Initial condition
#		param_init = self$find_init()
#		
#		## Optimization
#		self$optim_res  = stats::optim( param_init , fn = self$optim_function , gr = self$gradient_optim_function , method = self$method , hessian = TRUE )
#		self$loc_coef   = self$optim_res$par[1:self$nloc]
#		self$scale_coef = self$optim_res$par[(self$nloc+1):(self$nloc+self$nscale)]
#		self$shape_coef = self$optim_res$par[(self$nloc+self$nscale+1):(self$ncov)]
#		
#		## Set scale and shape
#		self$update_param( self$optim_res$par )
#	},
#	##}}}
#	
#	link = function( x ) ##{{{
#	{
#		if( self$use_phi )
#		{
#			return( base::exp(x) )
#		}
#		else
#		{
#			return(x)
#		}
#	},
#	##}}}
#	
#	link_inv = function( x ) ##{{{
#	{
#		if( self$use_phi )
#		{
#			return( base::log(x) )
#		}
#		else
#		{
#			return(x)
#		}
#	},
#	##}}}
#	
#	update_param = function( param ) ##{{{
#	{
#		## Extract coefficients from param
#		loc_coef   = param[1:self$nloc]
#		scale_coef = param[(self$nloc+1):(self$nloc+self$nscale)]
#		shape_coef = param[(self$nloc+self$nscale+1):(self$ncov)]
#		
#		## Set scale and shape
#		self$loc   = self$loc_design %*% loc_coef
#		self$scale = self$link( self$scale_design %*% scale_coef )
#		self$shape = self$shape_design %*% shape_coef
#	},
#	##}}}
#	
#	find_init_extRemes = function() ##{{{
#	{
#		## LMoments
#		lmom1 = SDFC::lmoments( self$Y , 1 )
#		lmom2 = SDFC::lmoments( self$Y , 2 )
#		lmom3 = SDFC::lmoments( self$Y , 3 )
#		
#		tau3  = lmom3 / lmom2
#		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
#		kappa = 7.8590 * co + 2.9554 * co**2
#		g     = base::gamma( 1. + kappa )
#		
#		init_loc   = numeric( self$nloc )
#		init_scale = numeric( self$nscale )
#		init_shape = numeric( self$nshape )
#		
#		init_scale[1] = lmom2 * kappa / ( ( 1 - 2^(- kappa) ) * g )
#		init_loc[1]   = lmom1 - init_scale[1] * (1 - g) / kappa
#		init_shape[1] = - kappa
#		init_scale[1] = self$link_inv(init_scale[1])
#		
#		param   = base::c(init_loc,init_scale,init_shape)
#		test_ll = self$optim_function(param)
#		
#		## Moments
#		m     = base::mean(self$Y)
#		s     = base::sqrt(6) * stats::sd(self$Y) / base::pi
#		
#		init_loc2   = numeric( self$nloc )
#		init_scale2 = numeric( self$nscale )
#		init_shape2 = numeric( self$nshape )
#		
#		init_loc2[1]   = m - 0.57722 * s
#		init_scale2[1] = base::log(s)
#		init_shape2[1] = 1e-8
#		init_scale2[1] = self$link_inv(init_scale2[1])
#		
#		param2   = base::c(init_loc2,init_scale2,init_shape2)
#		test_ll2 = self$optim_function(param2)
#		
#		if( !is.finite(test_ll) && !is.finite(test_ll2) )
#		{
#			init_loc3      = numeric( self$nloc )
#			init_scale3    = numeric( self$nscale )
#			init_shape3    = numeric( self$nshape )
#			init_loc3[1]   = 0
#			init_scale3[1] = 1
#			init_shape3[1] = 0.01
#			param3 = base::c(init_loc3,init_scale3,init_shape3)
#			test_ll3 = self$optim_function(param3)
#			return( list( param = param3 , ll = test_ll3 ) )
#		}
#		
#		if( test_ll < test_ll2 )
#			return( list( param = param , ll = test_ll ) )
#		else
#			return( list( param = param2 , ll = test_ll2 ) )
#	}, ##}}}
#	
#	find_init_quantiles = function() ##{{{
#	{
#		init_loc   = numeric( self$nloc )
#		init_scale = numeric( self$nscale )
#		init_shape = numeric( self$nshape )
#		
#		## Fit loc
#		if( self$nloc == 1)
#		{
#			rvY = SDFC::rv_histogram$new( self$Y )
#			init_loc[1] = rvY$icdf( base::exp(-1) )
#			loc = base::rep( init_loc[1] , length(self$Y) )
#		}
#		else
#		{
#			reg = SDFC::QuantileRegression$new( base::c( base::exp(-1) ) )
#			reg$fit( self$Y , self$loc_design[,2:self$nloc] )
#			init_loc = reg$coef()
#			loc      = reg$predict()
#		}
#		
#		## Fit shape	
#		lmom1 = SDFC::lmoments( self$Y - loc , 1 )
#		lmom2 = SDFC::lmoments( self$Y - loc , 2 )
#		lmom3 = SDFC::lmoments( self$Y - loc , 3 )
#		
#		tau3  = lmom3 / lmom2
#		co    = 2. / ( 3. + tau3 ) - base::log(2) / base::log(3)
#		kappa = 7.8590 * co + 2.9554 * co^2
#		init_shape[1] = - kappa
#		
#		## Fit scale
#		qscale = base::c(0.25,0.5,0.75)
#		if( self$nscale == 1 )
#		{
#			rvY = SDFC::rv_histogram$new( self$Y )
#			coef = - kappa / ( ( -np.log(qscale) )^kappa - 1 )
#			init_scale[1] = self$link_inv( base::mean( rvY$icdf( qscale ) * coef ) )
#		}
#		else
#		{
#			reg = QuantileRegression$new( qscale )
#			reg$fit( self$Y - loc , self$scale_design[,2:self$nscale] )
#			fshape = base::rep( -kappa , length(self$Y) )
#			coef = matrix( NA , nrow = length(qscale) , ncol = length(self$Y) )
#			for( i in 1:length(qscale) )
#			{
#				coef[i,] = fshape / ( (-base::log(qscale[i]))^(-fshape) -1 )
#			}
#			fscale = self$link_inv( base::apply( reg$predict() * base::t(coef) , 1 , base::mean ) )
#			lm = stats::lm( fscale ~ self$scale_design[,2:self$nscale] )
#			init_scale = as.vector( lm$coefficients )
#		}
#		
#		param   = base::c(init_loc,init_scale,init_shape)
#		test_ll = self$optim_function(param)
#		return( list( param = param , ll = test_ll ) )
#	}, ##}}}
#	
#	find_init = function() ##{{{
#	{
#		method_ext = self$find_init_extRemes()
#		method_qua = self$find_init_quantiles()
#		
#		if( method_ext$ll < method_qua$ll )
#		{
#			return(method_ext$param)
#		}
#		else
#		{
#			return(method_qua$param)
#		}
#	},
#	##}}}
#	
#	negloglikelihood = function() ##{{{
#	{
#		## Impossible scale
#		if( base::any( self$scale <= 0 ) )
#			return(Inf)
#		
#		## Fuck exponential case
#		zero_shape = ( base::abs(self$shape) < 1e-10 )
#		if( base::any(zero_shape) )
#			self$shape[zero_shape] = -1e-10
#		
#		##
#		Z = 1 + self$shape * ( self$Y - self$loc ) / self$scale
#		
#		if( base::any(Z <= 0) )
#			return( Inf )
#		
#		res = base::sum( ( 1. + 1. / self$shape ) * base::log(Z) + Z^( - 1. / self$shape ) + base::log(self$scale) )
#		
#		if( is.finite(res) )
#			return(res)
#		else
#			return(Inf)
#	},
#	##}}}
#	
#	optim_function = function( param )##{{{
#	{
#		self$update_param(param)
#		return( self$negloglikelihood() )
#	},
#	##}}}
#	
#	logZafun = function( Z , alpha ) ##{{{
#	{
#		return( alpha * base::log( 1. + self$shape * Z ) )
#	},
#	##}}}
#	
#	Zafun = function( Z , alpha ) ##{{{
#	{
#		return( base::exp( self$logZafun( Z , alpha ) ) )
#	},
#	##}}}
#	
#	gradient_optim_function = function( param ) ##{{{
#	{
#		self$update_param(param)
#		
#		## Impossible
#		if( base::any( 1. + self$shape * ( self$Y - self$loc ) / self$scale <= 0 ) )
#			return( numeric( self$ncov ) + NA )
#		
#		## Usefull values
#		Z      = ( self$Y - self$loc ) / self$scale
#		ishape = 1. / self$shape
#		Za1    = self$Zafun( Z , 1. )
#		Zamsi  = self$Zafun( Z , - ishape ) ## Za of Minus Shape Inverse
#		
#		## Vectors
#		phi_vect   = if( self$use_phi ) 1. else 1. / self$scale
#		loc_vect   = (Zamsi - 1 - self$shape) / ( self$scale * Za1 )
#		scale_vect = phi_vect * ( 1. + Z * ( Zamsi - 1 - self$shape ) / Za1 )
#		shape_vect = ( Zamsi - 1. ) * base::log(Za1) * ishape^2 + ( 1. + ishape - ishape * Zamsi ) * Z / Za1
#		
#		## Gradients
#		grad_loc   = base::t(self$loc_design)   %*% loc_vect 
#		grad_scale = base::t(self$scale_design) %*% scale_vect
#		grad_shape = base::t(self$shape_design) %*% shape_vect 
#		
#		return( base::c(grad_loc,grad_scale,grad_shape) )
#	}
	##}}}
	
	),
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	loc_   = NULL,
	scale_ = NULL,
	shape_ = NULL
	
	)
)
##}}}


