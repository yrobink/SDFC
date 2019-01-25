
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

#' GPDLaw (Generalized Pareto)
#'
#' Class to generate, fit, and use a generalized pareto law.
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
#'   \item{\code{new(loc,scale,shape,use_phi,method,verbose)}}{Initialize pareto law with code{GPDLaw}}
#'   \item{\code{rvs(size,loc,scale,shape)}}{Random number generator}.
#'   \item{\code{density(x,loc,scale,shape)}}{Density}.
#'   \item{\code{cdf(Y,loc,scale,shape)}}{Cumulative distribution function.}.
#'   \item{\code{icdf(p,loc,scale,shape)}}{Inverse of Cumulative distribution function.}.
#'   \item{\code{sf(Y,loc,scale,shape)}}{Survival function.}.
#'   \item{\code{isf(p,loc,scale,shape)}}{Inverse of survival function.}.
#'   \item{\code{fit(Y,loc,scale_cov,shape_cov)}}{Fit the GPD law}.
#' }
#' @examples
#' ## Data
#' size = 2000
#' data = SDFC::Dataset2(size)
#' t = data$t
#' X = data$X
#' Y = data$Y
#' 
#' ## Quantile regression for loc parameters
#' ltau = base::c( 0.05 , 0.95 )
#' qr = SDFC::QuantileRegression$new( ltau )
#' qr$fit( Y , X )
#' ## Yq[,1] is the lower loc of GPD, and Yq[,2] the upper loc.
#' Yq = if( qr$is_success() ) qr$predict() else NULL
#' 
#' ## Upper GPD
#' gpdU = SDFC::GPDLaw$new()
#' gpdU$fit( Y , loc = Yq[,2] , scale_cov = X )
#' Yu = gpdU$rvs(size)
#' print( gpdU$scale_coef_ ) ## Scale coef fitted
#' print( gpdU$shape_coef_ ) ## Shape coef fitted
#' 
#' ## Lower GPD
#' gpdL = SDFC::GPDLaw$new()
#' gpdL$fit( -Y , loc = -Yq[,1] , scale_cov = -X ) ## Minima of Y is the maxima of -Y
#' Yl = -gpdL$rvs(size)
#' @export
GPDLaw = R6::R6Class( "GPDLaw" , ##{{{
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	use_phi      = NULL,
	method       = NULL,
	verbose      = NULL,
	Y            = NULL,
	size         = NULL,
	loc          = NULL,
	scale        = NULL,
	shape        = NULL,
	scale_design = NULL,
	shape_design = NULL,
	ncov         = NULL,
	nscale       = NULL,
	nshape       = NULL,
	optim_res    = NULL,
	scale_coef   = NULL,
	shape_coef   = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( loc = 0 , scale = 1 , shape = -0.1 , use_phi = FALSE , method = "BFGS" , verbose = FALSE ) ##{{{
	{
		self$loc     = loc
		self$scale   = scale
		self$shape   = shape
		self$use_phi = use_phi
		self$method  = method
		self$verbose = verbose
	},
	##}}}
	
	rvs = function( size , loc = NULL , scale = NULL , shape = NULL ) ##{{{
	{
		Z   = stats::runif( size , min = 0 , max = 1 )
		
#		shape_zero  = ( shape == 0 )
#		cshape_zero = !shape_zero
#		
#		out = numeric(size) + NA
#		if( base::any(shape_zero) )
#		{
#			out[shape_zero] = loc[shape_zero] + stats::rexp( n , rate = 1. / scale[shape_zero] )
#		}
#		
#		if( base::any(cshape_zero) )
#		{
#			out[cshape_zero] = loc[cshape_zero] + scale[cshape_zero] * ( Z[cshape_zero]^( -shape[cshape_zero] ) - 1 ) / shape[cshape_zero]
#		}
		
		return(self$icdf( Z , loc = loc , scale = scale , shape = shape ))
	},
	##}}}
	
	density = function( x , loc = NULL , scale = NULL , shape = NULL ) ##{{{
	{
		loc   = if( is.null(loc) )   self$loc   else loc
		scale = if( is.null(scale) ) self$scale else scale
		shape = if( is.null(shape) ) self$shape else shape
		
		loc   = if( length(loc)   == 1 ) base::rep( loc   , length(x) ) else loc
		scale = if( length(scale) == 1 ) base::rep( scale , length(x) ) else scale
		shape = if( length(shape) == 1 ) base::rep( shape , length(x) ) else shape
		
		Z = ( x - loc ) / scale
		
		shape_zero  = ( shape == 0 )
		cshape_zero = !shape_zero
		valid       = (Z > 0) & (1 + shape * Z > 0)
		
		out = numeric(length(x)) + NA
		if( base::any(shape_zero) && base::any(valid) )
		{
			idx = valid[shape_zero]
			out[shape_zero][idx] = -base::log(scale[shape_zero][idx]) - Z[shape_zero][idx]
		}
		
		if( base::any(cshape_zero) && base::any(valid) )
		{
			idx = valid[cshape_zero]
			out[cshape_zero][idx] = -base::log(scale[cshape_zero][idx]) - base::log(1. + shape[cshape_zero][idx] * Z[cshape_zero][idx]) / shape[cshape_zero][idx]
		}
		
		if( base::any(!valid) )
		{
			out[!valid] = -Inf
		}
		return(base::exp(out))
	},
	##}}}
	
	cdf = function( Y , loc = NULL , scale = NULL , shape = NULL ) ##{{{
	{
		loc   = if( is.null(loc) )   self$loc   else loc
		scale = if( is.null(scale) ) self$scale else scale
		shape = if( is.null(shape) ) self$shape else shape
		
		loc   = if( length(loc)   == 1 ) base::rep( loc   , length(Y) ) else loc
		scale = if( length(scale) == 1 ) base::rep( scale , length(Y) ) else scale
		shape = if( length(shape) == 1 ) base::rep( shape , length(Y) ) else shape
		
		shape_zero  = ( shape == 0 )
		cshape_zero = !shape_zero
		
		Z = ( Y - loc ) / scale
		out = numeric(length(Y)) + NA
		if( base::any(shape_zero) )
		{
			out[shape_zero] = 1. - base::exp( - Z )
		}
		
		if( base::any(cshape_zero) )
		{
			out[cshape_zero] = 1. - ( 1. + shape * Z )^( -1. / shape )
		}
		
		return(out)
	},
	##}}}
	
	icdf = function( p , loc = NULL , scale = NULL , shape = NULL ) ##{{{
	{
		loc   = if( is.null(loc) )   self$loc   else loc
		scale = if( is.null(scale) ) self$scale else scale
		shape = if( is.null(shape) ) self$shape else shape
		
		loc   = if( length(loc)   == 1 ) base::rep( loc   , length(p) ) else loc
		scale = if( length(scale) == 1 ) base::rep( scale , length(p) ) else scale
		shape = if( length(shape) == 1 ) base::rep( shape , length(p) ) else shape
		
		shape_zero  = ( shape == 0 )
		cshape_zero = !shape_zero
		
		out = numeric(length(p)) + NA
		if( base::any(shape_zero) )
		{
			out[shape_zero] = loc - scale * base::log(1. - p)
		}
		
		if( base::any(cshape_zero) )
		{
			out[cshape_zero] = loc + scale * ( (1. - p)^(-shape) - 1 ) / shape
		}
		
		return(out)
	},
	##}}}
	
	sf = function( Y , loc = NULL , scale = NULL , shape = NULL ) ##{{{
	{
		return( 1 - self$cdf(Y) , loc = loc , scale = scale , shape = shape )
	},
	##}}}
	
	isf = function( p , loc = NULL , scale = NULL , shape = NULL ) ##{{{
	{
		return( self$icdf( 1. - p , loc = loc , scale = scale , shape = shape ) )
	},
	##}}}
	
	fit = function( Y , loc , scale_cov = NULL , shape_cov = NULL ) ##{{{
	{
		self$Y    = Y
		self$size = length(Y)
		self$loc  = if( length(loc) == length(Y) ) loc else base::rep(loc[1],length(Y))
		
		## Design matrix
		size = length(Y)
		self$scale_design = base::cbind( base::rep(1,self$size) , scale_cov )
		self$shape_design = base::cbind( base::rep(1,self$size) , shape_cov )
		self$nscale       = base::ncol(self$scale_design)
		self$nshape       = base::ncol(self$shape_design)
		self$ncov         = self$nscale + self$nshape
		
		## Initial condition
		param_init = self$find_init()
		
		## Optimization
		self$optim_res = stats::optim( param_init , fn = self$optim_function , gr = self$gradient_optim_function , method = self$method , hessian = TRUE )
		self$scale_coef = self$optim_res$par[1:self$nscale]
		self$shape_coef = self$optim_res$par[(self$nscale+1):(self$ncov)]
		
		## Set scale and shape
		self$update_param( self$optim_res$par )
	},
	##}}}
	
	link = function( x ) ##{{{
	{
		if( self$use_phi )
		{
			return( base::exp(x) )
		}
		else
		{
			return(x)
		}
	},
	##}}}
	
	link_inv = function( x ) ##{{{
	{
		if( self$use_phi )
		{
			return( base::log(x) )
		}
		else
		{
			return(x)
		}
	},
	##}}}
	
	update_param = function( param ) ##{{{
	{
		## Extract coefficients from param
		scale_coef = param[1:self$nscale]
		shape_coef = param[(self$nscale+1):(self$ncov)]
		
		## Set scale and shape
		self$scale = self$link( self$scale_design %*% scale_coef )
		self$shape = self$shape_design %*% shape_coef
	},
	##}}}
	
	find_init = function() ##{{{
	{
		## LMoments initial condition
		idx_excess = (self$Y > self$loc)
		excess = self$Y[idx_excess] - self$loc[idx_excess]
		lmo1     = LMoments1(excess)
		lmo2     = LMoments2(excess)
		itau     = lmo1 / lmo2
		scale_lm = lmo1 * ( itau - 1 )
		scale_lm = if( scale_lm > 0 ) scale_lm else 1e-8
		shape_lm = - ( itau - 2 )
		
		## Check negloglikelihood
		self$scale = base::rep( scale_lm , self$size )
		self$shape = base::rep( shape_lm , self$size )
		eval_lm = self$negloglikelihood()
		
		## MOMS initial condition
		scale_mm = base::sqrt( stats::var(excess) )
		scale_mm = if( scale_mm > 0 ) scale_mm else 1e-8
		shape_mm = -1e-8
		
		## Check negloglikelihood
		self$scale = base::rep( scale_mm , self$size )
		self$shape = base::rep( shape_mm , self$size )
		eval_mm = self$negloglikelihood()
		
		## Keep best
		init_scale = base::rep( 0 , self$nscale )
		init_shape = base::rep( 0 , self$nshape )
		if( is.finite(eval_mm) || is.finite(eval_lm) )
		{
			init_scale[1] = self$link_inv( if( eval_mm < eval_lm ) scale_mm else scale_lm )
			init_shape[1] = if( eval_mm < eval_lm ) shape_mm else shape_lm
		}
		return(base::c(init_scale,init_shape))
	},
	##}}}
	
	negloglikelihood = function() ##{{{
	{
		## Impossible scale
		if( base::any( self$scale <= 0 ) )
			return(Inf)
		
		## Fuck exponential case
		zero.shape = ( base::abs(self$shape) < 1e-10 )
		if( !is.null(zero.shape) )
		{
			self$shape[zero.shape] = -1e-10
		}
		
		##
		idx_excess = (self$Y > self$loc)
		loc   = self$loc[idx_excess]
		scale = self$scale[idx_excess]
		shape = self$shape[idx_excess]
		Z = 1. + shape * ( self$Y[idx_excess] - loc ) / scale
		
		if( base::any(Z <= 0) )
			return(Inf)
		
		res = base::sum( base::log( scale ) + base::log(Z) * ( 1 + 1. / shape ) )
		
		if( is.na(res) )
			return(Inf)
		return(res)
	},
	##}}}
	
	optim_function = function( param )##{{{
	{
		self$update_param(param)
		return( self$negloglikelihood() )
	},
	##}}}
	
	gradient_optim_function = function( param ) ##{{{
	{
		self$update_param(param)
		
		idx_excess = ( self$Y > self$loc )
		Y      = self$Y[idx_excess]
		loc    = self$loc[idx_excess]
		scale  = self$scale[idx_excess]
		shape  = self$shape[idx_excess]
		Y_zero   = Y - loc
		Z        = 1. + shape * Y_zero / scale
		exponent = 1. + 1. / shape
		
		grad_phi   = if( self$use_phi ) scale else 1.
		grad_scale = base::t(self$scale_design[idx_excess,]) %*% (grad_phi / scale - ( grad_phi * exponent * shape * Y_zero / (scale^2) ) / Z)
		grad_shape = if( base::all(Z>0) ) base::t(self$shape_design[idx_excess,]) %*% ( - base::log(Z) / (shape**2) + exponent * Y_zero / scale / Z) else base::rep(NaN,self$nshape)
		
#		phi.grad = if( use_phi ) scale else 1.
#		
#		coef = phi.grad / scale - ( phi.grad * exponent * shape * Y_zero / (scale^2) ) / Z
#		grad_scale = base::apply( scale_design , 2 , function(X) { return(base::sum(coef*X)) } )
#		
#		coef =  - base::log(Z) / (shape^2) + exponent * Y_zero / scale / Z 
#		grad_shape = base::apply( shape_design , 2 , function(X) { return(base::sum(coef*X)) } )
		
		return( base::c(grad_scale,grad_shape) )
	}
	##}}}
	
	)
)
##}}}
