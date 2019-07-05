
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
## Generalized Pareto function like R                                                                                 ##
##                                                                                                                    ##
########################################################################################################################

## dgpd {{{

#' dgpd
#'
#' Density function of Generalized Pareto Distribution
#'
#' @param x      [vector] Vector of values
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param log    [bool]   Return log of density if TRUE, default is FALSE
#'
#' @return [vector] Density of GPD at x
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' x = base::seq( -5 , 5 , length = 1000 )
#' y = dgpd( x , loc = loc , scale = scale , shape = shape )
#' @export
dgpd = function( x , loc = 0 , scale = 1 , shape = 0 , log = FALSE )
{
	size_x = length(x)
	loc   = if( length(loc)   == size_x ) loc   else base::rep( loc[1]   , size_x )
	scale = if( length(scale) == size_x ) scale else base::rep( scale[1] , size_x )
	shape = if( length(shape) == size_x ) shape else base::rep( shape[1] , size_x )
	
	Z = ( x - loc ) / scale
	
	shape_zero  = ( shape == 0 )
	cshape_zero = !shape_zero
	valid       = (Z > 0) & (1 + shape * Z > 0)
	
	out = numeric(length(x)) + NA
	if( base::any(shape_zero) )
	{
		idx = shape_zero & valid
		if( base::any(idx) )
		{
			out[idx] = -base::log(scale[idx]) - Z[idx]
		}
	}
	
	if( base::any(cshape_zero) )
	{
		idx = cshape_zero & valid
		if( base::any(idx) )
		{
			out[idx] = -base::log(scale[idx]) - base::log(1. + shape[idx] * Z[idx]) * ( 1 + 1 / shape[idx] )
		}
	}
	
	if( base::any(!valid) )
	{
		out[!valid] = -Inf
	}
	if( log )
		return(out)
	else
		return(base::exp(out))

}
##}}}

## pgpd {{{

#' pgpd
#'
#' Cumulative distribution function (or survival function) of Generalized Pareto distribution
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
#' cdfx = pgpd( x , loc = loc , scale = scale , shape = shape )
#' @export
pgpd = function( q , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE )
{
	if( !lower.tail )
	{
		return( 1 - pgpd( q , loc , scale , shape , lower.tail = TRUE ) )
	}
	
	size_q = length(q)
	loc   = if( length(loc)   == size_q ) loc   else base::rep( loc[1]   , size_q )
	scale = if( length(scale) == size_q ) scale else base::rep( scale[1] , size_q )
	shape = if( length(shape) == size_q ) shape else base::rep( shape[1] , size_q )
	
	shape_zero  = ( shape == 0 )
	cshape_zero = !shape_zero
	
	Z     = ( q - loc ) / scale
	valid = (Z > 0) & (1 + shape * Z > 0)
	
	out = numeric(size_q) + NA
	if( base::any(shape_zero) )
	{
		idx = shape_zero & valid
		if( base::any(idx) )
		{
			out[idx] = 1. - base::exp( - Z[idx] )
		}
	}
	
	if( base::any(cshape_zero) )
	{
		idx = cshape_zero & valid
		if( base::any(idx) )
		{
			out[idx] = 1. - ( 1. + shape[idx] * Z[idx] )^( -1. / shape[idx] )
		}
	}
	
	if( base::any(!valid) )
	{
		idx0 = (Z > 1) & !valid
		idx1 = !(Z > 1) & !valid
		if( base::any(idx0) )
			out[idx0] = 1
		if( base::any(idx1) )
			out[idx1] = 0
	}
	
	return(out)
}
##}}}

## qgpd {{{

#' qgpd
#'
#' Inverse of CDF (or SF) function of Generalized Pareto distribution
#'
#' @param p      [vector] Vector of probabilities
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param lower.tail [bool] Return inverse of CDF if TRUE, else return inverse of survival function
#'
#' @return [vector] Inverse of CDF or SF of GPD for probabilities p
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' p = base::seq( 0.01 , 0.99 , length = 100 )
#' q = qgpd( p , loc = loc , scale = scale , shape = shape )
#' @export
qgpd = function( p , loc = 0 , scale = 1 , shape = 0 , lower.tail = TRUE )
{
	if( !lower.tail )
		return( qgpd( 1 - p , loc , scale , shape , TRUE ) )
	
	## Test size
	size_p = length(p)
	loc   = if( length(loc)   == size_p ) loc   else base::rep( loc[1]   , size_p )
	scale = if( length(scale) == size_p ) scale else base::rep( scale[1] , size_p )
	shape = if( length(shape) == size_p ) shape else base::rep( shape[1] , size_p )
	
	## Difference shape == 0 and non zero
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
}
##}}}

## rgpd {{{

#' rgpd
#'
#' Random value generator of Generalized Pareto distribution
#'
#' @param n      [int]    Numbers of values generated
#' @param loc    [vector] Location parameter
#' @param scale  [vector] Scale parameter
#' @param shape  [vector] Shape parameter
#' @param log    [bool]   Return log of density if TRUE, default is FALSE
#'
#' @return [vector] Random value following a loc + GPD(scale,shape)
#'
#' @examples
#' ## Data
#' loc = 1
#' scale = 0.5
#' shape = -0.2
#' gev = rgpd( 100 , loc = loc , scale = scale , shape = shape )
#' @export
rgpd = function( n = 1 , loc = 0 , scale = 1 , shape = 0 ) 
{
	p = stats::runif( n = n )
	return( qgpd( p , loc , scale , shape ) )
}
##}}}


########################################################################################################################
##                                                                                                                    ##
## Generalized Pareto class                                                                                           ##
##                                                                                                                    ##
########################################################################################################################

## GPDLaw {{{

#' GPDLaw (Generalized Pareto Distribution)
#'
#' Class to fit a Generalized Pareto Distribution.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param method [string]
#'        Fit method, "moments", "lmoments", and "MLE" are available.
#' @param link_fct_scale [SDFC::LinkFct]
#'        Link function for scale parameter. Can be an element of SDFC, or a class based on SDFC::LinkFct
#' @param link_fct_shape [SDFC::LinkFct]
#'        Link function for shape parameter. Can be an element of SDFC, or a class based on SDFC::LinkFct
#' @param n_bootstrap [int]
#'        Number of bootstrap, default 0
#' @param alpha [float]
#'        Level of confidence interval, default 0.05
#' @param loc  [Vector]
#'        Location parameter, use quantile regression to fit it
#' @param scale_cov  [matrix or NULL]
#'        Scale covariate for fit
#' @param shape_cov  [matrix or NULL]
#'        Shape covariate for fit
#' @param fscale [vector or NULL]
#'        Value of scale if it is not necessary to fit
#' @param fshape [vector or NULL]
#'        Value of shape if it is not necessary to fit
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(method,link_fct_scale,link_fct_shape,n_bootstrap,alpha)}}{Initialize GPD law with code{GPDLaw}}
#'   \item{\code{fit(Y,loc,scale_cov,shape_cov,fscale,fshape)}}{Fit the GPD law}.
#' }
#' @examples
#' ## Data
#' size  = 2500
#' t = base::seq( 0 , 1 , length = size )
#' X0 = t^2
#' X2 = base::seq( -1 , 1 , length = size )
#' loc   = 0.5 + 1.5 * X0
#' scale = 0.1 + 0.1 * X0
#' shape = 0.3 * X2
#' 
#' Y = rgpd( n = size , loc = loc , scale = scale , shape = shape )
#' 
#' ## Fit
#' gpd = GPDLaw$new( method = "MLE" , n_bootstrap = 10 )
#' gpd$fit( Y , loc = loc , scale_cov = X0 , shape_cov = X2 )
#' 
#' gpd$loc   ## Loc fitted
#' gpd$scale ## Scale fitted
#' gpd$shape ## Shape fitted
#' @export
GPDLaw = R6::R6Class( "GPDLaw" , 
	
	inherit = AbstractLaw,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	loc    = NULL,
	scale  = NULL,
	shape  = NULL,
	scale_ = NULL,
	shape_ = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( method = "MLE" , link_fct_scale = SDFC::IdLinkFct$new() , link_fct_shape = SDFC::IdLinkFct$new() , n_bootstrap = 0 , alpha = 0.05 ) ##{{{
	{
		super$initialize( method , n_bootstrap , alpha )
		
		self$loc       = NULL
		self$scale     = NULL
		self$shape     = NULL
		
		self$scale_ = LawParam$new( linkFct = link_fct_scale , kind = "scale" )
		self$shape_ = LawParam$new( linkFct = link_fct_shape , kind = "shape" )
	},
	##}}}
	
	
	#############
	## Methods ##
	#############
	
	fit = function( Y , loc , scale_cov = NULL , shape_cov = NULL , fscale = NULL , fshape = NULL ) ##{{{
	{
		private$size_ = length(Y)
		
		## Bootstrap
		if( self$n_bootstrap > 0 )
		{
			if( !is.null(scale_cov) && !is.matrix(scale_cov) )
				scale_cov = matrix( scale_cov , nrow = private$size_ , ncol = 1 )
			if( !is.null(shape_cov) && !is.matrix(shape_cov) )
				shape_cov = matrix( shape_cov , nrow = private$size_ , ncol = 1 )
			
			self$coefs_bootstrap = base::c()
			
			for( i in 1:self$n_bootstrap )
			{
				idx = base::sample( 1:private$size_ , private$size_ , replace = TRUE )
				loc_bs       = loc[idx]
				scale_cov_bs = if( is.null(scale_cov) ) scale_cov else scale_cov[idx,]
				shape_cov_bs = if( is.null(shape_cov) ) shape_cov else shape_cov[idx,]
				fscale_bs    = if( is.null(fscale) || length(fscale) == 1 ) fscale    else fscale[idx]
				fshape_bs    = if( is.null(fshape) || length(fshape) == 1 ) fshape    else fshape[idx]
				
				private$fit_( Y[idx] , loc_bs , scale_cov_bs , shape_cov_bs , fscale_bs , fshape_bs )
				self$coefs_bootstrap = base::rbind( self$coefs_bootstrap , self$coef_ )
			}
			self$confidence_interval = base::apply( self$coefs_bootstrap , 2 , stats::quantile , probs = base::c( self$alpha / 2. , 1. - self$alpha / 2. ) )
		}
		
		private$fit_( Y , loc , scale_cov , shape_cov , fscale , fshape )
	}
	##}}}
	
	),
	
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	
	
	#############
	## Methods ##
	#############
	
	fit_ = function( Y , loc , scale_cov , shape_cov , fscale , fshape ) ##{{{
	{
		private$Y_ = as.vector(Y)
		self$loc = if( length(loc) == private$size_) as.vector(loc) else as.vector( base::rep( loc[1] , private$size_ ) )
		self$scale_$init( X = scale_cov , fix_values = fscale , size = private$size_ )
		self$shape_$init( X = shape_cov , fix_values = fshape , size = private$size_ )
		
		if( self$method == "moments" )
		{
			private$fit_moments()
		}
		else if( self$method == "lmoments" )
		{
			private$fit_lmoments()
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
		idx_excess = (private$Y_ > self$loc)
		excess = private$Y_[idx_excess] - self$loc[idx_excess]
		
		S = np_std( excess , self$scale_$design_wo1()[idx_excess,] , linkFct = self$scale_$linkFct , return_coef = TRUE )
		self$scale_$set_coef( S )
		self$shape_$set_intercept( self$shape_$linkFct$inverse(-1e-8) )
		
		
		self$scale_$update()
		self$shape_$update()
		
		self$scale = self$scale_$valueLf()
		self$shape = self$shape_$valueLf()
		
	},
	##}}}
	
	fit_lmoments = function() ##{{{
	{
		idx_excess = (private$Y_ > self$loc)
		excess   = private$Y_[idx_excess] - self$loc[idx_excess]
		lmo1     = SDFC::np_lmoments(excess,1)
		lmo2     = SDFC::np_lmoments(excess,2)
		itau     = lmo1 / lmo2
		scale_lm = lmo1 * ( itau - 1 )
		scale_lm = if( scale_lm > 0 ) scale_lm else 1e-8
		shape_lm = - ( itau - 2 )
		
		if( self$scale_$not_fixed() )
		{
			self$scale_$set_intercept( self$scale_$linkFct$inverse(scale_lm) )
		}
		if( self$shape_$not_fixed() )
		{
			self$shape_$set_intercept( self$shape_$linkFct$inverse(shape_lm) )
		}
		
		self$scale_$update()
		self$shape_$update()
		
		self$scale = self$scale_$valueLf()
		self$shape = self$shape_$valueLf()
	},
	##}}}
	
	fit_mle = function()##{{{
	{
		private$fit_lmoments()
		param_init = private$concat_param()
		
		## Test for initial value
		nll  = private$optim_function(param_init)
		gnll = private$gradient_optim_function(param_init)
		
		if( !is.finite(nll) || !is.finite(gnll) )
		{
			self$shape_$set_coef( numeric( self$shape_$size_ ) )
			param_init = private$concat_param()
		}
		
		optim_result = stats::optim( param_init , fn = private$optim_function , gr = private$gradient_optim_function , method = "BFGS" )
		private$update_param( optim_result$par )
	
	},
	##}}}
	
	
	split_param = function( param )##{{{
	{
		param_scale = NULL
		param_shape = NULL
		s1 = self$scale_$size_
		s2 = self$shape_$size_
		
		if( self$scale_$not_fixed() && self$shape_$not_fixed() )
		{
			param_scale = param[1:s1]
			param_shape = param[(s1+1):(s1+s2)]
		}
		else if( self$scale_$not_fixed() )
		{
			param_scale = param
		}
		else if( self$shape_$not_fixed() )
		{
			param_shape = param
		}
		
		return( list( scale = param_scale , shape = param_shape ) )
	},
	##}}}
	
	concat_param = function()##{{{
	{
		param = NULL
		param_scale = if( self$scale_$not_fixed() ) self$scale_$coef_ else NULL
		param_shape = if( self$shape_$not_fixed() ) self$shape_$coef_ else NULL
		
		param = base::c( param_scale , param_shape )
		
		return( param )
	},
	##}}}
	
	update_param = function( param ) ##{{{
	{
		param_sp = private$split_param(param)
		self$scale_$set_coef( param_sp$scale )
		self$shape_$set_coef( param_sp$shape )
		self$scale_$update()
		self$shape_$update()
		self$scale = self$scale_$valueLf()
		self$shape = self$shape_$valueLf()
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
		idx_excess = (private$Y_ > self$loc)
		loc   = self$loc[idx_excess]
		scale = self$scale[idx_excess]
		shape = self$shape[idx_excess]
		Z = 1. + shape * ( private$Y_[idx_excess] - loc ) / scale
		
		if( base::any(!(Z > 0)) )
			return(Inf)
		
		res = base::sum( base::log( scale ) + base::log(Z) * ( 1 + 1. / shape ) )
		
		if( !is.finite(res) )
			return(Inf)
		return(res)
	},
	##}}}
	
	optim_function = function( param )##{{{
	{
		private$update_param(param)
		return( private$negloglikelihood() )
	},
	##}}}
	
	gradient_optim_function = function( param ) ##{{{
	{
		private$update_param(param)
		
		idx      = ( private$Y_ > self$loc )
		Y        = private$Y_[idx]
		loc      = self$loc[idx]
		scale    = self$scale[idx]
		shape    = self$shape[idx]
		Z        = ( Y - loc ) / scale
		ZZ       = 1. + shape * Z
		exponent = 1. + 1. / shape
		
		grad = base::c()
		
		if( self$scale_$not_fixed() )
		{
			gr_scale   = self$scale_$valueGrLf()[idx]
			grad_scale = base::t(self$scale_$design_[idx,]) %*% ( gr_scale * ( - exponent * shape * Z / ZZ / scale + 1. / scale ) )
			grad       = base::c( grad , grad_scale )
		}
		
		if( self$shape_$not_fixed() )
		{
			gr_shape   = self$shape_$valueGrLf()[idx]
			grad_shape = base::rep( NaN , self$shape_$size_ )
			if( base::all(ZZ > 0) )
				grad_shape = base::t(self$shape_$design_[idx,]) %*% ( gr_shape * ( - base::log(ZZ) / shape^2 + exponent * Z / ZZ ) )
			grad       = base::c( grad , grad_shape )
		}
		
		return( grad )
	}
	##}}}
	
	
	)
	
)
##}}}



