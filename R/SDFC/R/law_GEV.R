
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


#' GEV (Generalized Extreme Value distribution)
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
#'   \item{\code{new(method,n_bootstrap,alpha)}}{Initialize GEV law with code{GEV}}
#'   \item{\code{fit(Y,...)}}{Fit the GEV law}.
#' }
#' @examples
#' ## Start by generate non-stationary GEV dataset
#' size = 2000
#' c_data = dataset.covariates(size)
#' 
#' t       = c_data$t
#' X_loc   = c_data$X_loc
#' X_scale = c_data$X_scale
#' X_shape = c_data$X_shape
#' 
#' loc   = 1.  + 0.8  * X_loc
#' scale = 0.2 + 0.08 * X_scale
#' shape = 0.  + 0.3  * X_shape
#' 
#' 
#' Y = SDFC::rgev( size , loc , scale , shape )
#' 
#' ## Regression with MLE
#' law = SDFC::GEV$new( "mle" )
#' law$fit( Y , c_loc = X_loc , c_scale = X_scale , c_shape = X_shape )
#' 
#' ## Assuming scale is known (available for any covariates)
#' law = SDFC::GEV$new( "mle" )
#' law$fit( Y , c_loc = X_loc , f_scale = scale , c_shape = X_shape )
#' 
#' ## And if we want a link function
#' law = SDFC::GEV$new( "mle" )
#' law$fit( Y , c_loc = X_loc , c_scale = X_scale , l_scale = SDFC::ExpLink$new() , c_shape = X_shape )
#' 
#' ## If we do not give a parameter, it is assumed constant
#' law = SDFC::GEV$new( "mle" )
#' law$fit( Y , c_scale = X_scale )
#' 
#' @export
GEV = R6::R6Class( "GEV" ,
	
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
			self$params$set_intercept( pscale$link$inverse(iscale) , "scale" )
		
		if( !ploc$is_fix() )
		{
			if( pscale$is_fix() )
			{
				iloc = m - 0.57722 * base::exp(pscale$value)
				self$params$update_coef( np_mean( iloc , ploc$design_wo1() , value = FALSE , link = ploc$link ) , "loc" )
			}
			else
			{
				self$params$set_intercept( ploc$link$inverse(iloc) , "loc" )
			}
		}
		
		if( !pshape$is_fix() )
			self$params$set_intercept( pshape$link$inverse(ishape) , "shape" )
		
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
			self$params$set_intercept( pscale$link$inverse(iscale) , "scale" )
		
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
				self$params$set_intercept( ploc$link$inverse(iloc) , "loc" )
			}
		}
		
		## Fit shape
		if( !pshape$is_fix() )
			self$params$set_intercept( pshape$link$inverse(ishape) , "shape" )
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
	
	fit_quantiles = function()##{{{
	{
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		if( !ploc$is_fix() )
		{
			loc = np_quantile( private$Y , base::c(base::exp(-1)) , c_Y = ploc$design_wo1() , value = TRUE )
			self$params$update_coef( np_mean( loc , ploc$design_wo1() , link = ploc$link , value = FALSE ) , "loc" )
		}
		
		if( !pscale$is_fix() )
		{
			qscale = base::c( 0.25 , 0.5 , 0.75 )
			coef   = - 1. / base::log( - base::log(qscale) )
			qreg = matrix( np_quantile( private$Y - self$loc , qscale , pscale$design_wo1() ) , ncol = 3 )
			fscale = base::t( base::apply( qreg , 1 , function(x) { coef * x } ) )
			fscale = base::apply( fscale , 1 , base::mean )
			fscale[!(fscale > 0)] = 0.1
			self$params$update_coef( np_mean( fscale , pscale$design_wo1() , link = pscale$link , value = FALSE ) , "scale" )
		}
		
		if( !pshape$is_fix() )
		{
			p = base::c( 0.1 , 0.9 )
			qval = matrix( np_quantile( (private$Y - self$loc) / self$scale , p , pshape$design_wo1() ) , ncol = 2 )
			kappa = qval[,1] / qval[,2]
			llp = base::log( - base::log(p) )
			shape = ( 2 * (llp[1] - kappa * llp[2] ) / ( llp[1]^2 - kappa * llp[2]^2 ) )
			self$params$update_coef( np_mean( shape , pshape$design_wo1() , link = pshape$link , value = FALSE ) , "shape" )
		}
	},
	##}}}
	
	initialization_mle = function()##{{{
	{
		private$fit_lmoments_experimental()
		nlll = private$negloglikelihood(self$coef_)
		grad = private$gradient_nlll(self$coef_)
		
		f_scale = 1
		f_shape = 1
		while( !is.finite(nlll) || !base::any(base::is.finite(grad)) )
		{
			pscale = self$params$dparams_[["scale"]]
			pshape = self$params$dparams_[["shape"]]
			
			if( pshape$is_fix() && !pscale$is_fix() )
			{
				coef_ = base::rep( 0 , pscale$n_features )
				coef_[1] = pscale$link$inverse(1. * f_scale)
				self$params$update_coef( coef_ , "scale" )
			}
			else if( !pshape$is_fix() )
			{
				coef_ = base::rep( 0 , pshape$n_features )
				coef_[1] = pshape$link$inverse(1e-1 / f_shape)
				self$params$update_coef( coef_ , "shape" )
			}
			else
			{
				self$fit_quantiles()
			}
			
			f_scale = f_scale * 2
			f_shape = f_shape * 2
			nlll = private$negloglikelihood(self$coef_)
			grad = private$gradient_nlll(self$coef_)
		}
	},
	##}}}
	
	fit_ = function()##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
		else if( self$method == "lmoments" )
			private$fit_lmoments()
		else if( self$method == "quantiles" )
			private$fit_quantiles()
		else
			private$fit_lmoments_experimental()
	},
	##}}}
	
	negloglikelihood = function( coef )##{{{
	{
		self$coef_ = coef
		## Impossible scale
		if( !base::all( self$scale > 0 ) )
			return(Inf)
		
		## Remove exponential case
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


