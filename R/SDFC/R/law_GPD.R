
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



#' @export
GPD = R6::R6Class( "GPD" ,
	
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
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		idx = (private$Y > self$loc)
		Y   = (private$Y[idx] - self$loc[idx])
		
		
		if( !pscale$is_fix() )
		{
			c_scale = pscale$design_wo1()
			if( !is.null(c_scale) )
				c_scale = c_scale[idx,]
			
			self$params$update_coef( n_std( Y , c_Y = c_scale , m_Y = 0 , value = FALSE , link = pscale$link ) , "scale" )
		}
		
		if( !pshape$is_fix() )
			self$params$set_intercept( -1e-8 , "shape" )
		
	},
	##}}}
	
	fit_lmoments = function() ##{{{
	{
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		idx = (private$Y > self$loc)
		Y   = (private$Y[idx] - self$loc[idx])
		lmom = np_lmoments( Y )
		
		if( !pscale$is_fix() && !pshape$is_fix() )
		{
			itau     = lmom[1] / lmom[2]
			scale_lm = lmom[1] * ( itau - 1 )
			scale_lm = if( scale_lm > 0 ) scale_lm else 1e-8
			shape_lm = 2 - itau
			self$params$set_intercept( scale_lm , "scale" )
			self$params$set_intercept( shape_lm , "shape" )
		}
		else if( !pscale$is_fix() )
		{
			scale = lmom[1] * ( 1 - self$shape )
			scale[ !(scale > 0) ] = 1e-8
			self$params$update_coef( np_mean( scale , pscale$design_wo1() , value = FALSE , link = pscale$link ) , "scale" )
		}
		else if( !pshape$is_fix() )
		{
			Y = Y / self$scale[idx]
			lmom = lmoments(Y)
			itau     = lmom[1] / lmom[2]
			self$params$set_intercept( 2 - itau , "shape" )
		}
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
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		idx = (private$Y > self$loc)
		Y   = (private$Y[idx] - self$loc[idx])
		c_Y = matrix( c_Y[idx,] , ncol = base::ncol(c_Y) )
		
		## First step, find lmoments
		lmom = np_lmoments( Y , c_Y = c_Y )
		
		if( !pscale$is_fix() && !pshape$is_fix() )
		{
			itau  = lmom[,1] / lmom[,2]
			scale = lmom[,1] * ( itau - 1 )
			shape = 2 - itau
			scale_design = pscale$design_wo1()
			if( !is.null(scale_design) )
				scale_design = matrix( scale_design[idx,] , ncol = pscale$n_features - 1 )
			self$params$update_coef( np_mean( scale , scale_design , link = pscale$link , value = FALSE ) , "scale" )
			
			shape_design = pshape$design_wo1()
			if( !is.null(shape_design) )
				shape_design = matrix( shape_design[idx,] , ncol = pshape$n_features - 1 )
			self$params$update_coef( np_mean( shape , shape_design , link = pshape$link , value = FALSE ) , "shape" )
		}
		else if( !pscale$is_fix() )
		{
			scale = lmom[,1] * ( 1 - self$shape[idx] )
			scale_design = pscale$design_wo1()
			if( !is.null(scale_design) )
				scale_design = matrix( scale_design[idx,] , ncol = pscale$n_features - 1 )
			self$params$update_coef( np_mean( scale , scale_design , link = pscale$link , value = FALSE ) , "scale" )
		}
		else if( !pshape$is_fix() )
		{
			Y     = Y / self$scale[idx]
			lmom  = np_lmoments( Y , pshape$design_wo1()[idx,] )
			shape = 2 - lmom[,1] / lmom[,2]
			shape_design = pshape$design_wo1()
			if( !is.null(shape_design) )
				shape_design = matrix( shape_design[idx,] , ncol = pshape$n_features - 1 )
			self$params$update_coef( np_mean( shape , shape_design , link = pshape$link , value = FALSE ) , "shape" )
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
		idx_excess = (private$Y > self$loc)
		loc   = self$loc[idx_excess]
		scale = self$scale[idx_excess]
		shape = self$shape[idx_excess]
		Z = 1. + shape * ( private$Y[idx_excess] - loc ) / scale
		
		if( !base::all(Z > 0) )
			return(Inf)
		
		res = base::sum( base::log( scale ) + base::log(Z) * ( 1 + 1. / shape ) )
		
		if( is.finite(res) )
			return(res)
		else
			return(Inf)
	},
	##}}}
	
	gradient_nlll = function( coef ) ##{{{
	{
		self$coef_ = coef
		
		## Remove exponential case
		shape = self$shape
		zero_shape = ( base::abs(shape) < 1e-10 )
		if( base::any(zero_shape) )
			shape[zero_shape] = -1e-10
		
		
		## Usefull values
		idx      = ( private$Y > self$loc )
		Y        = private$Y[idx]
		loc      = self$loc[idx]
		scale    = self$scale[idx]
		shape    = shape[idx]
		Z        = ( Y - loc ) / scale
		ZZ       = 1. + shape * Z
		exponent = 1. + 1. / shape
		
		
		grad = base::c()
		
		pscale = self$params$dparams_[["scale"]]
		if( !pscale$is_fix() )
		{
			scale_vect = pscale$gradient() * ( - exponent * shape * Z / ZZ / scale + 1. / scale )
			grad_scale = base::t(pscale$design_[idx,]) %*% scale_vect
			grad       = base::c( grad , grad_scale )
		}
		
		pshape = self$params$dparams_[["shape"]]
		if( !pshape$is_fix() )
		{
			grad_shape = base::rep( NaN , pshape$n_features )
			if( base::all(ZZ > 0) )
			{
				shape_vect = pshape$gradient() * ( - base::log(ZZ) / shape^2 + exponent * Z / ZZ )
				grad_shape = base::t(pshape$design_[idx,]) %*% shape_vect
			}
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


