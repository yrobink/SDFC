
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

#' Gamma distribution
#'
#' Class to fit a Gamma law.
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
#'   \item{\code{new(method,n_bootstrap,alpha)}}{Initialize Gamma law with code{Gamma}}
#'   \item{\code{fit(Y,...)}}{Fit the Gamma law}.
#' }
#' @examples
#' ## Start by generate non-stationary Gamma dataset
#' size = 2000
#' c_data = dataset.covariates(size)
#' 
#' t       = c_data$t
#' X_scale = c_data$X_scale
#' X_shape = c_data$X_shape
#' 
#' scale = 0.2 + 0.08 * X_scale
#' shape = 0.4 + 0.3  * X_shape
#' 
#' 
#' Y = rgamma( size , shape = shape , scale = scale )
#' 
#' ## Regression with MLE
#' law = SDFC::Gamma$new( "mle" )
#' law$fit( Y , c_scale = X_scale , c_shape = X_shape )
#' 
#' ## Assuming scale is known (available for any covariates)
#' law = SDFC::Gamma$new( "mle" )
#' law$fit( Y , f_scale = scale , c_shape = X_shape )
#' 
#' ## And if we want a link function
#' law = SDFC::Gamma$new( "mle" )
#' law$fit( Y , c_scale = X_scale , l_scale = SDFC::ExpLink$new() , c_shape = X_shape )
#' 
#' ## If we do not give a parameter, it is assumed constant
#' law = SDFC::Gamma$new( "mle" )
#' law$fit( Y , c_scale = X_scale )
#' 
#' @export
Gamma = R6::R6Class( "Gamma" ,
	
	inherit = AbstractLaw,
	
	##################
	## Private list ##
	##{{{
	
	private = list(
	
	## Arguments
	##==========
	
	## Methods
	##========
	
	fit_moments = function()##{{{
	{
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		
		if( !pscale$is_fix() && !pshape$is_fix() )
		{
			mX = numeric( self$params$n_samples ) + 1
			vX = numeric( self$params$n_samples ) + 1
			for( i in 2:pscale$n_features )
			{
				for( j in 1:pshape$n_features )
				{
					mX = base::cbind( mX , pscale$design_[,i]   * pshape$design_[,j] )
					vX = base::cbind( vX , pscale$design_[,i]^2 * pshape$design_[,j] )
				}
			}
			
			m   = np_mean( private$Y , mX[,2:base::ncol(mX)] )
			v   = np_var(  private$Y , vX[,2:base::ncol(vX)] )
			idx = ( abs(m) < 1e-8 ) | v < 1e-8
			
			scale = numeric(self$params$n_samples)
			shape = numeric(self$params$n_samples)
			scale[!idx] = v[!idx] / m[!idx]
			shape[!idx] = m[!idx]^2 / v[!idx]
			
			if( base::any(idx) )
			{
				scale[idx] = base::min(scale[!idx])
				shape[idx] = base::min(shape[!idx])
			}
			self$params$update_coef( np_mean( scale , pscale$design_wo1() , value = FALSE , link = pscale$link ) , "scale" )
			self$params$update_coef( np_mean( shape , pshape$design_wo1() , value = FALSE , link = pshape$link ) , "shape" )
		}
	},
	##}}}
	
	initialization_mle = function()##{{{
	{
		private$fit_moments()
	},
	##}}}
	
	fit_ = function()##{{{
	{
		private$fit_moments()
	},
	##}}}
	
	negloglikelihood = function( coef )##{{{
	{
		self$coef_ = coef
		## Impossible scale
		if( !base::all( self$scale > 0 ) || !base::all( self$shape > 0 ) || !base::all( private$Y > 0 ) )
			return(Inf)
		
		res = base::sum( private$Y / self$scale + base::lgamma(self$shape) + self$shape * base::log(self$scale) - (self$shape-1) * base::log(private$Y) )
		
		if( is.finite(res) )
			return(res)
		else
			return(Inf)
	},
	##}}}
	
	gradient_nlll = function( coef ) ##{{{
	{
		self$coef_ = coef
		
		if( !base::all( self$scale > 0 ) || !base::all( self$shape > 0 ) || !base::all( private$Y > 0 ) )
		{
			return( base::rep( NaN , length(coef) ) )
		}
		
		pscale = self$params$dparams_[["scale"]]
		pshape = self$params$dparams_[["shape"]]
		
		
		grad = base::c()
		if( !pscale$is_fix() )
		{
			scale_vect = pscale$gradient() * ( self$shape / self$scale - private$Y / self$shape^2 )
			grad_scale = base::t(pscale$design_) %*% scale_vect
			grad       = base::c( grad , grad_scale )
		}
		if( !pshape$is_fix() )
		{
			shape_vect = pshape$gradient() * ( base::digamma(self$shape) + base::log(self$scale) - base::log(private$Y) )
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
		super$initialize( base::c( "scale" , "shape" ) , method , n_bootstrap , alpha )
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




























