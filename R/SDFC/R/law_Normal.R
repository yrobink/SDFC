
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

## Normal {{{

#' Normal (Gaussian Law)
#'
#' Class to fit a Normal law.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param method [string]
#'        Fit method, "moments", "MLE" and "Bayesian" are available.
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
Normal = R6::R6Class( "Normal" ,##{{{
	
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
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		if( !ploc$is_fix() )
			self$params$update_coef( np_mean( private$Y , ploc$design_wo1() , value = FALSE , link = ploc$link ) , "loc" )
		
		if( !pscale$is_fix() )
			self$params$update_coef( np_std( private$Y , pscale$design_wo1() , m_Y = self$loc , value = FALSE , link = pscale$link ) , "scale" )
	},
	##}}}
	
	initialization_mle = function()##{{{
	{
		private$fit_moments()
	},
	##}}}
	
	fit_ = function()##{{{
	{
		if( self$method == "moments" )
			private$fit_moments()
	},
	##}}}
	
	negloglikelihood = function( coef )##{{{
	{
		self$coef_ = coef
		if( !base::all( self$scale > 0 ) )
		{
			return(Inf)
		}
		
		scale2 = self$scale^2
		res =  base::sum( base::log( scale2 ) ) / 2. + base::sum( ( private$Y - self$loc )^2 / scale2 ) / 2. 
		return(res)
	},
	##}}}
	
	gradient_nlll = function( coef ) ##{{{
	{
		self$coef_ = coef
		ploc   = self$params$dparams_[["loc"]]
		pscale = self$params$dparams_[["scale"]]
		
		Yc = private$Y - self$loc
		grad = base::c()
		
		if( !ploc$is_fix() )
		{
			grad_loc = - base::t( ploc$design_  ) %*% ( Yc / self$scale^2 * ploc$gradient() )
			grad     = base::c( grad , grad_loc )
		}
		if( !pscale$is_fix() )
		{
			grad_scale = base::t( pscale$design_) %*% ( ( 1. / self$scale - Yc^2 / self$scale^3 ) * pscale$gradient() )
			grad       = base::c( grad , grad_scale )
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
		super$initialize( base::c( "loc" , "scale" ) , method , n_bootstrap , alpha )
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
	}
	##}}}
	
	)
	##}}}
	#################

)
##}}}
##}}}

