
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

## NormalLaw {{{

#' NormalLaw (Gaussian Law)
#'
#' Class to fit a Normal law.
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param method [string]
#'        Fit method, "moments" and "MLE" are available.
#' @param link_fct_loc [SDFC::LinkFct]
#'        Link function for loc parameter. Can be an element of SDFC, or a class based on SDFC::LinkFct
#' @param link_fct_scale [SDFC::LinkFct]
#'        Link function for scale parameter. Can be an element of SDFC, or a class based on SDFC::LinkFct
#' @param n_bootstrap [int]
#'        Number of bootstrap, default 0
#' @param alpha [float]
#'        Level of confidence interval, default 0.05
#' @param loc_cov  [matrix or NULL ]
#'        Location covariate for fit
#' @param scale_cov  [matrix or NULL]
#'        Scale covariate for fit
#' @param floc [vector or NULL]
#'        Value of loc if it is not necessary to fit
#' @param fscale [vector or NULL]
#'        Value of scale if it is not necessary to fit
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(method,link_fct_loc,link_fct_scale,n_bootstrap,alpha)}}{Initialize Normal law with code{NormalLaw}}
#'   \item{\code{fit(Y,loc_cov,scale_cov,floc,fscale)}}{Fit the Normal law}.
#' }
#' @examples
#' ## Data
#' size = 2500
#' t    = base::seq( 0 , 1 , length = size )
#' X0    = t^2
#' X1    = base::cos( 2 * base::pi * t )
#' loc   = 1. + 2 * X0
#' scale = 0.6 + 0.5 * X1
#' Y    = stats::rnorm( n = size , mean = loc , sd = scale )
#' 
#' 
#' ## Fit
#' law = SDFC::NormalLaw$new( method = "MLE" ,  n_bootstrap = 10 )
#' law$fit( Y , loc_cov = X0 , scale_cov = X1 )
#' law$loc   ## Loc fitted
#' law$scale ## Scale fitted
#' @export
NormalLaw = R6::R6Class( "NormalLaw" ,
	
	inherit = AbstractLaw,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	loc      = NULL,
	scale    = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( method = "MLE" , link_fct_loc = IdLinkFct$new() , link_fct_scale = IdLinkFct$new() , n_bootstrap = 0 , alpha = 0.05 ) ##{{{
	{
		super$initialize( method , n_bootstrap , alpha )
		self$loc       = NULL
		self$scale     = NULL
		private$loc_   = LawParam$new( linkFct = link_fct_loc   , kind = "loc"   )
		private$scale_ = LawParam$new( linkFct = link_fct_scale , kind = "scale" )
	},
	##}}}
	
	
	#############
	## Methods ##
	#############
	
	fit = function( Y , loc_cov = NULL , scale_cov = NULL , floc = NULL , fscale = NULL )##{{{
	{
		Y = as.vector(Y)
		private$size_ = length(Y)
		
		## Bootstrap
		if( self$n_bootstrap > 0 )
		{
			if( !is.null(loc_cov) && !is.matrix(loc_cov) )
				loc_cov = matrix( loc_cov , nrow = private$size_ , ncol = 1 )
			if( !is.null(scale_cov) && !is.matrix(scale_cov) )
				scale_cov = matrix( scale_cov , nrow = private$size_ , ncol = 1 )
			
			self$coefs_bootstrap = base::c()
			
			for( i in 1:self$n_bootstrap )
			{
				idx = base::sample( 1:private$size_ , private$size_ , replace = TRUE )
				loc_cov_bs   = if( is.null(loc_cov) )   loc_cov   else loc_cov[idx,]
				scale_cov_bs = if( is.null(scale_cov) ) scale_cov else scale_cov[idx,]
				floc_bs      = if( is.null(floc) || length(floc) == 1 )     floc      else floc[idx]
				fscale_bs    = if( is.null(fscale) || length(fscale) == 1 ) fscale    else fscale[idx]
				
				private$fit_( Y[idx] , loc_cov_bs , scale_cov_bs , floc_bs , fscale_bs )
				self$coefs_bootstrap = base::rbind( self$coefs_bootstrap , self$coef_ )
			}
			self$confidence_interval = base::apply( self$coefs_bootstrap , 2 , stats::quantile , probs = base::c( self$alpha / 2. , 1. - self$alpha / 2. ) )
		}
		
		
		## Good fit
		private$fit_( Y , loc_cov , scale_cov , floc , fscale )
	}
	##}}}
	
	
	),
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	loc_     = NULL,
	scale_   = NULL,
	
	
	#############
	## Methods ##
	#############
	
	fit_ = function( Y , loc_cov , scale_cov , floc , fscale ) ##{{{
	{
		private$Y_    = as.vector(Y)
		
		private$loc_$init(   X = loc_cov   , fix_values = floc   , size = private$size_ )
		private$scale_$init( X = scale_cov , fix_values = fscale , size = private$size_ )
		
		if( self$method == "moments" )
		{
			private$fit_moments()
		}
		else
		{
			private$fit_mle()
		}
		
		self$coef_ = private$concat_param()
	},
	##}}}
	
	fit_moments = function()##{{{
	{
		if( private$loc_$not_fixed() )
		{
			lX = private$loc_$design_wo1()
			private$loc_$set_coef( np_mean( private$Y_ , lX , return_coef = TRUE , linkFct = private$loc_$linkFct ) )
			private$loc_$update()
		}
		
		if( private$scale_$not_fixed() )
		{
			sX = private$scale_$design_wo1()
			private$scale_$set_coef( np_std( private$Y_ , sX , m = private$loc_$valueLf() , return_coef = TRUE , linkFct = private$scale_$linkFct ) )
			private$scale_$update()
		}
		self$loc   = private$loc_$valueLf()
		self$scale = private$scale_$valueLf()
	},
	##}}}
	
	fit_mle = function()##{{{
	{
		private$fit_moments()
		param_init = private$concat_param()
		optim_result = stats::optim( param_init , fn = private$optim_function , gr = private$gradient_optim_function , method = "BFGS" )
		private$update_param( optim_result$par )
	},
	##}}}
	
	split_param = function( param )##{{{
	{
		param_loc   = NULL
		param_scale = NULL
		
		if( private$loc_$not_fixed() )
		{
			param_loc   = param[1:private$loc_$size_]
			if( private$scale_$not_fixed() )
			{
				param_scale = param[(private$loc_$size_+1):length(param)]
			}
		}
		else if( private$scale_$not_fixed() )
		{
			param_scale = param[1:private$scale_$size_]
		}
		
		return( list( loc = param_loc , scale = param_scale ) )
	},
	##}}}
	
	concat_param = function()##{{{
	{
		param = NULL
		if( private$loc_$not_fixed() && private$scale_$not_fixed() )
		{
			param = base::c( private$loc_$coef_ , private$scale_$coef_ )
		}
		else if( private$loc_$not_fixed() )
		{
			param = private$loc_$coef_
		}
		else if( private$scale_$not_fixed() )
		{
			param = private$scale_$coef_
		}
		return( param )
	},
	##}}}
	
	negloglikelihood = function() ##{{{
	{
		if( !base::all( self$scale > 0 ) )
		{
			return(Inf)
		}
		
		scale2 = self$scale^2
		res =  base::sum( base::log( scale2 ) ) / 2. + base::sum( ( private$Y_ - self$loc )^2 / scale2 ) / 2. 
		return(res)
	},
	##}}}
	
	update_param = function( param ) ##{{{
	{
		param_sp = private$split_param(param)
		private$loc_$set_coef(   param_sp$loc )
		private$scale_$set_coef( param_sp$scale )
		private$loc_$update()
		private$scale_$update()
		self$loc   = private$loc_$valueLf()
		self$scale = private$scale_$valueLf()
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
		
		Yc = private$Y_ - self$loc
		grad = base::c()
		
		if( private$loc_$not_fixed() )
		{
			grad_loc = - base::t( private$loc_$design_  ) %*% ( Yc / self$scale^2 * private$loc_$valueGrLf() )
			grad     = base::c( grad , grad_loc )
		}
		if( private$scale_$not_fixed() )
		{
			grad_scale = base::t( private$scale_$design_) %*% ( ( 1. / self$scale - Yc^2 / self$scale^3 ) * private$scale_$valueGrLf() )
			grad       = base::c( grad , grad_scale )
		}
		
		return( grad )
	}
	##}}}
	
	
	),
	
)
##}}}

