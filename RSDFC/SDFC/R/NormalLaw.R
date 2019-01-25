
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


#' NormalLaw (Gaussian)
#'
#' Class to fit a Normal law with covariates
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param loc  [vector]
#'        Location parameters
#' @param scale  [vector]
#'        Scale parameters
#' @param loc_cov  [matrix]
#'        Location covariate for fit
#' @param scale_cov  [matrix]
#'        Scale covariate for fit
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
#'   \item{\code{new(use_phi,method,verbose)}}{Initialize normal law with code{NormalLaw}}
#'   \item{\code{fit(Y,loc_cov,scale_cov)}}{Fit the Normal law}.
#' }
#' @examples
#' ## Data
#' size = 2000
#' data = SDFC::Dataset2(size)
#' t = data$t
#' X = data$X
#' Y = data$Y
#' 
#' ## Normal Law
#' norm = SDFC::NormalLaw$new()
#' norm$fit( Y , loc_cov = X , scale_cov = X )
#' print(norm$loc_coef_) ## Location coef fitted
#' print(norm$scale_coef_) ## Scale coef fitted
#'
#' ## In fact, it is better to fit:
#' norm_best = SDFC::NormalLaw$new()
#' norm_best$fit( Y , loc_cov = X[,1] , scale_cov = X[,2] )
#' @export
NormalLaw = R6::R6Class( "NormalLaw" ,
	
	public = list(
	
	use_phi      = NULL,
	method       = NULL,
	verbose      = NULL,
	Y            = NULL,
	size         = NULL,
	Nlog2pi      = NULL,
	loc          = NULL,
	scale        = NULL,
	loc_design   = NULL,
	scale_design = NULL,
	nloc         = NULL,
	nscale       = NULL,
	ncov         = NULL,
	optim_result = NULL,
	loc_coef_    = NULL,
	scale_coef_  = NULL,
	
	
	###############
	## Arguments ##
	###############
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( use_phi = FALSE , method = "BFGS" , verbose = FALSE ) ##{{{
	{
		self$use_phi = use_phi
		self$method  = method
		self$verbose = verbose
	},
	##}}}

	fit = function( Y , loc_cov = NULL , scale_cov = NULL )##{{{
	{
		self$Y    = Y
		self$size = base::length(Y)
		self$Nlog2pi = self$size * base::log( 2 * base::pi ) / 2.
		
		## Design matrix
		self$loc_design   = if( !is.null(loc_cov)   ) base::cbind( 1 , loc_cov)   else matrix( 1 , nrow = self$size , ncol = 1 )
		self$scale_design = if( !is.null(scale_cov) ) base::cbind( 1 , scale_cov) else matrix( 1 , nrow = self$size , ncol = 1 )
		
		if( base::qr(self$loc_design)$rank < base::ncol(self$loc_design) )
		{
			if( self$verbose )
			{
				print( "SFDC::NormalLaw: singular design matrix for loc, co-variable coefficients are set to 0" )
			}
			self$loc_design = matrix( 1 , nrow = self$size , ncol = 1 )
		}
		if( base::qr(self$scale_design)$rank < base::ncol(self$scale_design) )
		{
			if( self$verbose )
			{
				print( "SFDC::NormalLaw: singular design matrix for scale, co-variable coefficients are set to 0" )
			}
			self$scale_design = matrix( 1 , nrow = self$size , ncol = 1 )
		}
		
		self$nloc   = base::ncol(self$loc_design)
		self$nscale = base::ncol(self$scale_design)
		self$ncov   = self$nloc + self$nscale
		
		## Initial condition
		param_init = self$find_init()
		
		## Optimization
		self$optim_result = stats::optim( param_init , fn = self$optim_function , gr = self$gradient_optim_function , method = self$method )
		
		## Set result
		self$loc_coef_   = self$optim_result$par[1:self$nloc]
		self$scale_coef_ = self$optim_result$par[(self$nloc+1):self$ncov]
		self$update_param( self$optim_result$par )
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
		self$loc   = self$loc_design %*% param[1:self$nloc]
		self$scale = self$link( self$scale_design %*% param[(self$nloc+1):self$ncov] )
	},
	##}}}
	
	find_init = function() ##{{{
	{
		## Estimate mu
		lm_res = stats::lm( self$loc_design ~ self$Y )
		init_loc = lm_res$coefficients[1:self$nloc]
		init_scale = numeric(self$nscale)
		init_scale[1] = stats::sd(self$Y)
		return( base::c(init_loc,init_scale) )
	},
	##}}}
	
	negloglikelihood = function() ##{{{
	{
		if( !base::all( self$scale > 0 ) )
		{
			return(Inf)
		}
		
		scale2 = self$scale^2
		return( self$Nlog2pi + base::sum( base::log( scale2 ) ) / 2. + base::sum( ( self$Y - self$loc )^2 / scale2 ) / 2. )
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
		
		Yc = self$Y - self$loc
		grad_loc   = - base::t(self$loc_design) %*% (Yc / self$scale^2)
		grad_phi   = if( self$use_phi ) 1. else self$scale
		grad_scale = base::t(self$scale_design) %*% ( 1. / grad_phi - (Yc / self$scale)^2 / grad_phi )
		
		return( base::c(grad_loc,grad_scale) )
	}
	##}}}

	)
)
