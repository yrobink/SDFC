
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


#' QuantileRegression
#'
#' Class to perform a Quantile Regression with covariates
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param ltau [vector]
#'        Vector of quantiles where we want to perform the fit
#' @param method [string]
#'        Method used to fit, currently, only "Frish-Newton" is available
#' @param verbose [bool]
#'        Print warning and error message
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(ltau,method,verbose)}}{Initialize Quantile Regression with code{QuantileRegression}}
#'   \item{\code{fit(Y,X)}}{Fit the quantile regression}.
#'   \item{\code{predict()}}{Return the quantile fitted}.
#'   \item{\code{coef()}}{Return coefficients fitted}.
#'   \item{\code{is_fitted()}}{TRUE if fit is already called}.
#'   \item{\code{is_success()}}{TRUE if fit is a success}.
#'   \item{\code{is_unfeasible()}}{TRUE if fit is not feasible}.
#' }
#' @examples
#' ## Data
#' size = 2000
#' c_data = dataset.covariates(size)
#' 
#' loc   = 0.5 + 2 * c_data$X_loc
#' scale = 1 + 2 * c_data$X_scale
#' Y = stats::rnorm( size , mean = loc , sd = scale )
#' 
#' ## Quantile regression
#' ltau  = base::seq( 0.01 , 0.99 , length = 100 )
#' qr = SDFC::QuantileRegression$new(ltau)
#' qr$fit( Y , base::cbind( c_data$X_loc , c_data$X_scale ) )
#' Yq = if( qr$is_success() ) qr$predict() else NULL
#' @export
QuantileRegression = R6::R6Class( "QuantileRegression" ,
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	ltau    = NULL,
	maxit   = NULL,
	tol     = NULL,
	beta    = NULL,
	method  = NULL,
	verbose = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( ltau , method = "Frish-Newton" , verbose = FALSE ) ##{{{
	{
		self$ltau        = ltau
		self$maxit       = 50
		self$tol         = 1e-6
		self$beta        = 0.99995
		self$method      = method
		self$verbose     = verbose
		private$qrmethod = FrishNewton$new( ltau , self$maxit , self$tol , self$beta )
	},
	##}}}
	
	
	###############
	## Accessors ##
	###############
	
	coef = function() ##{{{
	{
		return( private$qrmethod$coef() )
	},
	##}}}
	
	
	###########
	## State ##
	###########
	
	is_fitted = function() ##{{{
	{
		return( private$qrmethod$is_fitted() )
	},
	##}}}
	
	is_success = function() ##{{{
	{
		return( private$qrmethod$is_success() )
	},
	##}}}
	
	is_unfeasible = function() ##{{{
	{
		return( private$qrmethod$is_unfeasible() )
	},
	##}}}
	
	#############
	## Methods ##
	#############
	
	fit = function( Y , X ) ##{{{
	{
		if( !is.matrix(X) )
		{
			X = matrix( X , nrow = length(X) , ncol = 1 ) 
		}
		private$qrmethod$fit( Y , X )
		if( private$qrmethod$is_unfeasible() && self$verbose )
		{
			print( "SDFC::QuantileRegression : Unfeasible problem" )
		}
	},
	##}}}
	
	predict = function() ##{{{
	{
		return( private$qrmethod$predict() )
	}
	##}}}
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	qrmethod = NULL
	
	
	#############
	## Methods ##
	#############
	
	
	)
)

