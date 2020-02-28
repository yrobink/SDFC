
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

#' AbstractLaw (Base class)
#'
#' Base class inherited by other laws, do not use it
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param kinds_params [vector] Vector of kind of params from herited law
#' @param method [string] Fit method
#' @param n_bootstrap [integer] Number of boostrap for confidence interval
#' @param alpha [float] Level of confidence interval
#'
#' @return Object of \code{\link{R6Class}} 
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(method,n_bootstrap,alpha)}}{Initialize abstract law with code{AbstractLaw}}
#' }
#' @examples
#' ## Data
#' @export
AbstractLaw = R6::R6Class( "AbstractLaw" ,##{{{

	##################
	## Private list ##
	##{{{
	
	private = list(
	
	## Arguments
	##==========
	
	kinds_params_ = NULL,
	method_       = NULL,
	Y             = NULL,
	
	
	## Methods
	##========
	
	fit_global = function( Y , ... )##{{{
	{
		self$params = LawParams$new( kinds = self$kinds_params )
		kwargs = list(...)
		kwargs$n_samples = length(Y)
		
		base::do.call( self$params$add_params , kwargs )
		private$Y = matrix( Y , ncol = 1 )
		
		if( self$method == "mle" )
			private$fit_mle()
		else if( self$method == "bayesian" )
			base::do.call( private$fit_bayesian , kwargs )
		else
			private$fit_()
		
		private$Y = NULL
	},
	##}}}
	
	fit_mle = function()##{{{
	{
		private$initialization_mle()
		self$info$optim_result = stats::optim( self$coef_ , fn = private$negloglikelihood , gr = private$gradient_nlll , method = "BFGS" , hessian = TRUE )
		self$coef_ = self$info$optim_result$par
		self$info$cov = base::solve(self$info$optim_result$hessian)
	},
	##}}}
	
	fit_bayesian = function(...)##{{{
	{
		kwargs = list(...)
		
		## Find numbers of features
		##=========================
		n_features = 0
		for( k in self$params$dparams_ )
		{
			n_features = n_features + k$n_features
		}
		
		## Define prior, default is multivariate normal law
		##=================================================
		prior = kwargs$prior
		if( is.null(prior) )
		{
			prior = list()
			prior$rvs    = function( n = 1 ) { return( matrix( stats::rnorm( n * n_features , mean = 0 , sd = 10 ) , nrow = n ) ) }
			prior$logpdf = function(x) { base::sum(base::log(stats::dnorm( x , mean = 0 , sd = 10 ))) }
		}
		
		## Define transition
		##==================
		transition = kwargs$transition
		if( is.null(transition) )
		{
			transition = function(x) { return( x + stats::rnorm( n = n_features , mean = 0 , sd = 0.1 ) ) }
		}
		
		## Define numbers of iterations of MCMC algorithm
		##===============================================
		n_mcmc_drawn = kwargs$n_mcmc_drawn
		if( is.null(n_mcmc_drawn) )
		{
			n_mcmc_drawn = 10000
		}
		
		## MCMC algorithm
		##===============
		draw   = matrix( NA , nrow = n_mcmc_drawn , ncol = n_features )
		accept = base::rep( TRUE , n_mcmc_drawn )
		
		## Init values
		##============
		init = kwargs$mcmc_init
		if( is.null(init) )
		{
			init = prior$rvs()
		}
		
		repeat
		{
			lll_current   = - private$negloglikelihood(init)
			prior_current = base::sum( prior$logpdf(init) )
			p_current     = prior_current + lll_current
			
			if( is.finite(p_current) )
				break
			init = prior$rvs()
		}
		draw[1,]  = init
		
		## Main loop on MCMC algorithm
		##============================
		for( i in 2:n_mcmc_drawn )
		{
			draw[i,] = transition(draw[i-1,])
			
			## Likelihood and probability of new points
			lll_next   = - private$negloglikelihood(draw[i,])
			prior_next = base::sum(prior$logpdf(draw[i,]))
			p_next     = prior_next + lll_next
			
			## Accept or not ?
			p_accept = base::exp( p_next - p_current )
			if( stats::runif(1) < p_accept )
			{
				lll_current   = lll_next
				prior_current = prior_next
				p_current     = p_next
				accept[i] = TRUE
			}
			else
			{
				draw[i,] = draw[i-1,]
				accept[i] = FALSE
			}
		}
		
		
		## Exclude 10% of first elements to find a "best estimate"
		##========================================================
		start = as.integer( 0.1 * n_mcmc_drawn )
		self$coef_ = base::apply( draw[start:n_mcmc_drawn,] , 2 , base::mean )
		
		## Update information
		##===================
		self$info$draw         = draw
		self$info$accept       = accept
		self$info$n_mcmc_drawn = n_mcmc_drawn
		self$info$rate_accept  = base::sum(accept) / n_mcmc_drawn
		self$info$cov          = stats::cov(draw)
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
	
	params    = NULL,
	info      = NULL,
	bootstrap = NULL,
	
	## Constructor
	##============
	initialize = function( kinds_params , method , n_bootstrap , alpha )##{{{
	{
		private$kinds_params_ = kinds_params
		self$method           = method
		self$bootstrap = list( n_bootstrap = n_bootstrap , alpha = alpha )
	},
	##}}}
	
	## Methods
	##========
	
#	print = function(...)
#	{
#		invisible(self)
#	},
	
	fit = function( Y , ... ) ##{{{
	{
		kwargs = list(...)
		
		## Bootstrap
		##==========
		if( self$bootstrap$n_bootstrap > 0 )
		{
			self$bootstrap$coef_ = NULL
			n_sample = base::length(Y)
			for( i in 1:self$bootstrap$n_bootstrap )
			{
				kwargs$resample = base::sample( n_sample , n_sample , TRUE )
				kwargs$Y        = Y[kwargs$resample]
				base::do.call( private$fit_global , kwargs )
				self$bootstrap$coef_ = base::rbind( self$bootstrap$coef_ , self$coef_ )
			}
			self$bootstrap$confidence_interval = base::apply( self$bootstrap$coef_ , 2 , quantile , probs = base::c( self$bootstrap$alpha / 2 , 1 - self$bootstrap$alpha / 2 ) )
		}
		
		## Classic fit
		##============
		kwargs$Y        = Y
		kwargs$resample = NA
		base::do.call( private$fit_global , kwargs )
		
		invisible(self)
	}
	##}}}
	
	),
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	kinds_params = function( kp ) ##{{{
	{
		if( missing(kp) )
			return(private$kinds_params_)
	},
	##}}}
	
	method = function(m) ##{{{
	{
		if( missing(m) )
			return( private$method_ )
		
		private$method_ = base::tolower(m)
	},
	##}}}
	
	coef_ = function(coef) ##{{{
	{
		if( missing(coef) )
			return(self$params$coef_)
		self$params$update_coef(coef)
	}
	##}}}
	
	)
	##}}}
	#################

)
##}}}


