
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

## Link Functions
##===============

## Base class for Link {{{

#' Base class for Link function
#'
#' Base class used to define generic link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{LinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of LinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of LinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of LinkFct at x}
#' }
#' @examples
#'
#' @export
AbstractLink = R6::R6Class( "AbstractLink" ,
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
	},

	eval = function(x)
	{
	},
	
	inverse = function(x)
	{
	},
	
	gradient = function(x)
	{
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## IdLink {{{

#' IdLink
#'
#' Identity link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{IdLinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of IdLinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of IdLinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of IdLinkFct at x}
#' }
#' @examples
#'
#' @export
IdLink = R6::R6Class( "IdLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		super$initialize()
	},
	
	eval = function(x)
	{
		return(x)
	},
	
	inverse = function(x)
	{
		return(x)
	},
	
	gradient = function(x)
	{
		return( base::rep( 1. , length(x) ) )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}

## ExpLink {{{

#' ExpLink
#'
#' Exponential link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new()}}{This method is used to create object of this class with \code{ExpLinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of ExpLinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of ExpLinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of ExpLinkFct at x}
#' }
#' @examples
#'
#' @export
ExpLink = R6::R6Class( "ExpLink" ,
	
	inherit = SDFC::AbstractLink,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		super$initialize()
	},
	
	eval = function(x)
	{
		return( base::exp(x) )
	},
	
	inverse = function(x)
	{
		return( base::log(x) )
	},
	
	gradient = function(x)
	{
		return( base::exp(x) )
	}
	
	
	),
	
	
	######################
	## Private elements ##
	######################
	
	private = list(
	
	###############
	## Arguments ##
	###############
	
	)
)
##}}}



## Parameters
##===========

AbstractParam = R6::R6Class( "AbstractParam" , ##{{{
	
	#################
	## Public list ##
	##{{{
	
	public = list(
	
	## Arguments
	##==========
	kind       = NULL,
	link       = NULL,
	n_samples  = NULL,
	n_features = NULL,
	coef_      = NULL,
	fit_       = NULL,
	
	## Constructor
	##============
	initialize = function( kind , n_samples , ... )
	{
		self$kind = kind
		self$n_samples = n_samples
		self$n_features = 0
		self$coef_      = NULL
		self$fit_       = NULL
		
		kwargs = list(...)
		name_link = base::paste( "l_" , self$kind , sep = "" )
		if( name_link %in% base::names(kwargs) )
			self$link = kwargs[[name_link]]
		else
			self$link = SDFC::IdLink$new()
	},
	
	
	## Methods
	##========
	is_fix = function()
	{},
	
	gradient = function()
	{
		return( self$link$gradient(self$fit_) )
	},
	
	set_coef = function( coef_ )
	{}
	
	),
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	value = function(value_)
	{
		if( missing(value_) )
			return( self$link$eval(self$fit_) )
	}
	
	)
	##}}}
	#################
)
##}}}

CovariateParam = R6::R6Class( "CovariateParam" , ##{{{
	
	inherit = SDFC::AbstractParam,
	
	#################
	## Public list ##
	##{{{
	
	public = list(
	
	## Arguments
	##==========
	design_ = NULL,
	
	## Constructor
	##============
	initialize = function( kind , n_samples , resample , ... )
	{
		## Initialize mother class
		kwargs = list(...)
		kwargs$kind      = kind
		kwargs$n_samples = n_samples
		base::do.call( super$initialize , kwargs )
		
		## Build design matrix
		name_cov = base::paste( "c_" , self$kind , sep = "" )
		X = kwargs[[name_cov]]
		
		if( !is.matrix(X) )
			X = matrix( X , nrow = length(X) , ncol = 1 )
		self$n_features = base::ncol(X) + 1
		
		
		self$design_ = base::cbind( 1 , X )
		if( !base::any(is.na(resample)) )
			self$design_ = self$design_[resample,]
		
		self$coef_   = base::rep( 0 , self$n_features )
		
		if( base::qr(self$design_)$rank < self$n_features )
		{
			self$design_ = matrix( 1 , nrow = self$n_samples , ncol = 1 )
			self$coef_   = 0
		}
	},
	
	
	## Methods
	##========
	is_fix = function()
	{
		return(FALSE)
	},
	
	update = function()
	{
		self$fit_ = matrix( self$design_ %*% self$coef_ , ncol = 1 )
	},
	
	set_intercept = function( coef_ )
	{
		self$coef_[1] = coef_
		self$update()
	},
	
	set_coef = function( coef_ )
	{
		self$coef_ = as.vector(coef_)
		self$update()
	},
	
	design_wo1 = function()
	{
		return(self$design_[,2:self$n_features])
	}
	
	),
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	
	)
	##}}}
	#################
)
##}}}

StationaryParam = R6::R6Class( "StationaryParam" , ##{{{
	
	inherit = SDFC::AbstractParam,
	
	#################
	## Public list ##
	##{{{
	
	public = list(
	
	## Arguments
	##==========
	design_ = NULL,
	
	## Constructor
	##============
	initialize = function( kind , n_samples , resample , ... )
	{
		## Initialize mother class
		kwargs = list(...)
		kwargs$kind      = kind
		kwargs$n_samples = n_samples
		base::do.call( super$initialize , kwargs )
		
		## Build design matrix
		self$n_features = 1
		self$design_ = matrix( 1 , nrow = self$n_samples , ncol = 1 )
		self$coef_   = 0
		self$fit_ = matrix( 0 , nrow = self$n_samples , ncol = 1 )
		
	},
	
	
	## Methods
	##========
	is_fix = function()
	{
		return(FALSE)
	},
	
	update = function()
	{
		self$fit_ = matrix( base::rep( self$coef_ , self$n_samples ) , ncol = 1 )
	},
	
	set_intercept = function( coef_ )
	{
		self$set_coef(coef_)
	},
	
	set_coef = function( coef_ )
	{
		self$coef_ = as.vector(coef_)
		self$update()
	},
	
	design_wo1 = function()
	{
		return(NULL)
	}
	
	),
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	
	)
	##}}}
	#################
)
##}}}

FixParam = R6::R6Class( "FixParam" , ##{{{
	
	inherit = SDFC::AbstractParam,
	
	#################
	## Public list ##
	##{{{
	
	public = list(
	
	## Arguments
	##==========
	
	## Constructor
	##============
	initialize = function( kind , n_samples , resample , ... )
	{
		## Initialize mother class
		kwargs = list(...)
		kwargs$kind      = kind
		kwargs$n_samples = n_samples
		base::do.call( super$initialize , kwargs )
		
		f_par = kwargs[[base::paste( "f_" , self$kind , sep = "" )]]
		self$fit_ = as.vector( self$link$inverse(f_par) )
		if( length(self$fit_) == 1 )
			self$fit_ = base::rep( self$fit_ , self$n_samples )
		
		if( !is.na(resample) )
			self$fit_ = self$fit_[resample]
		
		self$fit_ = matrix( self$fit_ , ncol = 1 )
		
	},
	
	
	## Methods
	##========
	is_fix = function()
	{
		return(TRUE)
	}
	
	),
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	
	)
	##}}}
	#################
)
##}}}

LawParams = R6::R6Class( "LawParams" , ##{{{
	
	#################
	## Public list ##
	##{{{
	
	public = list(
	
	## Arguments
	##==========
	
	kinds    = NULL,
	dparams_ = NULL,
	coef_    = NULL,
	
	## Constructor
	##============
	initialize = function( kinds , ... )
	{
		self$kinds    = kinds
		self$dparams_ = list()
	},
	
	
	## Methods
	##========
	
	add_params = function( n_samples , resample , ... )##{{{
	{
		kwargs = list(...)
		for( kind in self$kinds )
		{
			kwargs$kind = kind
			k_param = base::do.call( self$filter , kwargs )
			config  = base::do.call( self$infer_configuration , k_param )
			k_param$kind      = kind
			k_param$n_samples = n_samples
			k_param$resample  = resample
			if( self$is_covariate(config)  ) self$dparams_[[kind]] = base::do.call( CovariateParam$new  , k_param ) 
			if( self$is_stationary(config) ) self$dparams_[[kind]] = base::do.call( StationaryParam$new , k_param ) 
			if( self$is_fix(config)        ) self$dparams_[[kind]] = base::do.call( FixParam$new        , k_param ) 
		}
	},
	##}}}
	
	merge_covariate = function()##{{{
	{
		l_c = list()
		for( k in base::names(self$dparams_) )
		{
			if( "CovariateParam" %in% class(self$dparams_[[k]]) )
			{
				l_c[[k]] = matrix( self$dparams_[[k]]$design_wo1() , nrow = self$n_samples )
			}
		}
		
		if( length(l_c) == 0 )
			return(NULL)
		
		C = matrix( 1 , nrow = self$n_samples , ncol = 1 )
		
		for( c in l_c )
		{
			Cnext = base::cbind( C , c )
			if( base::qr(Cnext) == base::ncol(Cnext) )
				C = Cnext
		}
		
		if( base::ncol(C) == 1 )
			return(NULL)
		else
			return( C[,2:base::ncol(C)] )
	},
	##}}}
	
	merge_coef = function()##{{{
	{
		self$coef_ = base::c()
		for( k in base::names(self$dparams_) )
		{
			if( !self$dparams_[[k]]$is_fix() )
				self$coef_ = base::c( self$coef_ , self$dparams_[[k]]$coef_ )
		}
	},
	##}}}
	
	split_coef = function( coef )##{{{
	{
		tcoef = list()
		a = 1
		b = 1
		for( k in base::names(self$dparams_) )
		{
			if( self$dparams_[[k]]$is_fix() )
			{
				tcoef[[k]] = NULL
			}
			else
			{
				b = a + self$dparams_[[k]]$n_features - 1
				tcoef[[k]] = coef[a:b]
				a = b + 1
			}
		}
		return(tcoef)
	},
	##}}}
	
	update_coef = function( coef , kind = NULL )##{{{
	{
		if( is.null(kind) )
		{
			lcoef = self$split_coef(coef)
			for( k in base::names(self$dparams_) )
				self$dparams_[[k]]$set_coef(lcoef[[k]])
		}
		else
		{
			self$dparams_[[kind]]$set_coef(coef)
		}
		self$merge_coef()
	},
	##}}}
	
	set_intercept = function( coef , kind )##{{{
	{
		self$dparams[[kind]]$set_intercept(coef)
		self$merge_coef()
	},
	##}}}
	
	
	
	## Methods to infer configuration
	##==============================={{{
	
	infer_configuration = function(...)
	{
		kwargs = list(...)
		has_c = FALSE
		has_s = FALSE
		has_f = FALSE
		has_l = FALSE
		for( k in base::names(kwargs) )
		{
			if( base::substr( k , 1 , 2 ) == "c_" ) has_c = TRUE
			if( base::substr( k , 1 , 2 ) == "f_" ) has_f = TRUE
			if( base::substr( k , 1 , 2 ) == "l_" ) has_l = TRUE
		}
		
		has_s = length(kwargs) == 0 || ( length(kwargs) == 1 && has_l )
		
		return( list( has_c = has_c , has_s = has_s , has_f = has_f , has_l = has_l ) )
	},
	
	is_covariate = function( c )
	{
		return( c$has_c && !( c$has_s || c$has_f ) )
	},
	
	is_stationary = function( c )
	{
		return( c$has_s && !( c$has_c || c$has_f ) )
	},
	
	is_fix = function( c )
	{
		return( c$has_f && !( c$has_c || c$has_s ) )
	},
	
	filter = function( kind , ... )
	{
		kwargs = list(...)
		out = list()
		for( k in base::names(kwargs) )
		{
			if( base::grepl( kind , k , TRUE ) )
				out[[k]] = kwargs[[k]]
		}
		return(out)
	}
	##}}}
	
	),
	
	##}}}
	#################
	
	#################
	## Active list ##
	##{{{
	
	active = list(
	
	
	)
	##}}}
	#################
)
##}}}


## Abstract Law
##=============

AbstractLaw2 = R6::R6Class( "AbstractLaw2" ,##{{{

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

Normal = R6::R6Class( "Normal" ,##{{{
	
	inherit = AbstractLaw2,
	
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
	initialize = function( method , n_bootstrap = 0 , alpha = 0.05 )
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




