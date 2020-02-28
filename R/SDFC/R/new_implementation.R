
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

## ChainLinkFct {{{

#' ChainLinkFct
#'
#' Chain link function to chain two link functions
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param x [vector]
#' @param linFct1 [LinkFct] Second link function to apply
#' @param linFct0 [LinkFct] First link function to apply
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(linkFct1,linkFct0)}}{This method is used to create object of this class with \code{IdLinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of IdLinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of IdLinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of IdLinkFct at x}
#' }
#' @examples
#'
#' @export
ChainLinkFct = R6::R6Class( "ChainLinkFct" ,
	
	inherit = SDFC::LinkFct,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	linkFct0 = NULL,
	linkFct1 = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function( linkFct1 , linkFct0 )
	{
		super$initialize()
		self$linkFct0 = linkFct0
		self$linkFct1 = linkFct1
	},
	
	eval = function(x)
	{
		return( self$linkFct1$eval( self$linkFct0$eval(x) ) )
	},
	
	inverse = function(x)
	{
		return( self$linkFct0$inverse( self$linkFct1$inverse(x) ) )
	},
	
	gradient = function(x)
	{
		return( self$linkFct0$gradient(x) * self$linkFct1$gradient( self$linkFct0$eval(x) ) )
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

## InverseLinkFct {{{

#' InverseLinkFct
#'
#' Inverse link function
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
#'   \item{\code{new()}}{This method is used to create object of this class with \code{InverseLinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of InverseLinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of InverseLinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of InverseLinkFct at x}
#' }
#' @examples
#'
#' @export
InverseLinkFct = R6::R6Class( "InverseLinkFct" ,
	
	inherit = SDFC::LinkFct,
	
	
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
		return( 1. / x )
	},
	
	inverse = function(x)
	{
		return( 1. / x )
	},
	
	gradient = function(x)
	{
		return( - 1. / x^2 )
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

## LogitLinkFct {{{

#' LogitLinkFct
#'
#' Logit link function
#'
#' @docType class
#' @importFrom R6 R6Class
#'
#' @param a [float] Lower bound of logit
#' @param b [float] Upper bound of logit
#' @param s [float] Speed of logit between a and b
#' @param x [vector]
#'
#' @return Object of \code{\link{R6Class}}
#' @format \code{\link{R6Class}} object.
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(a,b,s)}}{This method is used to create object of this class with \code{LogitLinkFct}}
#'   \item{\code{eval(x)}}{Evaluation of LogitLinkFct at x}
#'   \item{\code{inverse(x)}}{Inverse of LogitLinkFct at x}
#'   \item{\code{gradient(x)}}{Gradient of LogitLinkFct at x}
#' }
#' @examples
#'
#' @export
LogitLinkFct = R6::R6Class( "LogitLinkFct" ,
	
	inherit = SDFC::LinkFct,
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	a = 0,
	b = 1,
	s = 1,
	
	#################
	## Constructor ##
	#################
	
	initialize = function( a = 0 , b = 1 , s = 1 )
	{
		super$initialize()
		self$a = a
		self$b = b
		self$s = s
	},
	
	eval = function(x)
	{
		return( (self$b - self$a) / ( 1. + base::exp(- self$s * x) ) + self$a )
	},
	
	inverse = function(x)
	{
		idx_lo = x < self$a
		idx_up = x > self$b
		x[idx_lo] = self$a + 1e-3
		x[idx_up] = self$b - 1e-3
		return( - base::log( (self$b - self$a) / (x - self$a) - 1 ) / self$s )
	},
	
	gradient = function(x)
	{
		e = base::exp( - self$s * x )
		return( self$s * (self$b - self$a) * e / ( 1 + e )^2 )
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


