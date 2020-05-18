
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



#' Tabular class
#'
#' Class used to print tabular in terminal
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
#'   \item{\code{new()}}{This method is used to create object of this class with \code{Tabular}}
#'   \item{\code{add_row(row)}}{add a row to tabular}
#'   \item{\code{draw()}}{Return a character containing tabular}
#' }
#' @examples
#'
#' A = matrix( rnorm(20) , nrow = 5 , ncol = 4 )
#' tab = SDFC::Tabular$new()
#' tab$header = base::c( "A","B","C","D")
#' for( i in 1:5 )
#'     tab$add_row( A[i,] )
#' print(tab$draw())
#' 
#' @export
Tabular = R6::R6Class( "Tabular" ,
	
	## Private list
	##=========={{{
	
	private = list(
	
	.header = NULL,
	.rows   = NULL,
	.ncol   = NULL,
	
	build_size_cells = function() ##{{{
	{
			size_cells = base::rep( 0 , self$ncol )
		
		if( !is.null(self$header) )
		{
			i = 1
			for( h in self$header )
			{
				size_cells[i] = max( size_cells[i] , base::nchar(as.character(h)) )
				i = i + 1
			}
		}
		
		for( row in private$.rows )
		{
			i = 1
			for( r in row )
			{
				size_cells[i] = max( size_cells[i] , base::nchar(as.character(r)) )
				i = i + 1
			}
		}
		return(size_cells)
	},
	##}}}
	
	build_separator = function( size_cells , kind )##{{{
	{
		separator = base::c("+")
		for( k in size_cells )
			separator = base::c( separator , base::rep( kind , k + 2 ) , "+")
		separator = base::paste( separator , collapse = "" )
		return(separator)
	},
	##}}}
	
	build_line = function( row , size_cells , align = "right" )##{{{
	{
		line = base::c("|")
		i = 1
		for( r in row )
		{
			if( align == "right" )
				line = base::c( line , base::rep( " " , size_cells[i] - nchar(r) + 1 ) , as.character(r) , " |" )
			if( align == "center" )
			{
				nblanks = size_cells[i] - nchar(r)
				nblanks_left  = as.integer(base::round(nblanks/2))
				nblanks_right = nblanks - nblanks_left
				line = base::c( line , base::rep( " " , nblanks_left + 1 ) , as.character(r) , base::rep( " " , nblanks_right + 1 ) ,  "|" )
			}
			if( align == "left" )
				line = base::c( line , " " , as.character(r) , base::rep( " " , size_cells[i] - nchar(r) ) , " |" )
			i = i + 1
		}
		line = base::paste( line , collapse = "" )
		return(line)
	}
	##}}}
	
	),
	##}}}
	
	## Active list
	##========={{{
	
	active  = list(
	
	header = function(value)##{{{
	{
		if( base::missing(value) )
		{
			return(private$.header)
		}
		else
		{
			if( is.null(self$ncol) )
			{
				private$.ncol   = length(value)
				private$.header = value
			}
			else if( length(value) == self$ncol )
			{
				private$.header = row
			}
		}
	},
	##}}}
	
	ncol = function(value)##{{{
	{
		if( base::missing(value) )
		{
			return(private$.ncol)
		}
	},
	##}}}
	
	nrow = function(value)##{{{
	{
		if( base::missing(value) )
			return( length(private$.rows) )
	}
	##}}}
	
	),
	##}}}
	
	## Public list
	##========={{{
	
	public  = list(
	
	initialize = function()##{{{
	{
		private$.rows  = list()
	},
	##}}}
	
	add_row = function( row )##{{{
	{
		if( is.null(self$ncol) )
		{
			private$.ncol      = length(row)
			private$.rows[[1]] = row
		}
		else if( length(row) == self$ncol )
		{
			private$.rows[[self$nrow + 1]] = row
		}
		invisible(NULL)
	},
	##}}}
	
	draw = function() ##{{{
	{
		tabular = ""
		size_cells = private$build_size_cells()
		sep_header = private$build_separator(size_cells,"=")
		sep_lines  = private$build_separator(size_cells,"-")
		
		if( !is.null(self$header) )
		{
			line = private$build_line( self$header , size_cells , align = "center" )
			tabular = base::paste( sep_header , line , sep_header , sep = "\n" )
		}
		
		if( self$nrow > 0 )
		{
			for( row in private$.rows )
			{
				line = private$build_line( row , size_cells )
				tabular = base::paste( tabular , line , sep_lines , sep = "\n" )
			}
		}
		
		tabular = base::paste( tabular , "\n" , sep = "" )
		
		invisible(tabular)
	}
	##}}}
	
	)
	##}}}
)
