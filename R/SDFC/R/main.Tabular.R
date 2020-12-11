
## Copyright(c) 2020 Yoann Robin
## 
## This file is part of SDFC.
## 
## SDFC is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## SDFC is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with SDFC.  If not, see <https://www.gnu.org/licenses/>.


#' Tabular class
#'
#' @description
#' Class printing generating a tabular
#'
#' @details
#' Class used to generate a tabular
#'
#' @importFrom R6 R6Class
#' @importFrom methods new
#'
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
	
	#' @field header [vector of string] Header of the tabular
	.header = NULL,
	.rows   = NULL,
	#' @field ncol [integer] Number of columns
	.ncol   = NULL,
	#' @field nrow [integer] Number of rows
	.nrow   = NULL,
	
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
	
	## initialize ##{{{
	#' @description
    #' Create a new Tabular object.
	#' @return A new `Tabular` object.
	initialize = function()
	{
		private$.rows  = list()
	},
	##}}}
	
	## add_row ##{{{
	#' @description
    #' Add row to tabular
    #' @param row [vector] Vector of element to print
    #' @return NULL
	add_row = function( row )
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
	
	## draw ##{{{ 
	#' @description
    #' Generate the tabular
    #' @return [string] Return the tabular in string form
	draw = function()
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
