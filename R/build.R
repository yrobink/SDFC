
library(devtools)
library(roxygen2)


## Read command line arguments
##============================

args = commandArgs( trailingOnly = TRUE )
verbose = FALSE
install = FALSE
build   = TRUE
check   = FALSE
if( length(args) > 0 )
{
	for( a in args )
	{
		if( a == "-v" || a == "--verbose" )
		{
			verbose = TRUE
		}
		if( a == "-i" || a == "--install" )
		{
			install = TRUE
		}
		if( a == "-nb" || a == "--not-build" )
		{
			build = FALSE
		}
		if( a == "-c" || a == "--check" )
		{
			check = TRUE
		}
	}
}


## Building
##=========
cdsk = ""
if( build )
{
	if( verbose ) cat( "Generation of Rd files with roxygen" )
	roxygen2::roxygenize("SDFC")
	if( verbose ) cat( "Load of SDFC to generate Rd files with Rcpp" )
	devtools::load_all("SDFC")
	if( verbose ) cat( "Generation of Rd files for cpp with roxygen" )
	roxygen2::roxygenize("SDFC")
	if( verbose ) cat( "Final build" )
	cdsk = devtools::build("SDFC")
}


## Check
##======
if( check )
{
	if( verbose ) cat( "Check SDFC" )
	devtools::check( "SDFC" )
}


## Installation
##=============

if( install )
{
	if( verbose ) cat( "Installation" )
	if( cdsk == "" )
	{
		files = base::list.files()
		cdsk = ""
		
		for( f in files )
		{
			f_split = base::unlist(base::strsplit( f , "[.]" ))
			if( length(f_split) > 2 )
			{
				f_split = base::rev(f_split)
				if( f_split[1] == "gz" && f_split[2] == "tar" )
				{
					cdsk = f
				}
			}
		}
	}
	if( cdsk == "" )
	{
		cat( "Error, SDFC not build, so can not be install" )
	}
	else
	{
		install.packages(cdsk)
	}
}

