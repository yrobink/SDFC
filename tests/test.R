
#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

base::rm( list = base::ls() )


###############
## Libraries ##
###############

library(R6)
library(devtools)

try(roxygen2::roxygenize("../R/SDFC"))
devtools::load_all("../R/SDFC")
roxygen2::roxygenize("../R/SDFC")
devtools::load_all("../R/SDFC")
#library(SDFC)


###########################
## Useful plot functions ##
###########################

PlotTools = R6::R6Class( "PlotTools" , ##{{{
	
	
	public = list(
	
	###############
	## Arguments ##
	###############
	
	os = NULL,
	
	
	#################
	## Constructor ##
	#################
	
	initialize = function()
	{
		self$os = self$get_os()
	},
	
	
	#############
	## Methods ##
	#############
	
	get_os = function()
	{
		sysinf = base::Sys.info()
		if( !is.null(sysinf) )
		{
			os = sysinf['sysname']
			if( os == 'Darwin' ) os = "osx"
		}
		else
		{
			## mystery machine
			os = .Platform$OS.type
			if( base::grepl( "^darwin"   , R.version$os ) ) os = "osx"
			if( base::grepl( "linux-gnu" , R.version$os ) ) os = "linux"
		}
		invisible(tolower(os))
	},
	
	new_screen = function( nrow = 1 , ncol = 1 )
	{
		if( self$os == "osx" )
		{
			grDevices::quartz()
		}
		if( self$os == "linux" )
		{
			grDevices::X11()
		}
		
		graphics::par( mfrow = base::c( nrow , ncol ) )
	},
	
	wait = function()
	{
		while( base::names(grDevices::dev.cur()) !='null device' ) base::Sys.sleep(1)
	}
	
	)
)
##}}}

plt = PlotTools$new()


###############
## Functions ##
###############

test_normal = function( show = FALSE ) ##{{{
{
	size = 2000
	c_data = dataset.covariates(size)
	
	t       = c_data$t
	X_loc   = c_data$X_loc
	X_scale = c_data$X_scale
	loc   = 0.5 + 2 * X_loc
	scale =   1 + 2 * X_scale
	Y = stats::rnorm( size , mean = loc , sd = scale )
	
	
	law = SDFC::Normal$new( "mle" )
	law$fit( Y , c_loc = X_loc , c_scale = X_scale )
	
	if(show)
	{
		plt$new_screen()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Subplot 1
		graphics::plot(  t , Y       , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) )
		graphics::lines( t , law$loc , col = "red" )
		graphics::lines( t , law$loc + law$scale , col = "red" )
		graphics::lines( t , law$loc - law$scale , col = "red" )
		
		## Subplot 2
		graphics::plot( loc , law$loc , type = "p" , col = "blue" )
		
		## Subplot 3
		graphics::plot( scale , law$scale , type = "p" , col = "blue" )
	}
}
##}}}

test_exponential = function( show = FALSE )##{{{
{
	size = 10000
	c_data = dataset.covariates(size)
	
	t       = c_data$t
	X_scale = c_data$X_scale
	
	scale = 0.5 + 1.5 * X_scale
	
	
	Y = rexp( size , rate = 1. / scale )
	
	law = SDFC::Exponential$new( method = "mle" )
	law$fit( Y , c_scale = X_scale )
	
	if(show)
	{
		plt$new_screen( 1 , 2 )
		
		## Subplot 1
		graphics::plot(  t , Y       , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) )
		graphics::lines( t , law$scale , col = "red" )
		
		## Subplot 2
		graphics::plot( scale , law$scale , type = "p" , col = "blue" )
	}
}
##}}}

test_gamma = function( show = FALSE )##{{{
{
	size = 2000
	c_data = dataset.covariates(size)
	
	t       = c_data$t
	X_scale = c_data$X_scale
	X_shape = c_data$X_shape
	
	scale = 0.2 + 0.08 * X_scale
	shape = 0.4 + 0.3  * X_shape
	
	
	Y = rgamma( size , shape = shape , scale = scale )
	
	# Regression with MLE
	law = SDFC::Gamma$new( "mle" )
	law$fit( Y , c_scale = X_scale , c_shape = X_shape )
	
	if( show )
	{
		plt$new_screen(2,2)
		
		## Subplot 1
		graphics::plot(  t , Y       , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) )
		graphics::lines( t ,  law$scale , col = "red" )
		graphics::lines( t , -law$scale , col = "red" )
		
		## Subplot 2
		graphics::plot( scale , law$scale , type = "p" , col = "blue" )
		
		## Subplot 3
		graphics::plot( shape , law$shape , type = "p" , col = "blue" )
		
	}
}
##}}}

test_gev = function( show = FALSE ) ##{{{
{
	size = 2000
	c_data = dataset.covariates(size)
	
	t       = c_data$t
	X_loc   = c_data$X_loc
	X_scale = c_data$X_scale
	X_shape = c_data$X_shape
	
	loc   = 1.  + 0.8  * X_loc
	scale = 0.2 + 0.08 * X_scale
	shape = 0.  + 0.3  * X_shape
	
	
	Y = SDFC::rgev( size , loc , scale , shape )
	
	law = SDFC::GEV$new( method = "mle" )
	law$fit( Y , c_loc = X_loc , c_shape = X_shape , c_scale = X_scale )
	
	if(show)
	{
		plt$new_screen()
		graphics::par( mfrow = base::c( 2 , 2 ) )
		
		## Subplot 1
		graphics::plot(  t , Y       , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) )
		graphics::lines( t , law$loc , col = "red" )
		graphics::lines( t , law$loc + law$scale , col = "red" )
		graphics::lines( t , law$loc - law$scale , col = "red" )
		graphics::lines( t , law$loc - law$scale / law$shape , col = "black" )
		
		## Subplot 2
		graphics::plot( loc , law$loc , type = "p" , col = "blue" )
		
		## Subplot 3
		graphics::plot( scale , law$scale , type = "p" , col = "blue" )
		
		## Subplot 4
		graphics::plot( shape , law$shape , type = "p" , col = "blue" )
	}
}
##}}}

test_gpd = function( show = FALSE )##{{{
{
	size = 2000
	c_data = dataset.covariates(size)
	
	t       = c_data$t
	X_loc   = c_data$X_loc
	X_scale = c_data$X_scale
	X_shape = c_data$X_shape
	
	loc   = 1.  + 0.8  * X_loc
	scale = 0.2 + 0.08 * X_scale
	shape = 0   - 0.2  * X_shape
	
	
	Y = SDFC::rgpd( size , loc , scale , shape )
	
	law = GPD$new( method = "mle" )
	law$fit( Y , f_loc = loc , c_scale = X_scale , c_shape = X_shape )
	
	if( show )
	{
		plt$new_screen(2,2)
		
		## Subplot 1
		graphics::plot(  t , Y       , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) )
		graphics::lines( t , law$loc , col = "red" )
		graphics::lines( t , law$loc + law$scale , col = "red" )
		graphics::lines( t , law$loc - law$scale , col = "red" )
		graphics::lines( t , law$loc - law$scale / law$shape , col = "black" )
		
		## Subplot 2
		graphics::plot( loc , law$loc , type = "p" , col = "blue" )
		
		## Subplot 3
		graphics::plot( scale , law$scale , type = "p" , col = "blue" )
		
		## Subplot 4
		graphics::plot( shape , law$shape , type = "p" , col = "blue" )
		
	}
}
##}}}

test_qr = function( show = FALSE ) ##{{{
{
	size = 2000
	c_data = dataset.covariates(size)
	
	loc   = 0.5 + 2 * c_data$X_loc
	scale = 1 + 2 * c_data$X_scale
	Y = stats::rnorm( size , mean = loc , sd = scale )
	
	q  = base::seq( 0.01 , 0.99 , length = 100 )
	Yq = np_quantile( Y , q , base::cbind( c_data$X_loc , c_data$X_scale ) )
	
	if(show)
	{
		plt$new_screen()
		graphics::par( mfrow = base::c( 1 , 1 ) )
		plot( c_data$t , Y , col = "blue" , type = "p" )
		for( i in 1:length(q) )
		{
			lines( c_data$t , Yq[,i] , col = "red" )
		}
	}
}
##}}}

test_lmoments = function() ##{{{
{
	size = 2000
	c_data = dataset.covariates(size)
	
	t       = c_data$t
	X_loc   = c_data$X_loc
	X_scale = c_data$X_scale
	loc   = 0.5 + 2 * X_loc
	scale =   1 + 2 * X_scale
	Y = stats::rnorm( size , mean = loc , sd = scale )
	
	c_Y = base::cbind( X_loc , X_scale )
	
	lmom = np_lmoments( Y , c_Y = c_Y )
}
##}}}

run_all_tests = function( show = TRUE )##{{{
{
	test_normal(show)
	test_exponential(show)
	test_gamma(show)
	test_gev(show)
	test_gpd(show)
	test_qr(show)
	test_lmoments()
}
##}}}


##########
## main ##
##########

#run_all_tests()
#plt$wait()




print("Done")

