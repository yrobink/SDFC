
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
	
	new_screen = function()
	{
		if( self$os == "osx" )
		{
			grDevices::quartz()
		}
		if( self$os == "linux" )
		{
			grDevices::X11()
		}
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
	c_data = Dataset$covariates(size)
	
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

test_gev = function() ##{{{
{
	size  = 2500
	t = base::seq( 0 , 1 , length = size )
	X0 = t^2
	X2 = base::seq( -1 , 1 , length = size )
	loc   = 0.5 + 1.5 * X0
	scale = 0.1 + 0.1 * X0
	shape = 0.3 * X2
	
	Y = SDFC::rgev( n = size , loc = loc , scale = scale , shape = shape )
	
	gev = GEVLaw$new( method = "MLE" , n_bootstrap = 10 )
	gev$fit( Y , loc_cov = X0  , scale_cov = X0    , shape_cov = X2 )
	#gev$fit( Y , floc    = loc , scale_cov = X0    , shape_cov = X2 )
	#gev$fit( Y , loc_cov = X0  , fscale    = scale , shape_cov = X2 )
	#gev$fit( Y , loc_cov = X0  , scale_cov = X0    , fshape    = shape )
	
	plt$new_screen()
	graphics::par( mfrow = base::c( 2 , 2 ) )
	
	plot( t , Y , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( loc , gev$loc , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( scale , gev$scale , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( shape , gev$shape , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
}
##}}}

test_gpd = function()##{{{
{
	size  = 2500
	t = base::seq( 0 , 1 , length = size )
	X0 = t^2
	X2 = base::seq( -1 , 1 , length = size )
	loc   = 0.5 + 1.5 * X0
	scale = 0.1 + 0.1 * X0
	shape = 0.3 * X2
	
	Y = rgpd( n = size , loc = loc , scale = scale , shape = shape )
	
	gpd = GPDLaw$new( method = "MLE" , n_bootstrap = 10 )
	gpd$fit( Y , loc = loc , scale_cov = X0 , shape_cov = X2 )
	
	
	plt$new_screen()
	graphics::par( mfrow = base::c( 2 , 2 ) )
	
	plot( t , Y , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( loc , gpd$loc , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( scale , gpd$scale , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( shape , gpd$shape , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
}
##}}}

test_qr = function( show = FALSE ) ##{{{
{
	size = 2000
	c_data = Dataset$covariates(size)
	
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
	c_data = Dataset$covariates(size)
	
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

run_all_tests = function()##{{{
{
	test_normal()
	test_gev()
	test_gpd()
	test_qr()
}
##}}}


##########
## main ##
##########

#run_all_tests()

#test_normal(TRUE)

size = 2000
c_data = Dataset$covariates(size)

t       = c_data$t
X_loc   = c_data$X_loc
X_scale = c_data$X_scale
X_shape = c_data$X_shape

loc   = 1.  + 0.8  * X_loc
scale = 0.2 + 0.08 * X_scale
shape = 0.  + 0.3  * X_shape


Y = SDFC::rgev( size , loc , scale , shape )

gev = SDFC::GEV$new( method = "lmoments-experimental" )
gev$fit( Y , c_loc = X_loc , c_shape = X_shape , c_scale = X_scale )
print(gev$coef_)

plt$wait()
print("Done")

