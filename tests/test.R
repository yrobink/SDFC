
#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

base::rm( list = base::ls() )


###############
## Libraries ##
###############

library(roxygen2)
library(devtools)


roxygenize("../R/SDFC")
load_all("../R/SDFC")
#check("SDFC")


###############
## Functions ##
###############

test_normal = function() ##{{{
{
	size = 2500
	t    = base::seq( 0 , 1 , length = size )
	X0    = t^2
	X1    = base::cos( 2 * base::pi * t )
	loc   = 1. + 2 * X0
	scale = 0.6 + 0.5 * X1
	Y    = stats::rnorm( n = size , mean = loc , sd = scale )
	
	
	law = SDFC::NormalLaw$new( method = "MLE" ,  n_bootstrap = 10 )
	law$fit( Y , loc_cov = X0 , scale_cov = X1 )
	
	graphics::par( mfrow = base::c( 2 , 2 ) )
	
	## Subplot 1
	graphics::plot( t , Y , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) )
	graphics::lines( t , law$loc , col = "red" )
	graphics::lines( t , law$loc + law$scale , col = "red" )
	graphics::lines( t , law$loc - law$scale , col = "red" )
	
	## Subplot 2
	graphics::plot( loc , law$loc , type = "p" , col = "blue" )
	
	## Subplot 3
	graphics::plot( scale , law$scale , type = "p" , col = "blue" )
	
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
	
	
	graphics::par( mfrow = base::c( 2 , 2 ) )
	
	plot( t , Y , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( loc , gpd$loc , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( scale , gpd$scale , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
	plot( shape , gpd$shape , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	
}
##}}}

test_qr = function() ##{{{
{
	## Generate data
	size = 2500
	t    = base::seq( 0 , 1 , length = size )
	X0   = t^2
	X1   = base::cos( 2 * base::pi * t )
	
	loc   = 1. + 2 * X0
	scale = 0.6 + 0.5 * X1
	Y     = stats::rnorm( n = size , mean = loc , sd = scale )
	
	## QR
	ltau = base::seq( 0.01 , 0.99 , length = 100 )
	qr   = SDFC::np_quantile( Y , ltau = ltau , X = base::cbind( X0 , X1 ) )
	
	
	## Plot
	graphics::par( mfrow = base::c(1,1) )
	graphics::plot( t , Y , col = grDevices::rgb( 0 , 0 , 1 , 0.5 ) , type = "p" )
	for( i in 1:100 )
	{
		graphics::lines( t , qr[,i] , col = grDevices::rgb( 0.5 , 0.5 , 0.5 , 0.5 ) )
	}
}
##}}}

##########
## main ##
##########

test_normal()
test_gev()
test_gpd()
test_qr()
