
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


roxygenize("SDFC")
load_all("SDFC")
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



##########
## main ##
##########

test_normal()


