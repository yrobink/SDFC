
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

try(roxygen2::roxygenize("../SDFC"))
devtools::load_all("../SDFC")
roxygen2::roxygenize("../SDFC")
devtools::load_all("../SDFC")
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

##########
## main ##
##########

## Build data
n_samples = 2500
data = SDFC::dataset(n_samples)
coef_ = base::c(0.5,2.,0.5,-0.2)
loc   = coef_[1] + coef_[2] * data$X_loc
scale = base::exp( coef_[3] + coef_[4] * data$X_scale )
Y     = stats::rnorm( n = n_samples , mean = loc , sd = scale )

## And now the script
kwargs = list( c_loc = data$X_loc , c_scale = data$X_scale , l_scale = SDFC::ULExponential$new() )

lhs = LHS$new( base::c("loc","scale") , 2500 )
rhs = RHS$new(lhs)

base::do.call( rhs$build , kwargs )
rhs$coef_ = coef_

law = SDFC::Normal$new()


plt$new_screen( ncol = 2 )
plot( loc , rhs$lhs$values$loc , type = "p" )
plot( scale , rhs$lhs$values$scale , type = "p" )


## End
plt$wait()
print("Done")

