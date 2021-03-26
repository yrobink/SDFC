
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

base::rm( list = base::ls() )


###############
## Libraries ##
###############

library(R6)
library(devtools)
library(ROOPSD)

try(roxygen2::roxygenize("../SDFC"))
devtools::load_all("../SDFC")
roxygen2::roxygenize("../SDFC")
devtools::load_all("../SDFC")


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


#############
## Classes ##
#############

SDFCLawTest = R6::R6Class( "SDFCLawTest" , ##{{{
	
	public = list(
	
	## Attributs ##{{{
	name      = NULL,
	n_samples = NULL,
	sd_law    = NULL,
	shape_p   = NULL,
	t         = NULL,
	X_loc     = NULL,
	X_scale   = NULL,
	X_shape   = NULL,
	has_loc   = NULL,
	has_scale = NULL,
	has_shape = NULL,
	n_params  = NULL,
	kwargs    = NULL,
	coef_     = NULL,
	loc       = NULL,
	scale     = NULL,
	shape     = NULL,
	Y         = NULL,
	law       = NULL,
	
	##}}}
	
	initialize = function(...)##{{{
	{
		kwargs = list(...)
		self$name      = kwargs[["name"]]
		self$n_samples = kwargs[["n_samples"]]
		self$sd_law    = kwargs[["sd_law"]]
		self$shape_p   = kwargs[["shape_p"]]
		data           = SDFC::dataset(self$n_samples)
		self$t         = data$t
		self$X_loc     = data$X_loc
		self$X_scale   = data$X_scale
		self$X_shape   = data$X_shape
		self$has_loc   = "loc"   %in% kwargs[["params"]]
		self$has_scale = "scale" %in% kwargs[["params"]]
		self$has_shape = "shape" %in% kwargs[["params"]]
		self$n_params  = self$has_loc + self$has_scale + self$has_shape
		self$kwargs  = list()
		self$coef_   = base::c()
	},
	##}}}
	
	
	build_loc = function( code )##{{{
	{
		coef_    = base::c(0.5,1.)
		if( code == 1)
		{
			coef_      = coef_[1]
			self$loc   = base::rep( coef_[1] , self$n_samples )
			self$coef_ = base::c(self$coef_,coef_)
		}
		else
		{
			self$loc = coef_[1] + coef_[2] * self$X_loc
			if( code == 0 )
			{
				self$kwargs[["c_loc"]] = self$X_loc
				self$coef_ = base::c(self$coef_,coef_)
			}
			else
			{
				self$kwargs[["f_loc"]] = self$loc
			}
		}
	},
	##}}}
	
	build_scale = function( code )##{{{
	{
		coef_    = base::c(0.3,-0.9)
		if( code == 1)
		{
			coef_      = coef_[1]
			self$scale = base::rep( base::exp(coef_[1]) , self$n_samples )
			self$coef_ = base::c(self$coef_,coef_)
			self$kwargs[["l_scale"]] = ULExponential$new()
		}
		else
		{
			self$scale = base::exp( coef_[1] + coef_[2] * self$X_scale )
			if( code == 0 )
			{
				self$kwargs[["c_scale"]] = self$X_scale
				self$kwargs[["l_scale"]] = ULExponential$new()
				self$coef_ = base::c(self$coef_,coef_)
			}
			else
			{
				self$kwargs[["f_scale"]] = self$scale
			}
		}
	},
	##}}}
	
	build_shape = function( code )##{{{
	{
		if( self.shape_p )
		{
			coef_ = base::c(1,-0.2)
		}
		else
		{
			coef_    = base::c(0.,0.2)
		}
		if( code == 1 )
		{
			coef_      = -coef_[1]
			self$shape = base::rep( coef_[1] , self$n_samples )
			self$coef_ = base::c(self$coef_,coef_)
		}
		else
		{
			self$shape = coef_[1] + coef_[2] * self$X_shape
			if( code == 0 )
			{
				self$kwargs[["c_shape"]] = self$X_shape
				self$coef_ = base::c(self$coef_,coef_)
			}
			else
			{
				self$kwargs[["f_shape"]] = self$shape
			}
		}
	},
	##}}}
	
	
	testXXX = function( code , method = "MLE" )##{{{
	{
		self$kwargs = list()
		self$coef_  = base::c()
		i = 1
		if( self$has_loc )
		{
			self$build_loc( code[i] )
			i = i + 1
		}
		if( self$has_scale )
		{
			self$build_scale( code[i] )
			i = i + 1
		}
		if( self$has_shape )
		{
			self$build_shape( code[i] )
		}
		
		self$Y = self$rvs()
		
		if( length(self$coef_) > 0 )
		{
			self$kwargs[["n_mcmc_drawn"]] = 100
			if( length(self$coef_) == 1 )
			{
				self$kwargs[["prior"]] = ROOPSD::Normal$new( self$coef_ , 0.1 )
			}
			else
			{
				m = self$coef_
				S = base::diag(0.1 + numeric(length(self$coef_)))
				self$kwargs[["prior"]] = MultivariateNormal$new(m,S)
			}
		}
		
		self$law = self$sd_law$new( method = method )
		self$kwargs[["Y"]] = self$Y
		base::do.call( self$law$fit , self$kwargs )
	},
	##}}}
	
	
	run_all = function( method = "MLE" ) ##{{{
	{
		tab = Tabular$new()
		tab$header = base::c( base::paste0(self$name," law test (",method,")"),"Status","Max diff","True value","Estimated value" )
		
		idxkwargs = list()
		for( i in 1:self$n_params )
		{
			idxkwargs = base::c( idxkwargs , list(0:2) )
		}
		indexes = do.call( expand.grid , idxkwargs )
		colnames(indexes) = NULL
		rownames(indexes) = NULL
		indexes = as.matrix(indexes)
		
		for( i in 1:nrow(indexes) )
		{
			idx = indexes[i,]
			test = try(self$testXXX( idx , method ),silent=TRUE)
			str_idx = gsub(", ","",toString(idx))
			if( "try-error" %in% class(test) )
			{
				if( base::min(idx) == 2 )
				{
					row = base::c( paste0("Test ",str_idx),"OK", "/" , "/" , "/" )
				}
				else
				{
					row = base::c( paste0("Test ",str_idx),"Fail", "/" , "/" , "/" )
				}
			}
			else
			{
				row = base::c( paste0("Test ",str_idx),"OK",round(max(abs(self$coef_ - self$law$coef_)),3) , toString(self$coef_) , toString(round(self$law$coef_,2)) )
			}
			tab$add_row( row )
		}
		
		cat(tab$draw())
	}
	##}}}
	
	)
	
)
##}}}

NormalTest = R6::R6Class( "NormalTest" , ##{{{
	
	inherit = SDFCLawTest,
	
	public = list(
	
	initialize = function( n_samples = 2500 )
	{
		kwargs = list( n_samples = n_samples,
		               name    = "Normal",
		               sd_law  = SDFC::Normal,
		               params  = base::c("loc","scale"),
		               shape_p = FALSE
		)
		base::do.call( super$initialize , kwargs )
	},
	
	rvs = function()
	{
		return(stats::rnorm( length(self$loc) , self$loc , self$scale ))
	}
	
	)
)
##}}}


###############
## Functions ##
###############

##########
## main ##
##########

for( nt in list(NormalTest) )
{
	t = nt$new()
	t$run_all( "MLE" )
	base::cat("\n")
	t$run_all( "bayesian" )
	base::cat("\n")
}


## End
plt$wait()
base::cat("Done\n")

