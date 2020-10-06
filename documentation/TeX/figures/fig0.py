# -*- coding: utf-8 -*-

#############################
## Yoann Robin             ##
## yoann.robin.k@gmail.com ##
#############################

###############
## Libraries ##
##{{{

import sys,os
import pickle as pk
import multiprocessing as mp


## Scientific libraries
##=====================

import numpy as np
import scipy.stats as sc
import pandas as pd
import xarray as xr


## Plot libraries ##
##==================

import matplotlib as mpl
import matplotlib.pyplot as plt
## from mpl_toolkits.mplot3d import Axes3D
## import cartopy.crs as ccrs
## import cartopy.feature as cfeature

mpl.rcParams['font.size'] = 15
#plt.rc('text',usetex=True)
#plt.rcParams['text.latex.unicode'] = True

##}}}
###############

###############
## Fonctions ##
###############


#############
## Classes ##
#############


##########
## main ##
##########

if __name__ == "__main__":
	
	X = np.random.normal( size = 100000 )
	bins = np.linspace( X.min() - 1 , X.max() + 1 , 100 )
	kde  = sc.gaussian_kde(X)
	
	fig = plt.figure( figsize = (10,10 ) )
	
	ax = fig.add_subplot(1,1,1)
	ax.hist( X , bins = bins , color = "blue" , alpha = 0.5 , density = True , label = "Empirical histogram" )
	ax.plot( bins , kde(bins) , color = "blue" , label = "Gaussian kernel" )
	ax.legend( loc = "upper left" )
	ax.set_xlabel( r"$x$" )
	ax.set_ylabel( "Density" )
	
	fig.set_tight_layout(True)
	fig.savefig( os.path.join( os.path.dirname(os.path.abspath(__file__)) , "fig0.png" ) )
	
	print("Done")
