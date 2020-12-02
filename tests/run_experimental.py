# -*- coding: utf-8 -*-

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


###############
## Libraries ##
##{{{

import sys,os
import texttable as tt

## Scientific libraries
##=====================

import numpy as np
import scipy.stats as sc
import scipy.optimize as sco
import scipy.linalg as scl
import pandas as pd
import xarray as xr
import SDFC as sd
import SDFC.NonParametric as sdnp

## Plot libraries ##
##==================

import matplotlib as mpl
import matplotlib.pyplot as plt

#mpl.rcParams['font.size'] = 30

##}}}
###############

from experimental.core.__RHS import LHS
from experimental.core.__RHS import RHS
from experimental.link.__Univariate import ULIdentity
from experimental.link.__Univariate import ULExponential
from experimental.link.__Multivariate import MultivariateLink
from experimental.link.__Multivariate import MLLinear
from experimental.link.__Multivariate import MLTensor
from experimental.link.__Normal import NormalRatioLocScaleConstant
from experimental.link.__GEV    import    GEVRatioLocScaleConstant


from experimental import Normal
from experimental import GEV


###############
## Fonctions ##
###############


#############
## Classes ##
#############


## GEV part
##=========



## Tests
##======

class NormalTest: ##{{{
	
	def __init__( self , n_sample = 2000 ): ##{{{
		self.n_samples     = n_sample
		t,X_loc,X_scale,_ = sd.tools.Dataset.covariates(self.n_samples)
		self.t       = t
		self.X_loc   = X_loc.reshape(-1,1)
		self.X_scale = X_scale.reshape(-1,1)
	##}}}
	
	def test0( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		
		kwargs = { "c_loc" : self.X_loc , "c_scale" : self.X_scale , "l_scale" : ULExponential() }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def test1( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		
		l_global = MLTensor( [MLLinear( c = self.X_loc , l = ULIdentity() ) , MLLinear( c = self.X_scale , l = ULExponential() ) ] , [2,2] , n_samples = self.n_samples , n_features = 4 )
		kwargs = { "c_global" : [self.X_loc,self.X_scale] , "l_global" : l_global }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def test2( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		self.coef_ = self.coef_[:2]
		
		kwargs = { "c_loc" : self.X_loc , "f_scale" : self.scale }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def test3( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		self.coef_ = self.coef_[2:]
		
		kwargs = { "f_loc" : self.loc , "c_scale" : self.X_scale , "l_scale" : ULExponential() }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def test4( self , method = "MLE" ):##{{{
		self.coef_  = np.array([0.8,1.5,2])
		l_global    = NormalRatioLocScaleConstant( self.n_samples )
		self.loc,self.scale = l_global.transform( self.coef_ , self.X_loc )
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		
		kwargs = { "c_global" : [self.X_loc] , "l_global" : l_global }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def test5( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,0.3,-0.9] )
		self.loc   = np.repeat( self.coef_[0] , self.n_samples ).reshape(-1,1)
		self.scale = np.exp(self.coef_[1] + self.coef_[2] * self.X_scale)
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		
		kwargs = { "c_scale" : self.X_scale , "l_scale" : ULExponential() }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def test6( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.3,-0.9] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.repeat( np.exp(self.coef_[2]) , self.n_samples ).reshape(-1,1)
		self.Y     = np.random.normal( loc = self.loc , scale = self.scale )
		
		kwargs = { "c_loc" : self.X_loc , "l_scale" : ULExponential() }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = np.identity(self.coef_.size) )
		self.norm = Normal( method = method )
		self.norm.fit( self.Y , **kwargs )
	##}}}
	
	def summary( self , show = False ): ##{{{
		print( "## => {} / {} / {}".format( np.max(np.abs(self.coef_ - self.norm.coef_)) , self.coef_ , self.norm.coef_ ) )
	##}}}
	
	def run_all( self , method = "MLE" , show = True ):##{{{
		tab = tt.Texttable( max_width = 0 )
		tab.header( ["Normal law test ({})".format(method),"Status","Max diff","True value","Estimated value"] )
		for i in range(7):
			try:
				eval( "self.test{}( method = \"{}\" )".format(i,method) )
				tab.add_row( ["Test {}".format(i),"OK",np.max(np.abs(self.coef_ - self.norm.coef_)) , self.coef_ , np.round(self.norm.coef_,2)] )
			except:
				tab.add_row( ["Test {}".format(i),"Fail","/","/","/"] )
		
		if show: print(tab.draw())
		return tab
	##}}}
	
##}}}

class GEVTest:##{{{
	
	def __init__( self , n_sample = 2000 ): ##{{{
		self.n_samples     = n_sample
		t,X_loc,X_scale,X_shape = sd.tools.Dataset.covariates(self.n_samples)
		self.t       = t
		self.X_loc   = X_loc.reshape(-1,1)
		self.X_scale = X_scale.reshape(-1,1)
		self.X_shape = X_shape.reshape(-1,1)
	##}}}
	
	def test0( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9,-0.2] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.shape = np.repeat( self.coef_[4] , self.n_samples ).reshape(-1,1)
		self.Y     = sc.genextreme.rvs( loc = self.loc , scale = self.scale , c = - self.shape )
		
		kwargs = { "c_loc" : self.X_loc , "c_scale" : self.X_scale , "l_scale" : ULExponential() }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = 0.1 * np.identity(self.coef_.size) )
		self.gev = GEV( method = method )
		self.gev.fit( self.Y , **kwargs )
	##}}}
	
	def test1( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9,-0.2] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.shape = np.repeat( self.coef_[4] , self.n_samples ).reshape(-1,1)
		self.Y     = sc.genextreme.rvs( loc = self.loc , scale = self.scale , c = - self.shape )
		self.coef_ = np.array( [0.3,-0.9,-0.2] )
		
		kwargs = { "f_loc" : self.loc , "c_scale" : self.X_scale , "l_scale" : ULExponential() }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = 0.1 * np.identity(self.coef_.size) )
		self.gev = GEV( method = method )
		self.gev.fit( self.Y , **kwargs )
	##}}}
	
	def test2( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9,-0.2] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.shape = np.repeat( self.coef_[4] , self.n_samples ).reshape(-1,1)
		self.Y     = sc.genextreme.rvs( loc = self.loc , scale = self.scale , c = - self.shape )
		self.coef_ = np.array( [0.5,1.,-0.2] )
		
		kwargs = { "c_loc" : self.X_loc , "f_scale" : self.scale }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = 0.1 * np.identity(self.coef_.size) )
		self.gev = GEV( method = method )
		self.gev.fit( self.Y , **kwargs )
	##}}}
	
	def test3( self , method = "MLE" ):##{{{
		self.coef_ = np.array( [0.5,1.,0.3,-0.9,-0.2] )
		self.loc   = self.coef_[0] + self.coef_[1] * self.X_loc
		self.scale = np.exp(self.coef_[2] + self.coef_[3] * self.X_scale)
		self.shape = np.repeat( self.coef_[4] , self.n_samples ).reshape(-1,1)
		self.Y     = sc.genextreme.rvs( loc = self.loc , scale = self.scale , c = - self.shape )
		self.coef_ = np.array( [0.5,1.,0.3,-0.9] )
		
		kwargs = { "c_loc" : self.X_loc , "c_scale" : self.X_scale , "l_scale" : ULExponential() , "f_shape" : self.shape }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = 0.1 * np.identity(self.coef_.size) )
		self.gev = GEV( method = method )
		self.gev.fit( self.Y , **kwargs )
	##}}}
	
	def test4( self , method = "MLE" ):##{{{
		self.coef_  = np.array([0.8,1.5,2,-0.2])
		l_global    = GEVRatioLocScaleConstant( self.n_samples )
		self.loc,self.scale,self.shape = l_global.transform( self.coef_ , self.X_loc )
		self.Y     = sc.genextreme.rvs( loc = self.loc , scale = self.scale , c = - self.shape )
		
		kwargs = { "c_global" : [self.X_loc] , "l_global" : l_global }
		kwargs["prior"] = sc.multivariate_normal( mean = self.coef_ , cov = 0.1 * np.identity(self.coef_.size) )
		self.gev = GEV( method = method )
		self.gev.fit( self.Y , **kwargs )
	##}}}
	
	def summary( self , show = False ): ##{{{
		print( "## => {} / {} / {}".format( np.max(np.abs(self.coef_ - self.gev.coef_)) , self.coef_ , np.round(self.gev.coef_,3) ) )
	##}}}
	
	def run_all( self , method = "MLE" , show = True ):##{{{
		tab = tt.Texttable( max_width = 0 )
		tab.header( ["GEV law test ({})".format(method),"Status","Max diff","True value","Estimated value"] )
		for i in range(5):
			try:
				eval( "self.test{}( method = \"{}\" )".format(i,method) )
				tab.add_row( ["Test {}".format(i),"OK",np.max(np.abs(self.coef_ - self.gev.coef_)) , self.coef_ , np.round(self.gev.coef_,2)] )
			except:
				tab.add_row( ["Test {}".format(i),"Fail","/","/","/"] )
		
		if show: print(tab.draw())
		return tab
	##}}}
	
##}}}

##########
## main ##
##########

if __name__ == "__main__":
	np.seterr( all = "ignore" )
#	np.random.seed(42)
	
#	nt = NormalTest()
#	nt.run_all("MLE")
#	nt.run_all("bayesian")
	
	gt = GEVTest()
	gt.run_all("MLE")
	gt.run_all("bayesian")
	
	
	print("Done")



