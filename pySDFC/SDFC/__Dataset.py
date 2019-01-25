# -*- coding: utf-8 -*-

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## This software is a computer program that is part of the SDFC (Statistical    ##
## Distribution Fit with Covariates) library. This library makes it possible    ##
## to regress the parameters of some statistical law with co-variates.          ##
##                                                                              ##
## This software is governed by the CeCILL-C license under French law and       ##
## abiding by the rules of distribution of free software.  You can  use,        ##
## modify and/ or redistribute the software under the terms of the CeCILL-C     ##
## license as circulated by CEA, CNRS and INRIA at the following URL            ##
## "http://www.cecill.info".                                                    ##
##                                                                              ##
## As a counterpart to the access to the source code and  rights to copy,       ##
## modify and redistribute granted by the license, users are provided only      ##
## with a limited warranty  and the software's author,  the holder of the       ##
## economic rights,  and the successive licensors  have only  limited           ##
## liability.                                                                   ##
##                                                                              ##
## In this respect, the user's attention is drawn to the risks associated       ##
## with loading,  using,  modifying and/or developing or reproducing the        ##
## software by the user in light of its specific status of free software,       ##
## that may mean  that it is complicated to manipulate,  and  that  also        ##
## therefore means  that it is reserved for developers  and  experienced        ##
## professionals having in-depth computer knowledge. Users are therefore        ##
## encouraged to load and test the software's suitability as regards their      ##
## requirements in conditions enabling the security of their systems and/or     ##
## data to be ensured and,  more generally, to use and operate it in the        ##
## same conditions as regards security.                                         ##
##                                                                              ##
## The fact that you are presently reading this means that you have had         ##
## knowledge of the CeCILL-C license and that you accept its terms.             ##
##                                                                              ##
##################################################################################
##################################################################################

##################################################################################
##################################################################################
##                                                                              ##
## Copyright Yoann Robin, 2019                                                  ##
##                                                                              ##
## yoann.robin.k@gmail.com                                                      ##
##                                                                              ##
## Ce logiciel est un programme informatique faisant partie de la librairie     ##
## SDFC (Statistical Distribution Fit with Covariates). Cette librairie         ##
## permet de calculer de regresser les parametres de lois statistiques selon    ##
## plusieurs co-variables                                                       ##
##                                                                              ##
## Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    ##
## respectant les principes de diffusion des logiciels libres. Vous pouvez      ##
## utiliser, modifier et/ou redistribuer ce programme sous les conditions       ##
## de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     ##
## sur le site "http://www.cecill.info".                                        ##
##                                                                              ##
## En contrepartie de l'accessibilité au code source et des droits de copie,    ##
## de modification et de redistribution accordés par cette licence, il n'est    ##
## offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    ##
## seule une responsabilité restreinte pèse sur l'auteur du programme, le       ##
## titulaire des droits patrimoniaux et les concédants successifs.              ##
##                                                                              ##
## A cet égard  l'attention de l'utilisateur est attirée sur les risques        ##
## associés au chargement,  à l'utilisation,  à la modification et/ou au        ##
## développement et à la reproduction du logiciel par l'utilisateur étant       ##
## donné sa spécificité de logiciel libre, qui peut le rendre complexe à        ##
## manipuler et qui le réserve donc à des développeurs et des professionnels    ##
## avertis possédant  des  connaissances  informatiques approfondies.  Les      ##
## utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       ##
## logiciel à leurs besoins dans des conditions permettant d'assurer la         ##
## sécurité de leurs systèmes et ou de leurs données et, plus généralement,     ##
## à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           ##
##                                                                              ##
## Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       ##
## pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    ##
## termes.                                                                      ##
##                                                                              ##
##################################################################################
##################################################################################


###############
## Libraries ##
###############

import numpy as np


#############
## Classes ##
#############

class Dataset:
	"""
	SDFC.Dataset
	============
	
	Dataset to test the quantile regression
	
	Attributes
	----------
	
	t  : numpy.array
		Time axis
	X  : numpy.array
		Co-variable
	Y  : numpy.array
		Dataset to fit
	Yq : numpy.array
		Quantile regression
	tic : float
		time just at the end of the method __init__
	tac : float
		time at the begin of the method quick_plot
	time_execution: float
		tac-tic
	
	example
	-------
	
	import SDFC as sd
	
	## Define data
	size = 1000
	kind = 0
	use_phi = False
	data = sd.Dataset(size,kind)
	
	## Fit with a normal law
	norm = sd.NormalLaw( use_phi = use_phi )
	norm.fit( data.Y , loc_cov = data.X , scale_cov = data.X )
	Y_norm = np.random.normal( loc = norm.loc , scale = norm.scale )
	
	## Quantile Regression
	ltau = np.arange( 0.05 , 0.96 , 0.01 )
	qr = sd.QuantileRegression( ltau )
	qr.fit( data.Y , data.X )
	if qr.is_success():
		quants = qr.predict()
	
	## GDP for extremes
	gpdU = sd.GPDLaw( use_phi = use_phi )
	gpdU.fit( data.Y , loc = quants[:,-1] , scale_cov = data.X )
	Y_upper = sc.genpareto.rvs( loc = quants[:,-1].ravel() , scale = gpdU.scale , c = gpdU.shape )
	
	gpdL = sd.GPDLaw( use_phi = use_phi )
	gpdL.fit( -data.Y , loc = -quants[:,0] , scale_cov = data.X )
	Y_lower = -sc.genpareto.rvs( loc = -quants[:,0].ravel() , scale = gpdL.scale , c = gpdL.shape )
	
	## Plot
	nrow,ncol = 2,2
	figscale  = 4
	yLim = (np.min( (data.Y,Y_norm,Y_lower) ),np.max( (data.Y,Y_norm,Y_lower) ))
	fig = plt.figure( figsize = ( 1.6 * figscale * ncol , figscale * nrow ) )
	
	ax = fig.add_subplot( nrow , ncol , 1 )
	ax.plot( data.t , data.Y , color = "blue" , linestyle = "" , marker = "." , alpha = 0.5 )
	ax.set_ylim( yLim )
	ax.set_title( "Original data" )
	
	ax = fig.add_subplot( nrow , ncol , 2 )
	ax.plot( data.t , Y_norm , color = "blue" , linestyle = "" , marker = "." , alpha = 0.5 )
	ax.set_ylim( yLim )
	ax.set_title( "Normal law fit" )
	
	ax = fig.add_subplot( nrow , ncol , 3 )
	for i in range(ltau.size):
		ax.plot( data.t , quants[:,i] , color = "black" , linestyle = "-" )
	ax.plot( data.t , Y_lower , color = "red" , linestyle = "" , marker = "." , alpha = 0.5 )
	ax.plot( data.t , Y_upper , color = "red" , linestyle = "" , marker = "." , alpha = 0.5 )
	ax.set_ylim( yLim )
	ax.set_title( "GPD and QuantileRegression fit" )
	
	plt.tight_layout()
	plt.show()
	
	"""
	
	
	def __init__( self , size = 500 , kind = 0 , alpha = 0.05 , beta = 0.6 ): ##{{{
		"""
		Initialization of Dataset
		
		Arguments
		---------
		size : int
			Length of data
		kind : int
			=> if 0 : Y is a Gaussian law where location move with time, but scale is fixed at 0.2, X is the location parameter
			=> if 1 : Y is a Gaussian law where location and scale move with time, X is the location parameter
			=> if 2 : Y is a Gaussian law where location and scale move with time, X is the location and the scale parameter
		"""
		self.size = size
		self.alpha = alpha
		self.beta = beta
		self.t = None
		self.X = None
		self.Y = None
		self.Yq = None
		
		if kind == 0:
			self.kind0()
		elif kind == 1:
			self.kind1()
		elif kind == 2:
			self.kind2()
		elif kind == 3:
			self.kind3()
		else:
			pass
		
	##}}}
	
	def kind0(self): ##{{{
		self.t = np.linspace( 0 , 1 , self.size )
		self.X = np.zeros( (self.size,1) )
		self.X[:,0] = self.t**2 + np.cos( 2* np.pi * self.t ) * 0.2
		self.Y = self.X[:,0] + np.random.normal( loc = 0 , scale = 0.1 , size = self.size )
	##}}}
	
	def kind1(self): ##{{{
		self.t = np.linspace( 0 , 1 , self.size )
		self.X = np.zeros( (self.size,1) )
		self.X[:,0] = (self.beta - self.alpha) * np.cos( 2* np.pi * self.t ) / 2.
		sigma = self.X[:,0] + ( self.beta + self.alpha ) / 2.
		self.X[:,0] += self.t**2
		self.Y = self.X[:,0] + np.random.normal( loc = np.zeros(self.size) , scale = sigma )
	##}}}
	
	def kind2(self): ##{{{
		self.t = np.linspace( 0 , 1 , self.size )
		self.X = np.zeros( (self.size,2) )
		self.X[:,0] = (self.beta - self.alpha) * np.cos( 2* np.pi * self.t ) / 2.
		self.X[:,1] = self.X[:,0] + ( self.beta + self.alpha ) / 2.
		self.X[:,0] += self.t**2
		self.Y = self.X[:,0] + np.random.normal( loc = np.zeros(self.size) , scale = self.X[:,1] )
	##}}}
	
	def kind3(self): ##{{{
		self.t = np.linspace( 0 , 1 , self.size )
		self.X = np.zeros( (self.size,1) )
		self.Y = np.random.normal( loc = 0 , scale = 0.2 , size = self.size )
	##}}}
	
