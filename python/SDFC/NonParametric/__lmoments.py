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

import numpy         as np
import scipy.special as scs
from SDFC.NonParametric.__quantile import quantile


###############
## Functions ##
###############

def _lmoments_stationary( Y ):##{{{
	Ys = np.sort(Y.squeeze())
	lmom = np.zeros(4)
	
	## Order 2
	C0 = scs.binom( range( Y.size ) , 1 )
	C1 = scs.binom( range( Y.size - 1 , -1 , -1 ) , 1 )
	
	## Order 3
	C2 = scs.binom( range( Y.size ) , 2 )
	C3 = scs.binom( range( Y.size - 1 , -1 , -1 ) , 2 )
	
	## Order 4
	C4 = scs.binom( range( Y.size ) , 3 )
	C5 = scs.binom( range( Y.size - 1 , -1 , -1 ) , 3 )
	
	
	lmom[0] = np.mean(Ys)
	lmom[1] = np.sum( ( C0 - C1 ) * Ys ) / ( 2 * scs.binom( Y.size , 2 ) )
	lmom[2] = np.sum( ( C2 - 2 * C0 * C1 + C3 ) * Ys ) / ( 3 * scs.binom( Y.size , 3 ) )
	lmom[3] = np.sum( (C4 - 3 * C2 * C1 + 3 * C0 * C3 - C5 ) * Ys ) / ( 4 * scs.binom( Y.size , 4 ) )
	
	return lmom
##}}}

def lmoments( Y , c_Y = None , order = None , lq = np.arange( 0.05 , 0.96 , 0.01 ) ):##{{{
	"""
		SDFC.NonParametric.lmoments
		===========================
		
		Estimate the lmoments of orders 1 to 4. If a covariate is given, a quantile regression is performed
		and the instantaneous L-Moments are estimated from the quantile fitted.
		
		Parameters
		----------
		Y     : np.array
			Dataset to fit the lmoments
		c_Y   : np.array or None
			Covariate
		order : integer, list of integer or None
			Integers between 1 and 4
		lq    : np.array
			Quantiles for quantile regression, only used if a covariate is given. Default is np.arange(0.05,0.96,0.01)
		
		Returns
		-------
		The lmoments.
	"""
	
	order = order if order is None else np.array( [order] , dtype = np.int ).squeeze() - 1
	
	if c_Y is None:
		lmom = _lmoments_stationary(Y)
		return lmom if order is None else lmom[order]
	else:
		Y = Y.reshape(-1,1)
		if c_Y.ndim == 1: c_Y = c_Y.reshape(-1,1)
		Yq = quantile( Y , lq , c_Y )
		lmom = np.zeros((Y.size,4))
		for i in range(Y.size):
			lmom[i,:] = _lmoments_stationary(Yq[i,:])
		if order is None:
			return lmom
		else:
			return lmom[:,order]
##}}}


