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


###########
## Class ##
###########

def plot_confidences_intervals( law , ax , color = "red" ):##{{{
	"""
	SDFC.tools.plot_confidence_intervals
	====================================
	Plot confidence intervals of a law if fitted
	
	Parameters
	----------
	
	law : SDFC law
		A law fitted by SDFC
	ax  : matplotlib axes
		Axes where figure is drawn
	color : matplotlib color
		Color used, default is red
	"""
	
	if law.confidence_interval is None or law.coef_ is None:
		return
	
	coef   = law.coef_
	n_coef = law.coef_.size
	
	i = 0
	lticks = []
	for p in law._lparams:
		if p.not_fixed():
			n_p = p.coef_.size
			for j in range(n_p):
				xb = i - 0.3
				xe = i + 0.3
				ax.fill_between( [xb,xe] , p.linkFct(law.confidence_interval[0,i] - coef[i]) , p.linkFct(law.confidence_interval[1,i] - coef[i]) , color = color , alpha = 0.5 )
				ax.text( xb , 0 , "{}".format(round(p.linkFct(coef[i]),2)) )
				lticks.append( "{}{}".format( p.kind , j ) )
				i += 1
	ax.hlines( 0 , -0.5 , n_coef - 0.5 , color = "black" )
	ax.set_xlim( (-0.5,n_coef-0.5) )
	ax.set_xticks( range(n_coef) )
	ax.set_xticklabels( lticks )
	ax.set_ylabel( "Parameters values" )
	ax.set_xlabel( "Parameters" )
	
	return ax
##}}}

def plot_confidences_intervals_scaled( law , ax , color = "red" ):##{{{
	"""
	SDFC.tools.plot_confidence_intervals_scaled
	===========================================
	Plot confidence intervals of a law if fitted, where values of all parameters are scaled into the same range.
	
	Parameters
	----------
	
	law : SDFC law
		A law fitted by SDFC
	ax  : matplotlib axes
		Axes where figure is drawn
	color : matplotlib color
		Color used, default is red
	"""
	
	if law.confidence_interval is None or law.coef_ is None:
		return
	
	coef   = law.coef_
	n_coef = law.coef_.size
	
	i = 0
	lticks = []
	for p in law._lparams:
		if p.not_fixed():
			n_p = p.coef_.size
			for j in range(n_p):
				xb = i - 0.3
				xe = i + 0.3
				val  = p.linkFct(coef[i])
				valu = p.linkFct(law.confidence_interval[1,i] - coef[i])
				vall = p.linkFct(law.confidence_interval[0,i] - coef[i])
				valu -= p.linkFct(0)
				vall -= p.linkFct(0)
				exp  = -np.log10( max(valu,abs(vall)) )
				exp  = 0 if not np.isfinite(exp) else int(exp)
				ax.fill_between( [xb,xe] , vall * 10**exp, valu * 10**exp , color = color , alpha = 0.5 )
				txt  = r"$" + r"{}".format(round(val * 10**exp,2))
				if exp != 0:
					txt += r"(\times 10^{" + str(-exp) + r"})$"
				else:
					txt += r"$"
				ax.text( xb , 0 , txt )
				lticks.append( "{}{}".format( p.kind , j ) )
				i += 1
	ax.hlines( 0 , -0.5 , n_coef - 0.5 , color = "black" )
	ax.set_xlim( (-0.5,n_coef-0.5) )
	ax.set_xticks( range(n_coef) )
	ax.set_xticklabels( lticks )
	ax.set_ylabel( "Parameters values" )
	ax.set_xlabel( "Parameters" )
	
	return ax
##}}}


