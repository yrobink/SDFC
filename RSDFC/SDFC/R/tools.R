
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


#' lmoments1
#'
#' Compute the L-Moments of order 1 (just the mean...)
#'
#' @param Y  [vector] Dataset
#'
#' @return [lmom1] L-Moments of order 1
#'
#' @examples
#' ## Data
#' size = 2000
#' data = Dataset0(size)
#' lmom1 = SDFC::lmoments1(data$Y)
#' @export
lmoments1 = function(Y)
{
	return( base::mean(Y) )
}


#' lmoments2
#'
#' Compute the L-Moments of order 2 (half mean of pairwise difference)
#'
#' @param Y  [vector] Dataset
#'
#' @return [lmom2] L-Moments of order 2
#'
#' @examples
#' ## Data
#' size = 2000
#' data = Dataset0(size)
#' lmom2 = SDFC::lmoments2(data$Y)
#' @export
lmoments2 = function(Y)
{
	size = length(Y)
	res = 0
	
	X = Y[order(Y)]
	
	for( i in 1:size )
	{
		res = res + ( base::choose( i - 1 , 1 ) - base::choose( size - i , 1 ) ) * X[i]
	}
	res = res / ( 2 * base::choose( size , 2 ) )
	
	return(res)
}


#' lmoments3
#'
#' Compute the L-Moments of order 3 
#'
#' @param Y  [vector] Dataset
#'
#' @return [lmom3] L-Moments of order 3
#'
#' @examples
#' ## Data
#' size = 2000
#' data = Dataset0(size)
#' lmom2 = SDFC::lmoments3(data$Y)
#' @export
lmoments3 = function(Y)
{
	size = length(Y)
	res = 0
	
	X = Y[order(Y)]
	
	for( i in 1:size )
	{
		res = res + ( base::choose( i - 1 , 2 ) - 2 * base::choose( i-1 , 1 ) * base::choose( size - i , 1 ) + base::choose( size - i , 2 ) ) * X[i]
	}
	res = res / ( 3 * base::choose( size , 3 ) )
	
	return(res)
}



#' lmoments
#'
#' Compute the L-Moments
#'
#' @param Y  [vector] Dataset
#'
#' @param order [int] order of moments
#'
#' @return [lmom] L-Moments
#'
#' @examples
#' ## Data
#' size = 2000
#' data = Dataset0(size)
#' lmom1 = SDFC::lmoments(data$Y,1)
#' lmom2 = SDFC::lmoments(data$Y,2)
#' lmom3 = SDFC::lmoments(data$Y,3)
#' @export
lmoments = function( Y , order )
{
	if( order == 1 )
		return( SDFC::lmoments1(Y) )
	if( order == 2 )
		return( SDFC::lmoments2(Y) )
	if( order == 3 )
		return( SDFC::lmoments3(Y) )
	return(NULL)
}



