
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

#' np_lmoments_stationary
#'
#' Compute the four first L-Moments
#'
#' @param Y  [vector] Dataset
#'
#' @return [lmom] L-Moments
#'
#' @examples
#' ## Data
#' Y    = stats::rexp( n = 10000 )
#' lmom = SDFC::np_lmoments_stationary(Y)
#' @export
np_lmoments_stationary = function(Y)
{
	Ys = Y[order(Y)]
	size = length(Ys)
	lmom = numeric(4)
	
	## Order 2
	C0 = base::choose( 1:size , 1 )
	C1 = base::choose( base::seq( size - 1 , 0 , -1 ) , 1 )
	
	## Order 3
	C2 = base::choose( 1:size , 2 )
	C3 = base::choose( base::seq( size - 1 , 0 , -1 ) , 2 )
	
	## Order 4
	C4 = base::choose( 1:size , 3 )
	C5 = base::choose( base::seq( size - 1 , 0 , -1 ) , 3 )
	
	
	lmom[1] = base::mean(Ys)
	lmom[2] = base::sum( ( C0 - C1 ) * Ys ) / ( 2 * base::choose( size , 2 ) )
	lmom[3] = base::sum( ( C2 - 2 * C0 * C1 + C3 ) * Ys ) / ( 3 * base::choose( size , 3 ) )
	lmom[4] = base::sum( ( C4 - 3 * C2 * C1 + 3 * C0 * C3 - C5 ) * Ys ) / ( 4 * base::choose( size , 4 ) )
	
	return(lmom)
}


#' np_lmoments
#'
#' Compute the L-Moments
#'
#' @param Y  [vector] Dataset
#' @param c_Y  [vector] Covariate. If NULL stationary L-moments are computed
#' @param order [int] order of moments. a vector with elements between 1 and 4
#' @param lq [vector] Vector of quantile used for quantile regression. Default is seq(0.05,0.95,0.01)
#'
#' @return [lmom] L-Moments
#'
#' @examples
#' ## Data
#' size = 2000
#' c_data = Dataset$covariates(size)
#' 
#' t       = c_data$t
#' X_loc   = c_data$X_loc
#' X_scale = c_data$X_scale
#' loc   = 0.5 + 2 * X_loc
#' scale =   1 + 2 * X_scale
#' Y = stats::rnorm( size , mean = loc , sd = scale )
#' 
#' c_Y = base::cbind( X_loc , X_scale )
#' 
#' lmom = np_lmoments( Y , c_Y = c_Y )
#' @export
np_lmoments = function( Y , c_Y = NULL , order = NULL , lq = base::seq( 0.05 , 0.95 , 0.01 ) )
{
	if( is.null(order) )
		order = 1:4
	
	if( is.null(c_Y) )
	{
		lmom = np_lmoments_stationary(Y)
		return(lmom[order])
	}
	
	if( !is.matrix(c_Y) )
		c_Y = matrix( c_Y , nrow = length(c_Y) , ncol = 1 )
	
	Yq = np_quantile( Y , lq , c_Y )
	lmom = base::t( base::apply( Yq , 1 , np_lmoments_stationary ) )
	
	return( lmom[,order] )
}



