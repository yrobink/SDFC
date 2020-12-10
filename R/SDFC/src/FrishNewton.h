
/********************************************************************************/
/********************************************************************************/
/*                                                                              */
/* Copyright Yoann Robin, 2019                                                  */
/*                                                                              */
/* yoann.robin.k@gmail.com                                                      */
/*                                                                              */
/* This software is a computer program that is part of the SDFC (Statistical    */
/* Distribution Fitted with Covariates) library. This library makes it possible */
/* to regress the parameters of some statistical law with co-variates.          */
/*                                                                              */
/* This software is governed by the CeCILL-C license under French law and       */
/* abiding by the rules of distribution of free software.  You can  use,        */
/* modify and/ or redistribute the software under the terms of the CeCILL-C     */
/* license as circulated by CEA, CNRS and INRIA at the following URL            */
/* "http://www.cecill.info".                                                    */
/*                                                                              */
/* As a counterpart to the access to the source code and  rights to copy,       */
/* modify and redistribute granted by the license, users are provided only      */
/* with a limited warranty  and the software's author,  the holder of the       */
/* economic rights,  and the successive licensors  have only  limited           */
/* liability.                                                                   */
/*                                                                              */
/* In this respect, the user's attention is drawn to the risks associated       */
/* with loading,  using,  modifying and/or developing or reproducing the        */
/* software by the user in light of its specific status of free software,       */
/* that may mean  that it is complicated to manipulate,  and  that  also        */
/* therefore means  that it is reserved for developers  and  experienced        */
/* professionals having in-depth computer knowledge. Users are therefore        */
/* encouraged to load and test the software's suitability as regards their      */
/* requirements in conditions enabling the security of their systems and/or     */
/* data to be ensured and,  more generally, to use and operate it in the        */
/* same conditions as regards security.                                         */
/*                                                                              */
/* The fact that you are presently reading this means that you have had         */
/* knowledge of the CeCILL-C license and that you accept its terms.             */
/*                                                                              */
/********************************************************************************/
/********************************************************************************/

/********************************************************************************/
/********************************************************************************/
/*                                                                              */
/* Copyright Yoann Robin, 2019                                                  */
/*                                                                              */
/* yoann.robin.k@gmail.com                                                      */
/*                                                                              */
/* Ce logiciel est un programme informatique faisant partie de la librairie     */
/* SDFC (Statistical Distribution Fitted with Covariates). Cette librairie      */
/* permet de calculer de regresser les parametres de lois statistiques selon    */
/* plusieurs co-variables                                                       */
/*                                                                              */
/* Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    */
/* respectant les principes de diffusion des logiciels libres. Vous pouvez      */
/* utiliser, modifier et/ou redistribuer ce programme sous les conditions       */
/* de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     */
/* sur le site "http://www.cecill.info".                                        */
/*                                                                              */
/* En contrepartie de l'accessibilité au code source et des droits de copie,    */
/* de modification et de redistribution accordés par cette licence, il n'est    */
/* offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    */
/* seule une responsabilité restreinte pèse sur l'auteur du programme, le       */
/* titulaire des droits patrimoniaux et les concédants successifs.              */
/*                                                                              */
/* A cet égard  l'attention de l'utilisateur est attirée sur les risques        */
/* associés au chargement,  à l'utilisation,  à la modification et/ou au        */
/* développement et à la reproduction du logiciel par l'utilisateur étant       */
/* donné sa spécificité de logiciel libre, qui peut le rendre complexe à        */
/* manipuler et qui le réserve donc à des développeurs et des professionnels    */
/* avertis possédant  des  connaissances  informatiques approfondies.  Les      */
/* utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       */
/* logiciel à leurs besoins dans des conditions permettant d'assurer la         */
/* sécurité de leurs systèmes et ou de leurs données et, plus généralement,     */
/* à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           */
/*                                                                              */
/* Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       */
/* pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    */
/* termes.                                                                      */
/*                                                                              */
/********************************************************************************/
/********************************************************************************/

/********************************************************************************/
/* This software is based on the fortran code developped by Roger Koenker,      */
/* coming from the R package "quantreg", see:                                   */
/* https://cran.r-project.org/web/packages/quantreg/index.html                  */
/********************************************************************************/

/********************************************************************************/
/* Ce programme est basé sur le code fortran développé par Roger Koenker,       */
/* venant du package R "quantreg", voir:                                        */
/* https://cran.r-project.org/web/packages/quantreg/index.html                  */
/********************************************************************************/


#ifndef SDFC_QUANTILEREGRESSION_FRISHNEWTON_INCLUDED
#define SDFC_QUANTILEREGRESSION_FRISHNEWTON_INCLUDED


//===========//
// Libraries //
//===========//


#include <random>
#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Core>

#include "QRState.h"
#include "FrishNewtonCore.h"

namespace QuantileRegression
{

//===========//
// Class(es) //
//===========//

//' @export
class FrishNewton //{{{
{
	public:
	
	//=========//
	// Typedef // 
	//=========//
	
	//{{{
	typedef unsigned int    size_type  ;
	typedef double          value_type ;
	typedef Eigen::ArrayXd  Array      ;
	typedef Eigen::MatrixXd Matrix     ;
	typedef tools::qrstate_t qrstate_t ;
	//}}}
	
	
	//=============//
	// Constructor //
	//=============//
	
	FrishNewton( size_type maxit = 50 , value_type tol = 1e-6 , value_type beta = 0.99995 ): //{{{
		m_ltau() ,
		m_coef() ,
		m_frishNewtonCore( 0.5 , maxit , tol , beta ),
		m_state(tools::not_fitted)
	{} //}}}
	
	FrishNewton( Array ltau , size_type maxit = 50 , value_type tol = 1e-6 , value_type beta = 0.99995 ): //{{{
		m_ltau(ltau) ,
		m_coef() ,
		m_frishNewtonCore( 0.5 , maxit , tol , beta ),
		m_state(tools::not_fitted)
	{} //}}}
	
	~FrishNewton() //{{{
	{} //}}}
	
	
	//===========//
	// Accessors //
	//===========//
	
	void set_ltau( Array ltau ) //{{{
	{
		m_ltau = ltau ;
	} //}}}
	
	Matrix coef() //{{{
	{ return m_coef ; } //}}}
	
	qrstate_t state() //{{{
	{ return m_state ; } //}}}
	
	
	//=======//
	// State //
	//=======//
	
	bool is_fitted()//{{{
	{ return ( m_state == tools::success || m_state == tools::unfeasible ) ; }//}}}
	
	bool is_success()//{{{
	{ return m_state == tools::success ; }//}}}
	
	bool is_unfeasible()//{{{
	{ return m_state == tools::unfeasible ; }//}}}
	
	
	//=========//
	// Methods //
	//=========//
	
	void fit( Array& Y , Matrix& X ) //{{{
	{
		size_type n_tau = m_ltau.size() ;
		size_type n_cov = X.cols() + 1 ;
		m_coef.resize( n_tau , n_cov ) ;
		
		for( size_type i = 0 ; i < n_tau ; ++i )
		{
			m_frishNewtonCore.set_tau( m_ltau[i] ) ;
			m_frishNewtonCore.fit( Y , X , false ) ;
			m_state = m_frishNewtonCore.state() ;
			if( m_state == tools::unfeasible )
				break ;
			m_coef.row(i) = m_frishNewtonCore.m_coef ;
		}
	} //}}}
	
	Matrix predict() //{{{
	{
		size_type n_tau = m_ltau.size() ;
		Matrix& A(m_frishNewtonCore.m_A) ;
		Matrix Yq( A.rows() , n_tau ) ;
		Yq = A * m_coef.transpose() ;
		return Yq ;
	} //}}}
	
	private:
	
	//=========//
	// Methods //
	//=========//
	
	//===========//
	// Arguments //
	//===========//
	
	//{{{
	Array           m_ltau            ;
	Matrix          m_coef            ;
	FrishNewtonCore m_frishNewtonCore ;
	qrstate_t       m_state           ;
	//}}}
	
} ;
//}}}

}

#endif





