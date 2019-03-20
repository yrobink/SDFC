
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

#ifndef SDFC_QUANTILEREGRESSION_PYFRISHNEWTON_INCLUDED
#define SDFC_QUANTILEREGRESSION_PYFRISHNEWTON_INCLUDED

#include <vector>
#include "src/FrishNewton.hpp"


class PyFrishNewton
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
	//}}}
	
	
	//=============//
	// Constructor //
	//=============//
	
	PyFrishNewton()://{{{
		m_frishNewton(),
		m_state(QuantileRegression::tools::not_fitted)
	{}
	//}}}
	
	PyFrishNewton( std::vector<value_type>& ltau , int maxit = 50 , value_type tol = 1e-6 , value_type beta = 0.99995 )://{{{
		m_frishNewton( static_cast<unsigned int>(maxit) , tol , beta ),
		m_state(QuantileRegression::tools::not_fitted)
	{
		size_type n_tau(ltau.size()) ;
		Array eltau( n_tau ) ;
		
		for( size_type i = 0 ; i < n_tau ; ++i )
		{
			eltau[i] = ltau[i] ;
		}
		m_frishNewton.set_ltau(eltau) ;
	}
	//}}}
	
	~PyFrishNewton()//{{{
	{}
	//}}}
	
	
	//===========//
	// Accessors //
	//===========//
	
	std::vector<std::vector<value_type>> coef() //{{{
	{
		Matrix ecoef = m_frishNewton.coef() ;
		size_type nrow(ecoef.rows()) , ncol(ecoef.cols()) ;
		std::vector<std::vector<value_type>> mcoef( nrow ) ;
		
		for( size_type i = 0 ; i < nrow ; ++i )
		{
			mcoef[i].resize(ncol) ;
			for( size_type j = 0 ; j < ncol ; ++j )
			{
				mcoef[i][j] = ecoef(i,j) ;
			}
		}
		return mcoef ;
	}
	//}}}
	
	int state() //{{{
	{
		switch(m_state)
		{
			case QuantileRegression::tools::success:
				return 0 ;
			case QuantileRegression::tools::unfeasible:
				return 1 ;
			case QuantileRegression::tools::not_fitted:
				return 2 ;
			default:
				return 2 ;
		}
	}
	//}}}
	
	
	//=========//
	// Methods //
	//=========//
	
	void fit( std::vector<value_type>& Y , std::vector<std::vector<value_type>>& X ) //{{{
	{
		Array eY(Y.size()) ;
		Matrix eX( Y.size() , X[0].size() ) ;
		for( size_type i = 0 ; i < Y.size() ; ++i )
		{
			eY[i] = Y[i] ;
			for( size_type j = 0 ; j < X[0].size() ; ++j )
			{
				eX(i,j) = X[i][j] ;
			}
		}
		m_frishNewton.fit( eY , eX ) ;
		m_state = m_frishNewton.state() ;
	}
	//}}}
	
	std::vector<std::vector<value_type>> predict() //{{{
	{
		Matrix eYq = m_frishNewton.predict() ;
		size_type nrow(eYq.rows()) , ncol(eYq.cols()) ;
		std::vector<std::vector<value_type>> Yq( eYq.rows() ) ;
		
		for( size_type i = 0 ; i < nrow ; ++i )
		{
			Yq[i].resize(ncol) ;
			for( size_type j = 0 ; j < ncol ; ++j )
			{
				Yq[i][j] = eYq(i,j) ;
			}
		}
		return Yq ;
	}
	//}}}
	
	
	private:
	
	//===========//
	// Arguments //
	//===========//
	
	//{{{
	QuantileRegression::FrishNewton      m_frishNewton ;
	QuantileRegression::tools::qrstate_t m_state       ;
	//}}}
} ;


#endif
