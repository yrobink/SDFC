
// Copyright(c) 2020 Yoann Robin
// 
// This file is part of SDFC.
// 
// SDFC is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SDFC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with SDFC.  If not, see <https://www.gnu.org/licenses/>.

//==============================================================================//
// This software is based on the fortran code developped by Roger Koenker,      //
// coming from the R package "quantreg", see:                                   //
// https://cran.r-project.org/web/packages/quantreg/index.html                  //
//==============================================================================//

//==============================================================================//
// Ce programme est basé sur le code fortran développé par Roger Koenker,       //
// venant du package R "quantreg", voir:                                        //
// https://cran.r-project.org/web/packages/quantreg/index.html                  //
//==============================================================================//


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





