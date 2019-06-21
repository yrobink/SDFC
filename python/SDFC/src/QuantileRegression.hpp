
//==============================================================================//
//==============================================================================//
//                                                                              //
// Copyright Yoann Robin, 2019                                                  //
//                                                                              //
// yoann.robin.k@gmail.com                                                      //
//                                                                              //
// This software is a computer program that is part of the SDFC (Statistical    //
// Distribution Fit with Covariates) library. This library makes it possible    //
// to regress the parameters of some statistical law with co-variates.          //
//                                                                              //
// This software is governed by the CeCILL-C license under French law and       //
// abiding by the rules of distribution of free software.  You can  use,        //
// modify and/ or redistribute the software under the terms of the CeCILL-C     //
// license as circulated by CEA, CNRS and INRIA at the following URL            //
// "http://www.cecill.info".                                                    //
//                                                                              //
// As a counterpart to the access to the source code and  rights to copy,       //
// modify and redistribute granted by the license, users are provided only      //
// with a limited warranty  and the software's author,  the holder of the       //
// economic rights,  and the successive licensors  have only  limited           //
// liability.                                                                   //
//                                                                              //
// In this respect, the user's attention is drawn to the risks associated       //
// with loading,  using,  modifying and/or developing or reproducing the        //
// software by the user in light of its specific status of free software,       //
// that may mean  that it is complicated to manipulate,  and  that  also        //
// therefore means  that it is reserved for developers  and  experienced        //
// professionals having in-depth computer knowledge. Users are therefore        //
// encouraged to load and test the software's suitability as regards their      //
// requirements in conditions enabling the security of their systems and/or     //
// data to be ensured and,  more generally, to use and operate it in the        //
// same conditions as regards security.                                         //
//                                                                              //
// The fact that you are presently reading this means that you have had         //
// knowledge of the CeCILL-C license and that you accept its terms.             //
//                                                                              //
//==============================================================================//
//==============================================================================//

//==============================================================================//
//==============================================================================//
//                                                                              //
// Copyright Yoann Robin, 2019                                                  //
//                                                                              //
// yoann.robin.k@gmail.com                                                      //
//                                                                              //
// Ce logiciel est un programme informatique faisant partie de la librairie     //
// SDFC (Statistical Distribution Fit with Covariates). Cette librairie         //
// permet de calculer de regresser les parametres de lois statistiques selon    //
// plusieurs co-variables                                                       //
//                                                                              //
// Ce logiciel est régi par la licence CeCILL-C soumise au droit français et    //
// respectant les principes de diffusion des logiciels libres. Vous pouvez      //
// utiliser, modifier et/ou redistribuer ce programme sous les conditions       //
// de la licence CeCILL-C telle que diffusée par le CEA, le CNRS et l'INRIA     //
// sur le site "http://www.cecill.info".                                        //
//                                                                              //
// En contrepartie de l'accessibilité au code source et des droits de copie,    //
// de modification et de redistribution accordés par cette licence, il n'est    //
// offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,    //
// seule une responsabilité restreinte pèse sur l'auteur du programme, le       //
// titulaire des droits patrimoniaux et les concédants successifs.              //
//                                                                              //
// A cet égard  l'attention de l'utilisateur est attirée sur les risques        //
// associés au chargement,  à l'utilisation,  à la modification et/ou au        //
// développement et à la reproduction du logiciel par l'utilisateur étant       //
// donné sa spécificité de logiciel libre, qui peut le rendre complexe à        //
// manipuler et qui le réserve donc à des développeurs et des professionnels    //
// avertis possédant  des  connaissances  informatiques approfondies.  Les      //
// utilisateurs sont donc invités à charger  et  tester  l'adéquation  du       //
// logiciel à leurs besoins dans des conditions permettant d'assurer la         //
// sécurité de leurs systèmes et ou de leurs données et, plus généralement,     //
// à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.           //
//                                                                              //
// Le fait que vous puissiez accéder à cet en-tête signifie que vous avez       //
// pris connaissance de la licence CeCILL-C, et que vous en avez accepté les    //
// termes.                                                                      //
//                                                                              //
//==============================================================================//
//==============================================================================//

#ifndef SDFC_NONPARAMETRIC_QUANTILEREGRESSION
#define SDFC_NONPARAMETRIC_QUANTILEREGRESSION

//-----------//
// Libraries //
//-----------//


#include <limits>
#include <cmath>
#include <random>
#include <pybind11/pybind11.h>
#include <Eigen/Dense>
#include <Eigen/Core>


#include "FrishNewton.hpp"

//============//
// namespaces //
//============//

namespace py = pybind11 ;

//=======//
// Class //
//=======//

struct QuantileRegression
{
	//=========//
	// Typedef // 
	//=========//
	
	//{{{
	typedef unsigned int    size_type  ;
	typedef double          value_type ;
	typedef Eigen::ArrayXd  Array      ;
	typedef Eigen::VectorXd Vector     ;
	typedef Eigen::MatrixXd Matrix     ;
	//}}}
	
	
	//=============//
	// Constructor //
	//=============//
	
	QuantileRegression(): //{{{
		m_ltau() ,
		m_coef() ,
		m_quantiles() ,
		m_frishNewton( 0.5 , 50 , 1e-6 , 0.99995 ),
		m_state(not_fitted)
	{} //}}}
	
	QuantileRegression( value_type ltau ): //{{{
		m_ltau(1) ,
		m_coef() ,
		m_quantiles() ,
		m_frishNewton( 0.5 , 50 , 1e-6 , 0.99995 ),
		m_state(not_fitted)
	{
		m_ltau[0] = ltau ;
	} //}}}
	
	QuantileRegression( py::list ltau ): //{{{
		m_ltau(py::len(ltau)) ,
		m_coef() ,
		m_quantiles() ,
		m_frishNewton( 0.5 , 50 , 1e-6 , 0.99995 ),
		m_state(not_fitted)
	{
		int s = 0 ;
		for( auto item : ltau )
			m_ltau[s++] = item.cast<value_type>() ;
	} //}}}
	
	QuantileRegression( Eigen::Ref<Vector> ltau ): //{{{
		m_ltau(ltau) ,
		m_coef() ,
		m_quantiles() ,
		m_frishNewton( 0.5 , 50 , 1e-6 , 0.99995 ),
		m_state(not_fitted)
	{} //}}}
	
	~QuantileRegression() //{{{
	{} //}}}
	
	std::string repr()//{{{
	{
		std::string True("True"), False("False") ;
		std::string _repr("") ;
		_repr += "SDFC.NonParametric.QuantileRegression\n" ;
		_repr += "=====================================\n" ;
		_repr += "* Method    : Frish-Newton\n" ;
		_repr += "* Fitted    : " + ( is_fitted()     ? True : False ) + "\n" ;
		_repr += "* Success   : " + ( is_success()    ? True : False ) + "\n" ;
		_repr += "* Unfeasible: " + ( is_unfeasible() ? True : False ) + "\n" ;
		return _repr ;
	}
	//}}}
	
	
	//===========//
	// Accessors //
	//===========//
	
	void set_fit_params( size_type maxit , value_type tol , value_type beta ) //{{{
	{
		m_frishNewton.m_maxit = maxit ;
		m_frishNewton.m_tol   = tol ;
		m_frishNewton.m_beta  = beta ;
	}
	//}}}
	
	void set_ltau( double ltau ) //{{{
	{
		m_ltau.resize(1) ;
		m_ltau[0] = ltau ;
	} //}}}
	
	void set_ltau( py::list ltau ) //{{{
	{
		m_ltau.resize( py::len(ltau) ) ;
		int s = 0 ;
		for( auto item : ltau )
			m_ltau[s++] = item.cast<value_type>() ;
	} //}}}
	
	void set_ltau( Eigen::Ref<Vector> ltau ) //{{{
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
	{ return ( m_state == success || m_state == unfeasible ) ; }//}}}
	
	bool is_success()//{{{
	{ return m_state == success ; }//}}}
	
	bool is_unfeasible()//{{{
	{ return m_state == unfeasible ; }//}}}
	
	
	//=========//
	// Methods //
	//=========//
	
	void fit( Vector& Y , Matrix& X ) //{{{
	{
		size_type n_tau = m_ltau.size() ;
		size_type n_cov = X.cols() + 1 ;
		m_coef.resize( n_tau , n_cov ) ;
		for( size_type i = 0 ; i < n_tau ; ++i )
		{
			m_frishNewton.set_tau( m_ltau[i] ) ;
			m_frishNewton.fit( Y , X ) ;
			m_state = m_frishNewton.state() ;
			if( m_state == unfeasible )
				break ;
			m_coef.row(i) = m_frishNewton.m_coef ;
		}
		m_quantiles = predict() ;
	} //}}}
	
	Matrix predict() //{{{
	{
		size_type n_tau = m_ltau.size() ;
		Matrix& A(m_frishNewton.m_A) ;
		Matrix Yq( A.rows() , n_tau ) ;
		Yq = A * m_coef.transpose() ;
		return Yq ;
	} //}}}
	
	
	//===========//
	// Arguments //
	//===========//
	
	//{{{
	Array       m_ltau        ;
	Matrix      m_coef        ;
	Matrix      m_quantiles   ;
	FrishNewton m_frishNewton ;
	qrstate_t   m_state       ;
	//}}}
	
	
	
};

#endif

