
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


#ifndef SDFC_NONPARAMETRIC_FRISHNEWTON
#define SDFC_NONPARAMETRIC_FRISHNEWTON

//-----------//
// Libraries //
//-----------//


#include <limits>
#include <cmath>
#include <random>

#include <Eigen/Dense>
#include <Eigen/Core>


//=======//
// Class //
//=======//

namespace SDFC {

enum qrstate_t //{{{
{
	success = 0,
	unfeasible = 1,
	not_fitted = 2
} ;
//}}}

struct Xi //{{{
{
	typedef unsigned int size_type  ;
	typedef double       value_type ;
	
	Xi():
		y(),
		z(),
		x(),
		s(),
		w()
	{}
	
	Xi( size_type p , size_type n ):
		y(p),
		z(n),
		x(n),
		s(n),
		w(n)
	{}
	
	~Xi()
	{}
	
	
	value_type mu()
	{
		return value_type(x.matrix().transpose() * z.matrix()) + value_type(s.matrix().transpose() * w.matrix()) ;
	}
	
	Eigen::ArrayXd y ;
	Eigen::ArrayXd z ;
	Eigen::ArrayXd x ;
	Eigen::ArrayXd s ;
	Eigen::ArrayXd w ;
} ;
// }}}

struct QuantileRegression ;

class FrishNewton
{
	public:
	friend QuantileRegression ;
	
	//=========//
	// Typedef // 
	//=========//
	
	//{{{
	typedef unsigned int    size_type          ;
	typedef double          value_type         ;
	typedef Eigen::ArrayXd  Array              ;
	typedef Eigen::VectorXd Vector             ;
	typedef Eigen::MatrixXd Matrix             ;
	typedef Eigen::FullPivLU<Matrix> PLUMatrix ;
	//}}}
	
	//=============//
	// Constructor //
	//=============//
	
	FrishNewton( value_type tau , size_type maxit = 50 , value_type tol = 1e-6 , value_type beta = 0.99995 ): //{{{
		m_tau(tau) ,
		m_maxit(maxit),
		m_tol(tol),
		m_beta(beta),
		m_p(),
		m_n(),
		m_A(),
		m_lu(),
		m_iAQA(),
		m_xi(),
		m_dxi(),
		m_alphaP(),
		m_alphaD(),
		m_mu(0),
		m_b(),
		m_c(),
		m_iq(),
		m_u(),
		m_dr(),
		m_rhs(),
		m_coef(),
		m_state(not_fitted)
	{} //}}}
	
	~FrishNewton() //{{{
	{} //}}}
	
	//===========//
	// Accessors //
	//===========//
	
	void set_tau( value_type tau ) //{{{
	{
		m_tau = tau ;
	}//}}}
	
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
		
		initialize( Y , X ) ;
		if( m_state == unfeasible )
			return ;
		
		size_type nit = 0 ;
		while( m_mu > m_tol && ++nit < m_maxit )
		{
			predictor_step() ;
			if( m_state == unfeasible )
				return ;
			infer_gap() ;
			corrector_step() ;
			update() ;
		}
		m_state = success ;
		m_coef = - m_xi.y ;
		
	} //}}}
	
	Array predict() //{{{
	{
		return m_A * m_coef.matrix() ;
	} //}}}
	
	
	protected:
	
	//=========//
	// Methods //
	//=========//
	
	void initialize( Vector& Y , Matrix& X ) //{{{
	{
		m_n = Y.size() ;
		m_p = X.cols() + 1 ;
		m_xi  = Xi(m_p,m_n) ;
		m_dxi = Xi(m_p,m_n) ;
		m_c   = -Y ;
		m_iq  = Array(m_n) ;
		m_u   = Array(m_n) ;
		m_dr  = Array(m_n) ;
		m_rhs = Array(m_p) ;
		
		// Design matrix
		m_A = Matrix(m_n,m_p) ;
		m_A.col(0).array() = 1. ;
		m_A.block( 0 , 1 , m_n , m_p - 1 ) = X ;
		
		//
		m_b = ( 1. - m_tau ) * m_A.colwise().sum() ;
		m_lu = PLUMatrix(m_A.transpose() * m_A) ;
		if( !m_lu.isInvertible() )
		{
			m_state = unfeasible ;
			return ;
		}
		m_iAQA = m_lu.inverse() ;
		m_iq = 1. ;
		m_u  = 1. ;
		m_xi.x = 1. - m_tau ;
		m_xi.y = m_A.transpose() * m_c.matrix() ;
		m_xi.y = m_iAQA * m_xi.y.matrix() ;
		m_xi.s = m_c - (m_A * m_xi.y.matrix()).array() ;
		
		for( size_type i = 0 ; i < m_n ; ++i )
		{
			m_xi.z[i] = std::max(  m_xi.s[i] , 0. ) + ( (std::abs(m_xi.s[i]) < m_tol) ? m_tol : 0 ) ;
			m_xi.w[i] = std::max( -m_xi.s[i] , 0. ) + ( (std::abs(m_xi.s[i]) < m_tol) ? m_tol : 0 ) ;
		}
		m_xi.s = m_u - m_xi.x ;
		m_mu = m_xi.mu() ;
	} //}}}
	
	void predictor_step() //{{{
	{
		for( size_type i = 0 ; i < m_n ; ++i )
		{
			m_iq[i] = 1. / ( m_xi.z[i] / m_xi.x[i] + m_xi.w[i] / m_xi.s[i] ) ;
			m_dxi.s[i] = m_xi.z[i] - m_xi.w[i] ;
			m_dxi.z[i] = m_iq[i] * m_dxi.s[i] ;
		}
		
		m_lu = PLUMatrix(m_A.transpose() * m_iq.matrix().asDiagonal() * m_A) ;
		if( !m_lu.isInvertible() )
		{
			m_state = unfeasible ;
			return ;
		}
		m_iAQA = m_lu.inverse() ;
		m_rhs = m_b - (m_A.transpose() * m_xi.x.matrix()).array() ;
		m_rhs += (m_A.transpose() * m_dxi.z.matrix()).array() ;
		m_dxi.y = m_iAQA * m_rhs.matrix() ;
		m_dxi.s = ( m_A * m_dxi.y.matrix() ).array() - m_dxi.s ;
		
		m_alphaP = std::numeric_limits<value_type>::max() ;
		m_alphaD = std::numeric_limits<value_type>::max() ;
		
		for( size_type i = 0 ; i < m_n ; ++i )
		{
			m_dxi.x[i] = m_iq[i] * m_dxi.s[i] ;
			m_dxi.s[i] = - m_dxi.x[i] ;
			m_dxi.z[i] = - m_xi.z[i] * ( m_dxi.x[i] / m_xi.x[i] + 1. ) ;
			m_dxi.w[i] = - m_xi.w[i] * ( m_dxi.s[i] / m_xi.s[i] + 1. ) ;
			
			if( m_dxi.x[i] < 0 )
				m_alphaP = std::min( m_alphaP , - m_xi.x[i] / m_dxi.x[i] ) ;
			if( m_dxi.s[i] < 0 )
				m_alphaP = std::min( m_alphaP , - m_xi.s[i] / m_dxi.s[i] ) ;
			if( m_dxi.z[i] < 0 )
				m_alphaD = std::min( m_alphaD , - m_xi.z[i] / m_dxi.z[i] ) ;
			if( m_dxi.w[i] < 0 )
				m_alphaD = std::min( m_alphaD , - m_xi.w[i] / m_dxi.w[i] ) ;
		}
		m_alphaP = std::min( m_beta * m_alphaP , 1. ) ;
		m_alphaD = std::min( m_beta * m_alphaD , 1. ) ;
	}//}}}
	
	void infer_gap() //{{{
	{
		m_mu = m_xi.mu() ;
		value_type mu_new = 0 ;
		for( size_type i = 0 ; i < m_n ; ++i )
		{
			mu_new += (m_xi.x[i] + m_alphaP * m_dxi.x[i]) * (m_xi.z[i] + m_alphaD * m_dxi.z[i]) + (m_xi.s[i] + m_alphaP * m_dxi.s[i]) * (m_xi.w[i] + m_alphaD * m_dxi.w[i]) ;
		}
		m_mu = std::pow( mu_new , 3 ) / std::pow( m_mu , 2 ) / ( 2. * static_cast<value_type>(m_n) ) ;
	}//}}}
	
	void corrector_step() //{{{
	{
		for( size_type i = 0 ; i < m_n ; ++i )
		{
			m_dr[i] = m_iq[i] * ( m_mu * ( 1. / m_xi.s[i] - 1. / m_xi.x[i] ) + m_dxi.x[i] * m_dxi.z[i] / m_xi.x[i] - m_dxi.s[i] * m_dxi.w[i] / m_xi.s[i] ) ;
		}
		std::swap( m_dxi.y , m_rhs ) ;
		m_dxi.y += ( m_A.transpose() * m_dr.matrix() ).array() ;
		m_dxi.y = m_iAQA * m_dxi.y.matrix() ;
		m_u = m_A * m_dxi.y.matrix() ;
		
		m_alphaP = std::numeric_limits<value_type>::max() ;
		m_alphaD = std::numeric_limits<value_type>::max() ;
		
		value_type dxdz, dsdw ;
		for( size_type i = 0 ; i < m_n ; ++i )
		{
			dxdz = m_dxi.x[i] * m_dxi.z[i] ;
			dsdw = m_dxi.s[i] * m_dxi.w[i] ;
			m_dxi.x[i] = m_iq[i] * ( m_u[i] - m_xi.z[i] + m_xi.w[i] ) - m_dr[i] ;
			m_dxi.s[i] = - m_dxi.x[i] ;
			m_dxi.z[i] = - m_xi.z[i] + ( m_mu - m_xi.z[i] * m_dxi.x[i] - dxdz ) / m_xi.x[i] ;
			m_dxi.w[i] = - m_xi.w[i] + ( m_mu - m_xi.w[i] * m_dxi.s[i] - dsdw ) / m_xi.s[i] ;
			
			if( m_dxi.x[i] < 0 )
				m_alphaP = std::min( m_alphaP , - m_xi.x[i] / m_dxi.x[i] ) ;
			if( m_dxi.s[i] < 0 )
				m_alphaP = std::min( m_alphaP , - m_xi.s[i] / m_dxi.s[i] ) ;
			if( m_dxi.z[i] < 0 )
				m_alphaD = std::min( m_alphaD , - m_xi.z[i] / m_dxi.z[i] ) ;
			if( m_dxi.w[i] < 0 )
				m_alphaD = std::min( m_alphaD , - m_xi.w[i] / m_dxi.w[i] ) ;
		}
		m_alphaP = std::min( m_beta * m_alphaP , 1. ) ;
		m_alphaD = std::min( m_beta * m_alphaD , 1. ) ;
	}//}}}
	
	void update() //{{{
	{
		m_xi.x += m_alphaP * m_dxi.x ;
		m_xi.s += m_alphaP * m_dxi.s ;
		m_xi.y += m_alphaD * m_dxi.y ;
		m_xi.z += m_alphaD * m_dxi.z ;
		m_xi.w += m_alphaD * m_dxi.w ;
		m_mu = m_xi.mu() ;
	}//}}}
	
	
	//===========//
	// Arguments //
	//===========//
	
	//{{{
	value_type m_tau              ;
	size_type  m_maxit            ;
	value_type m_tol              ;
	value_type m_beta             ;
	size_type  m_p                ;
	size_type  m_n                ;
	Matrix     m_A                ;
	PLUMatrix  m_lu               ;
	Matrix     m_iAQA             ;
	Xi         m_xi               ;
	Xi         m_dxi              ;
	value_type m_alphaP           ;
	value_type m_alphaD           ;
	value_type m_mu               ;
	Array      m_b                ;
	Array      m_c                ;
	Array      m_iq               ;
	Array      m_u                ;
	Array      m_dr               ;
	Array      m_rhs              ;
	Array      m_coef             ;
	qrstate_t  m_state            ;
	//}}}
	

} ;

}

#endif

