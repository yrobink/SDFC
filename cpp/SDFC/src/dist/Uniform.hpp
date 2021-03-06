
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

#ifndef SDFC_UNIFORMLAW
#define SDFC_UNIFORMLAW


//===========//
// Libraries //
//===========//

#include <iostream>
#include <random>
#include <Eigen/Dense>



namespace SDFC
{

// Matrix and vector typedef
//==========================

template <class Scalar> using Matrix = Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> ;
template <class Scalar> using Vector = Eigen::Matrix<Scalar,Eigen::Dynamic,1> ;


// Statistical functions
//======================

template <class Scalar>
Matrix<Scalar> rvs_uniform( size_t n , size_t ncol = 1 , Scalar min = static_cast<Scalar>(0) , Scalar max = static_cast<Scalar>(1) )
{
	Matrix<Scalar> rvs = Matrix<Scalar>::Random(n,ncol) ;
	rvs.array() *= static_cast<Scalar>( (max-min) / 2 ) ;
	rvs.array() += static_cast<Scalar>( (max+min) / 2 ) ;
	return rvs ;
}

//Vector<Scalar> density_uniform( const Vector<Scalar>& , Scalar min = 0 , Scalar max = 1 )





class UniformLaw
{
	public:
	
	// Typedef
	//--------
	//{{{
	typedef double          value_type ;
	typedef unsigned int    size_type  ;
	typedef Eigen::VectorXd Vector     ;
	//}}}
	
	// Constructor / Destructor
	//-------------------------
	
	UniformLaw( value_type min = 0 , value_type max = 1 )://{{{
		m_max(max),
		m_min(min),
		m_gen(),
		m_law(m_min,m_max)
	{}
	//}}}
	
	virtual ~UniformLaw() //{{{
	{}//}}}
	
	// Overloading
	//------------
	
	// Methods
	//--------
	
	Vector rvs( size_type n = 1 ) //{{{
	{
		Eigen::VectorXd out(n) ;
		for( size_type i = 0 ; i < n ; ++i )
			out[i] = m_law(m_gen) ;
		return out ;
	} //}}}
	
	Vector density( Vector x ) //{{{
	{
		return x ;
	}//}}}
	
	Vector cdf( Vector q )//{{{
	{
		Vector p(q.size()) ;
		for( size_type s = 0 ; s < p.size() ; ++s )
			p[s] = q[s] < m_min ? 0. : q[s] > m_max ? 1. : ( q[s] - m_min ) / ( m_max - m_min ) ;
		return p ;
	}//}}}
	
	Vector sf( Vector q )//{{{
	{
		return 1. - UniformLaw::cdf(q).array() ;
	}//}}}
	
	Vector icdf( Vector p ) //{{{
	{
		return (m_max - m_min) * p.array() + m_min ;
	}//}}}
	
	Vector isf( Vector p ) //{{{
	{
		return UniformLaw::icdf( 1. - p.array() ) ;
	}//}}}
	
	
	// Arguments
	//----------
	//{{{
	value_type m_max ;
	value_type m_min ;
	std::default_random_engine     m_gen ;
	std::uniform_real_distribution<double> m_law ;
	//}}}
} ;

} // namespace SDFC

#endif

