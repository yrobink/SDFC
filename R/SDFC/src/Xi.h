
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

#ifndef SDFC_QUANTILEREGRESSION_XI_INCLUDED
#define SDFC_QUANTILEREGRESSION_XI_INCLUDED


//===========//
// Libraries //
//===========//


#include <random>
#include <cmath>
#include <limits>

#include <Eigen/Dense>
#include <Eigen/Core>


namespace QuantileRegression
{

//===========//
// Class(es) //
//===========//

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


}

#endif

