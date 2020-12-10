// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

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


#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
#include "FrishNewton.h"


RCPP_MODULE(SDFC_cpp){
	// Classes
	
	Rcpp::class_<QuantileRegression::FrishNewton>("FrishNewton")
	.constructor<Eigen::ArrayXd,int,double,double>()
	.method( "set_ltau" , &QuantileRegression::FrishNewton::set_ltau , "" )
	.method( "coef"     , &QuantileRegression::FrishNewton::coef     , "" )
	.method( "is_fitted"     , &QuantileRegression::FrishNewton::is_fitted     , "" )
	.method( "is_success"    , &QuantileRegression::FrishNewton::is_success    , "" )
	.method( "is_unfeasible" , &QuantileRegression::FrishNewton::is_unfeasible , "" )
	.method( "fit"     , &QuantileRegression::FrishNewton::fit     , "" )
	.method( "predict" , &QuantileRegression::FrishNewton::predict , "" )

	;
}



