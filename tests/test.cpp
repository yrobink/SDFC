
#include <iostream>
#include <random>
#include <Eigen/Dense>

#include "../cpp/QuantileRegression.hpp"

int main()
{
	// Dataset
	std::size_t size(2000) ;
	Eigen::VectorXd t( Eigen::VectorXd::LinSpaced( size , 0 , 1 ) ) , Y(size) ;
	Eigen::MatrixXd X( size , 1 ) ;
	X.col(0) = t.array().pow(2) ;
	
	std::default_random_engine generator;
	std::normal_distribution<double> gauss(0.,0.1);
	for( std::size_t s = 0 ; s < size ; ++s )
		Y[s] = X(s,0) + gauss(generator) ;
	
	Eigen::VectorXd ltau( Eigen::VectorXd::LinSpaced( 100 , 0.1 , 0.9 ) ) ;
	
	// Fit
	SDFC::QuantileRegression qreg( ltau ) ;
	qreg.fit( Y , X ) ;
	Eigen::MatrixXd quantiles = qreg.predict() ;
	
	std::cout << qreg.repr() << std::endl ;
	
	return 0 ;
}
