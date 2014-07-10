#include <vector>
#include <Eigen/Dense>

//!  DerivativeMatrix class. 
/*!
Provides a wrap up for the return value of chebdif
*/
class DerivativeMatrix
{
public:

	/**
	* Chebyshev points.
	*/
	Eigen::MatrixXd x;
	
	/**
	* Vector of differentiation matrices.
	*/
	std::vector<Eigen::MatrixXd> dm;


};