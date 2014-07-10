#ifndef UTILITY_H_
#define UTILITY_H_

#include <Eigen/Dense>
#include <string>

namespace Utility
{
	//!  A utils class. 
	/*!
	 Provides various general purpose methods.
	*/
	class UtilityMethods
	{
	public:
		//! A static method taking an Eigen library matrices and write it to specified file.
		/*!
		\param x Eigen Matrix to be written.
		\param filename Filename to be written to.
		*/
		static void EigenToCSV(Eigen::MatrixXd& x, std::string filename);
	};

}

#endif