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

		// TODO: Comment
		static Eigen::MatrixXd Matlab_fliplr(const Eigen::MatrixXd& x);

		// TODO: Comment
		static Eigen::MatrixXd Matlab_flipud(const Eigen::MatrixXd& x);
	};

}

