#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <string>
#include <utility>

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
		static void eigenToCSV(const Eigen::MatrixXd& x, std::string filename);

		// TODO: Comment
		static Eigen::MatrixXd matlab_fliplr(const Eigen::MatrixXd& x);

		// TODO: Comment
		static Eigen::MatrixXd matlab_flipud(const Eigen::MatrixXd& x);

		// TODO: Comment
		// [V,D] = eig(A,B)
		static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> matlab_eig(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B);

		// TODO: Comment
		// [V,D] = eig(A,B)
		static std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> matlab_ceig(const Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B);
	};

	/*
	template <typename Derived>
	void UtilityMethods<Derived>::eigenToCSV(const Eigen::MatrixBase<Derived>& x, std::string filename)
	{
		const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
		std::ofstream of(filename);
		if (!of) {
			std::cout << "Cannot open file.\n";
		}
		std::cout.setf(ios::fixed);
		of << setprecision(std::numeric_limits<double>::digits10 + 2) << x.format(CSVFormat) << endl;
		of.close();
	}
	template <typename Derived>
	typename Derived UtilityMethods<Derived>::matlab_fliplr(Eigen::MatrixBase<Derived> const& x)
	{
		return x.rowwise().reverse().eval();
	}
	
	template<typename Derived>
	typename Derived UtilityMethods<Derived>::matlab_flipud(Eigen::MatrixBase<Derived> const& x)
	{
		return x.colwise().reverse().eval();
	}
	
	// [V,D] = eig(A,B)
	template <typename Derived>
	std::pair<typename Derived, typename Derived> UtilityMethods<Derived>::matlab_eig(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B)
	{
		const Eigen::GeneralizedSelfAdjointEigenSolver<Derived> solver(A, B);
	
		const Derived v = solver.eigenvectors().cast<Derived>();
		const d = solver.eigenvalues().cast<Derived>();
	
		return std::make_pair(v, d);
	
	}

	*/
}

#endif