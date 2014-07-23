#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <string>
#include <utility>
#include <fstream>
#include <iostream>
#include <iomanip>

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
		template <typename Derived>
		static void eigenToCSV(const Eigen::MatrixBase<Derived>& x, std::string filename);

		// TODO: Comment
		template <typename Derived>
		static Derived matlab_fliplr(const Eigen::MatrixBase<Derived>& x);

		// TODO: Comment
		template <typename Derived>
		static Derived matlab_flipud(const Eigen::MatrixBase<Derived>& x);

		// TODO: Comment
		// [V,D] = eig(A,B)
		template <typename Derived>
		static std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> matlab_ceig(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B);

		// TODO: Comment
		// [V,D] = eig(A,B)
		template <typename Derived>
		static std::pair<Eigen::MatrixXd, Eigen::MatrixXd> matlab_eig(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B);

	};

	template <typename Derived>
	void UtilityMethods::eigenToCSV(const Eigen::MatrixBase<Derived>& x, std::string filename)
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
	typename Derived UtilityMethods::matlab_fliplr(const Eigen::MatrixBase<Derived>& x)
	{
		return x.rowwise().reverse().eval();
	}
	
	template<typename Derived>
	typename Derived UtilityMethods::matlab_flipud(const Eigen::MatrixBase<Derived>& x)
	{
		return x.colwise().reverse().eval();
	}
	
	// [V,D] = eig(A,B)
	template <typename Derived>
	std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> UtilityMethods::matlab_ceig(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B)
	{
		const Eigen::GeneralizedSelfAdjointEigenSolver<Derived> solver(A, B);
	
		const Eigen::MatrixXcd v = solver.eigenvectors().cast<complex<double>>();
		const Eigen::MatrixXcd d = solver.eigenvalues().cast<complex<double>>();
	
		return std::make_pair(v, d);
	
	}

	// [V,D] = eig(A,B)
	template <typename Derived>
	std::pair<Eigen::MatrixXd, Eigen::MatrixXd> UtilityMethods::matlab_eig(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B)
	{
		const Eigen::GeneralizedSelfAdjointEigenSolver<Derived> solver(A, B);

		const Eigen::MatrixXd v = solver.eigenvectors();
		const Eigen::MatrixXd d = solver.eigenvalues();

		return std::make_pair(v, d);

	}

}

#endif