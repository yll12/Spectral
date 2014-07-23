/*
#include "utils.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace Eigen;
namespace Utility
{
*/

	/*!
	*  Write Eigen Matrices to files
	*/
	/*
	template <const typename Derived>
	void UtilityMethods::eigenToCSV(const MatrixBase<Derived>& x, string filename)
	{
		const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
		ofstream of(filename);
		if (!of) {
			cout << "Cannot open file.\n";
		}
		cout.setf(ios::fixed);
		of << setprecision(numeric_limits<double>::digits10 + 2) << x.format(CSVFormat) << endl;
		of.close();
	}
	
	MatrixXd UtilityMethods::matlab_fliplr(const MatrixXd& x)
	{
		return x.rowwise().reverse().eval();
	}

	MatrixXd UtilityMethods::matlab_flipud(const MatrixXd& x)
	{
		return x.colwise().reverse().eval();
	}
	
	// [V,D] = eig(A,B)
	pair<MatrixXd, MatrixXd> UtilityMethods::matlab_eig(const MatrixXd& A, const MatrixXd& B)
	{
		const Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(A, B);

		const MatrixXd v = solver.eigenvectors();
		const MatrixXd d = solver.eigenvalues();

		return std::make_pair(v, d);

	}

	// [V,D] = eig(A,B)
	pair<MatrixXcd, MatrixXcd> UtilityMethods::matlab_ceig(const MatrixXcd& A, const MatrixXcd& B)
	{
		const Eigen::GeneralizedSelfAdjointEigenSolver<MatrixXcd> solver(A, B);

		const MatrixXcd v = solver.eigenvectors().cast<complex<double>>();
		const MatrixXcd d = solver.eigenvalues().cast<complex<double>>();

		return std::make_pair(v, d);

	}

}
	*/

