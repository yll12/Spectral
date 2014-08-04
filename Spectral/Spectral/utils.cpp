/*
#include "utils.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace Eigen;
namespace Utility
{

	std::pair<MatrixXcd, MatrixXcd> UtilityMethods::matlab_ceig(const MatrixXcd& A, const MatrixXcd& B)
	{
		MatrixXcd v, lambda;

		int N = A.cols(); // Number of columns of A and B. Number of rows of v.

		char jobvl = 'N', jobvr = 'V';
		int n, lda, ldb, ldvl, ldvr, lwork, info;
		n = lda = A.rows();
		ldb = B.rows();
		ldvl = 1;
		ldvr = n;
		lwork = std::max(1, n*n + 64); // This may be choosen better!
		MatrixXcd work(lwork, 1);
		MatrixXd rwork(8 * n, 1); // This may be choosen better
		MatrixXcd alpha(n, 1), beta(n, 1);
		MatrixXcd vl(1, 1), vr(n, n);
		zggev_(&jobvl, &jobvr, &n, A.data(), &lda,
			B.data(), &ldb, alpha.data(), beta.data(), vl.data(),
			&ldvl, vr.data(), &ldvr, work.data(), &lwork,
			rwork.data(), &info);
		lambda = alpha.cwiseQuotient(beta);
		v = vr;

		return std::make_pair(v, lambda);
	}

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
