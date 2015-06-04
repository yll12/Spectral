#ifndef UTILS_H_
#define UTILS_H_

#include <Eigen/Dense>
#include <string>
#include <utility>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>

extern "C" void dggev_(const char* JOBVL, const char* JOBVR, const int* N,
	const double* A, const int* LDA, const double* B, const int* LDB,
	double* ALPHAR, double* ALPHAI, double* BETA,
	double* VL, const int* LDVL, double* VR, const int* LDVR,
	double* WORK, const int* LWORK, int* INFO);

extern "C" void zggev_(const char *jobvl, const char *jobvr, int *n, const std::complex<double> *a,
		int *lda, const std::complex<double> *b, int *ldb, std::complex<double> *alpha,
		std::complex<double> *beta, std::complex<double> *vl,
		int *ldvl, std::complex<double> *vr, int *ldvr,
		std::complex<double> *work, int *lwork, double *rwork, int *info);


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
		static std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> matlab_eig(const Derived& A, const Derived& B);

		// TODO: Comment
		static bool compare(const double& x, const double& y);

		//! A static method that computes the differentiation
		//!matrices D1, D2, ..., DM on Chebyshev nodes. \n
		//!The code implements two strategies for enhanced
		//!accuracy suggested by W.Don and S.Solomonoff in
		//!SIAM J.Sci.Comp.Vol. 6, pp. 1253--1268 (1994).
		//!The two strategies are(a) the use of trigonometric
		//!identities to avoid the computation of differences
		//!x(k) - x(j) and(b) the use of the "flipping trick"
		//!which is necessary since sin t can be computed to high
		//!relative precision when t is small whereas sin(pi - t) cannot. \n
		//!
		//!J.A.C.Weideman, S.C.Reddy 1998.
		/*!
		\param n Size of differentiation matrix.
		\param m Number of derivatives required (integer). \n
		Note: 0 < m <= n-1.
		\return pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> x
				x.first - Chebyshev points.
				x.second - Vector of differentiation matrices.
		*/
		static std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> chebdif(int n, int m);

		static Eigen::MatrixXd rotate_YToZ(Eigen::MatrixXd c);

		static Eigen::MatrixXd rotate_YToX(Eigen::MatrixXd c);

		static std::vector<Eigen::MatrixXcd> compute_EOM(Eigen::MatrixXcd D, 
			Eigen::MatrixXcd N, double h, double k);

		static Eigen::MatrixXd normalise(Eigen::MatrixXd c);
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
	
	template <typename Derived>
	typename Derived UtilityMethods::matlab_flipud(const Eigen::MatrixBase<Derived>& x)
	{
		return x.colwise().reverse().eval();
	}

	// [V,D] = eig(A,B)
	template <typename Derived>
	std::pair<Eigen::MatrixXcd, Eigen::MatrixXcd> UtilityMethods::matlab_eig(const Derived& A, const Derived& B)
	{
		Eigen::MatrixXcd v, lambda;

		Eigen::MatrixXcd temp_A = A.cast<complex<double>>();
		Eigen::MatrixXcd temp_B = B.cast<complex<double>>();

		char jobvl = 'N', jobvr = 'V';
		int n, lda, ldb, ldvl, ldvr, lwork, info;
		n = lda = A.rows();
		ldb = B.rows();
		ldvl = 1;
		ldvr = n;
		lwork = std::max(1, n*n + 64); // This may be choosen better!
		Eigen::MatrixXcd work(lwork, 1);
		Eigen::MatrixXd rwork(8 * n, 1); // This may be choosen better
		Eigen::MatrixXcd alpha(n, 1), beta(n, 1);
		Eigen::MatrixXcd vl(1, 1), vr(n, n);
		zggev_(&jobvl, &jobvr, &n, temp_A.data(), &lda,
			temp_B.data(), &ldb, alpha.data(), beta.data(), vl.data(),
			&ldvl, vr.data(), &ldvr, work.data(), &lwork,
			rwork.data(), &info);
		lambda = alpha.cwiseQuotient(beta);
		v = vr;

		return std::make_pair(v, lambda);
	}

	// OLD Impementation:
	/*
	// [V,D] = eig(A,B)
	template <typename Derived>
	std::pair<typename Derived, typename Derived> UtilityMethods::matlab_eig(const Derived& A, const Derived& B)
	{
	Derived v, lambda;

	int N = A.cols(); // Number of columns of A and B. Number of rows of v.

	v.resize(N, N);
	lambda.resize(N, 3);

	int LDA = A.outerStride();
	int LDB = B.outerStride();
	int LDV = v.outerStride();

	double WORKDUMMY;
	int LWORK = -1; // Request optimum work size.
	int INFO = 0;

	double * alphar = const_cast<double*>(lambda.col(0).data());
	double * alphai = const_cast<double*>(lambda.col(1).data());
	double * beta = const_cast<double*>(lambda.col(2).data());

	// Get the optimum work size.
	dggev_("N", "V", &N, A.data(), &LDA, B.data(), &LDB, alphar, alphai, beta, 0, &LDV, v.data(), &LDV, &WORKDUMMY, &LWORK, &INFO);

	LWORK = int(WORKDUMMY) + 32;
	VectorXd WORK(LWORK);

	dggev_("N", "V", &N, A.data(), &LDA, B.data(), &LDB, alphar, alphai, beta, 0, &LDV, v.data(), &LDV, WORK.data(), &LWORK, &INFO);

	Derived col1 = lambda.col(0);
	Derived col3 = lambda.col(2);
	Derived d = col1.array() / col3.array();

	return std::make_pair(v, d);
	}
	*/

}

#endif