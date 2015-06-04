#include "Utils.h"

#include <fstream>
#include <iostream>
#include <iomanip>

namespace Utility
{
	bool UtilityMethods::compare(const double& x, const double &y)
	{
		return !(std::isinf(x) || std::isnan(x)) && ((std::isinf(y) || std::isnan(y)) || x < y);
	}

	std::vector<Eigen::MatrixXcd> UtilityMethods::compute_EOM(Eigen::MatrixXcd D,
		Eigen::MatrixXcd N, double h, double k) 
	{
		int G = D.rows();
		//Q0
		Eigen::MatrixXcd Q12X(G, G);
		Eigen::MatrixXcd Q12Y(G, G);
		Eigen::MatrixXcd Q12Z(G, G);
		
		//Q1
		Eigen::MatrixXcd Q02X(G, G);
		Eigen::MatrixXcd Q02Y(G, G);
		Eigen::MatrixXcd Q02Z(G, G);
		
		std::complex<double> i(0, 1);
		
		Q12X = N(1, 4) * i * Eigen::MatrixXcd::Identity(G, G);
		Q12Y = N(1, 3) * i * Eigen::MatrixXcd::Identity(G, G);
		Q12Z = N(1, 2) * i * Eigen::MatrixXcd::Identity(G, G);
		
		Q02X = N(1, 5) * D;
		Q02Y = N(1, 1) * D;
		Q02Z = N(1, 3) * D;
		
		std::vector<Eigen::MatrixXcd> Q0(6);
		std::vector<Eigen::MatrixXcd> Q1(6);
		Q0[1] = Eigen::MatrixXcd::Zero(G, 3 * G);
		Q1[1] = Eigen::MatrixXcd::Zero(G, 3 * G);
		Q0[1] << Q02X, Q02Y, Q02Z;
		Q1[1] << Q12X, Q12Y, Q12Z;
		
		std::vector<Eigen::MatrixXcd> T(6);
		T[1] = (h * h * k) * Q1[1] + h * Q0[1];
		return T;
	}

	Eigen::MatrixXd UtilityMethods::normalise(Eigen::MatrixXd c) 
	{
		double denominator = c(5, 5);
		return c / denominator;
	}

	Eigen::MatrixXd UtilityMethods::rotate_YToZ(Eigen::MatrixXd c) 
	{
		Eigen::MatrixXd c_rotated(6,6);
		c_rotated(0, 0) = c(1, 1);
		c_rotated(0, 1) = c(1, 2);
		c_rotated(0, 2) = c(0, 1);
		c_rotated(0, 3) = c(1, 4);
		c_rotated(0, 4) = c(1, 5);
		c_rotated(0, 5) = c(1, 3);
		c_rotated(1, 1) = c(2, 2);
		c_rotated(1, 2) = c(0, 2);
		c_rotated(1, 3) = c(2, 4);
		c_rotated(1, 4) = c(2, 5);
		c_rotated(1, 5) = c(2, 3);
		c_rotated(2, 2) = c(0, 0);
		c_rotated(2, 3) = c(0, 4);
		c_rotated(2, 4) = c(0, 5);
		c_rotated(2, 5) = c(0, 3);
		c_rotated(3, 3) = c(4, 4);
		c_rotated(3, 4) = c(4, 5);
		c_rotated(3, 5) = c(3, 4);
		c_rotated(4, 4) = c(5, 5);
		c_rotated(4, 5) = c(3, 5);
		c_rotated(5, 5) = c(3, 3);
		for (int i = 1; i < c.rows(); i++)
		{
			for (int j = 0; j < i; j++)
			{
				c_rotated(i, j) = c_rotated(j, i);
			}
		}
		return c_rotated;
	}

	Eigen::MatrixXd UtilityMethods::rotate_YToX(Eigen::MatrixXd c)
	{
		Eigen::MatrixXd c_rotated(6, 6);
		c_rotated(0, 0) = c(2, 2);
		c_rotated(0, 1) = c(0, 2);
		c_rotated(0, 2) = c(1, 2);
		c_rotated(0, 3) = c(2, 5);
		c_rotated(0, 4) = c(2, 3);
		c_rotated(0, 5) = c(2, 4);
		c_rotated(1, 1) = c(0, 0);
		c_rotated(1, 2) = c(0, 1);
		c_rotated(1, 3) = c(0, 5);
		c_rotated(1, 4) = c(0, 3);
		c_rotated(1, 5) = c(0, 4);
		c_rotated(2, 2) = c(1, 1);
		c_rotated(2, 3) = c(1, 5);
		c_rotated(2, 4) = c(1, 3);
		c_rotated(2, 5) = c(1, 4);
		c_rotated(3, 3) = c(5, 5);
		c_rotated(3, 4) = c(3, 5);
		c_rotated(3, 5) = c(4, 5);
		c_rotated(4, 4) = c(3, 3);
		c_rotated(4, 5) = c(3, 4);
		c_rotated(5, 5) = c(4, 4);
		for (int i = 1; i < c.rows(); i++)
		{
			for (int j = 0; j < i; j++)
			{
				c_rotated(i, j) = c_rotated(j, i);
			}
		}
		return c_rotated;
	}

	std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> UtilityMethods::chebdif(int n, int m)
	{

		int n1 = (int)floor(n / 2.0); // Indices used for flipping trick.	
		int n2 = (int)ceil(n / 2.0);

		Eigen::MatrixXd th(n, 1); // Compute theta vector.

		for (int j = 0; j < n; j++)
		{
			th(j, 0) = j * M_PI / (n - 1);
		}

		Eigen::MatrixXd x(n, 1); // Compute Chebyshev points.

		for (int index = 0, j = n - 1; j >= 1 - n; index++, j -= 2)
		{
			x(index, 0) = sin(M_PI * j / (2 * (n - 1)));
		}

		Eigen::MatrixXd t(n, n);
		t = (th / 2.0).replicate(1, n);

		Eigen::MatrixXd t_transpose(n, n);
		t_transpose = t.transpose();

		Eigen::MatrixXd dx(n, n);

		// Trigonometric identity.
		// dx = (2 * (t_transpose + t).unaryExpr(ptr_fun(sin))).cwiseProduct((t_transpose - t).unaryExpr(ptr_fun(sin)));
		dx = 2 * (t_transpose + t).array().sin().cwiseProduct((t_transpose - t).array().sin());

		// Flipping trick
		// DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];
		Eigen::MatrixXd top = dx.topRows(n2);
		dx << dx.topRows(n1),
			-UtilityMethods::matlab_flipud(UtilityMethods::matlab_fliplr(top));

		// Put 1's on the main diagonal of dx.
		dx += Eigen::MatrixXd::Identity(n, n);

		Eigen::MatrixXd h(n, 1);
		h.fill(-1);

		for (int i = 0; i < n; i += 2)
		{
			h(i, 0) = -h(i, 0);
		}

		Eigen::MatrixXd c(n, n); // C is the matrix with entries c(k) / c(j)
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				c(i, j) = h(abs(i - j));
			}

		}

		c.row(0) *= 2;
		c.row(n - 1) *= 2;
		c.col(0) /= 2;
		c.col(n - 1) /= 2;

		Eigen::MatrixXd z(n, n); // Z contains entries 1 / (x(k) - x(j)) with zeros on the diagonal.
		z = dx.array().inverse();

		z -= Eigen::MatrixXd::Identity(n, n);

		Eigen::MatrixXd d = Eigen::MatrixXd::Identity(n, n); // D contains diff. matrices.
		std::vector<Eigen::MatrixXd> _dm;
		_dm.resize(m);
		for (int i = 0; i < m; i++){
			_dm[i] = Eigen::MatrixXd::Zero(n, n);
		}

		for (int i = 0; i < m; i++)
		{
			Eigen::MatrixXd j = c.cwiseProduct(d.diagonal().replicate(1, n)) - d;
			d = (i + 1) * z.cwiseProduct(j); // Off-diagonals.
			Eigen::MatrixXd b = (-(d.transpose().colwise().sum()));
			for (int q = 0; q < n; q++)
			{
				d(q, q) = b(0, q); // Correct main diagonal of D
			}
			_dm[i] = d; // Store current D in DM
		}

		return make_pair(x, _dm);

	}

	// OLD Implementations:
	/*
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
	*/
}

