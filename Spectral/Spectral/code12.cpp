#include "code12.h"
#include "spectral.h"
#include "Utils.h"
#include "derivativeMatrix.h"
#include <complex>
#include <Eigen\Dense>

//Delete later
#include <iostream>

using namespace Eigen;
using namespace Spectral;
using namespace std;
using namespace Utility;

namespace Code12
{
	void Code12Method::code12(void)
	{
		double h = 2.5e-3;

		double rho = 5300; 

		MatrixXcd c = MatrixXcd::Zero(6, 6);

		c(0, 0) = 23.9e10;
		c(0, 1) = 10.4e10;
		c(0, 2) = 5e10;
		c(1, 1) = 24.7e10;
		c(1, 2) = 5.2e10;
		c(2, 2) = 13.5e10;
		c(3, 3) = 6.5e10;
		c(4, 4) = 6.6e10;
		c(5, 5) = 7.6e10;

		complex<double> vt = sqrt(c(5, 5) / rho);

		int n = 90;

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXcd x = result.x.cast<complex<double>>();

		x = x * h;

		MatrixXcd d1 = pow(h, -1) * result.dm[0].cast<complex<double>>();

		MatrixXcd d2 = pow(h, -2) * result.dm[1].cast<complex<double>>();

		MatrixXcd o = MatrixXcd::Zero(n, n);

		int kmin = 0;
		int kmax = 35;

		const int steps = 4201;

		MatrixXcd k = MatrixXcd::Zero(1, steps);

		{
			double fixed = ((double)(kmax - kmin)) / (steps - 1);
			double i;
			int j;
			for (i = 0, j = 0; j <= 4.2e3; i += fixed, j++)
			{
				k(0, j) = i;
			}

		}

		k = k.cwiseProduct(k.cwiseProduct(k));

		MatrixXd W = MatrixXd::Zero(2 * n, k.cols());

		int m = 0;

		MatrixXcd l11 = d2 * c(1, 1) - pow(k(0, m), 2) * c(3, 3) * MatrixXcd::Identity(n, n);
		
		complex<double> i(0, 1);

		MatrixXcd l12 = d1 * i * k(0, m) * (c(1, 2) + c(3, 3));
		MatrixXcd l22 = -pow(k(0, m), 2) * c(2, 2) * MatrixXcd::Identity(n, n) + d2 * c(3, 3);

		MatrixXcd lp(l11.rows() + l12.rows(), l11.cols() + l12.cols());

		lp << l11, l12,
			l12, l22;

		MatrixXcd l = lp;

		MatrixXcd s11 = c(1, 1) * d1;
		MatrixXcd s12 = c(1, 2) * i * k(0, m)* MatrixXcd::Identity(n, n);
		MatrixXcd s21 = c(3, 3) * i * k(0, m)* MatrixXcd::Identity(n, n);
		MatrixXcd s22 = c(3, 3) * d1;

		MatrixXcd s(s11.rows() + s21.rows(), s11.cols() + s12.cols());

		s << s11, s12,
			s21, s22;

		l.row(0) = s.row(0);
		l.row(n - 1) = s.row(n - 1);
		l.row(n) = s.row(n);
		l.row(2 * n - 1) = s.row(2 * n - 1);

		MatrixXcd mp = -rho*MatrixXcd::Identity(n, n); 
		mp(0, 0) = 0; 
		mp(n - 1, n - 1) = 0;

		MatrixXcd m2(mp.rows() + o.rows(), mp.cols() + o.cols());
		
		m2 << mp, o,
			o, mp;

		MatrixXd p, e;

		pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_eig(l, m2);

		p = eig.first.real();
		e = eig.second.real();

		MatrixXd w = e.cwiseSqrt();

		std::sort(w.data(), w.data() + w.size(), UtilityMethods::compare);

		MatrixXd p_(2 * n, steps);
		p_.col(m) = w.col(0);
		int q = 0;

		for (int i = 0; i < 2 * n; i++)
		{
			if (p_(i, m) != 0)
			{
				W(q, m) = p_(i, m);
				q++;
			}
		}

		UtilityMethods::eigenToCSV(W.col(m), "../../bigw12cpp_n90.csv");

		//UtilityMethods::eigenToCSV(p.real(), "../../p12cpp_n120.csv");
		//UtilityMethods::eigenToCSV(e.real(), "../../e12cpp_n120.csv");
		//UtilityMethods::eigenToCSV(w_2, "../../w12cpp_n120.csv");

	}

}