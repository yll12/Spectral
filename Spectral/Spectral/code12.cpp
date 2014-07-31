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

		int n = 120;

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXcd x = result.x.cast<complex<double>>();

		x = x * h;

		MatrixXcd d1 = pow(h, -1) * result.dm[0].cast<complex<double>>();
		MatrixXcd d2 = pow(h, -2) * result.dm[1].cast<complex<double>>();

		MatrixXcd o = MatrixXcd::Zero(n, n);

		MatrixXcd k = MatrixXcd::Zero(1, 1);

		int m = 0;

		MatrixXcd l11 = d2 * c(2, 2) - pow(k(0, m), 2) * c(3, 3) * MatrixXcd::Identity(n, n);
		
		complex<double> i(0, 1);

		MatrixXcd l12 = d1 * i * k(0, m ) * (c(1, 2) + c(3, 3));
		MatrixXcd l22 = -pow(k(0, m), 2) * c(2, 2) * MatrixXcd::Identity(n, n) + d2 * c(4, 4);

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

		MatrixXcd p, e;

		pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_ceig(l, m2);

		p = eig.first;
		e = eig.second;

		MatrixXd w = e.cwiseSqrt().real();

		std::sort(w.data(), w.data() + w.size());

		UtilityMethods::eigenToCSV(p.real(), "../../p12cpp_n120.csv");
		UtilityMethods::eigenToCSV(e.real(), "../../e12cpp_n120.csv");
		UtilityMethods::eigenToCSV(w, "../../w12cpp_n120.csv");

	}

}