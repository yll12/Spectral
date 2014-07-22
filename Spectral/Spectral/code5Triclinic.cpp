#include "code5Triclinic.h"
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

namespace Code5Triclinic
{
	void Code5TriclinicMethods::code5Triclinic(void)
	{
		double h = 2.5e-3;

		double rho = 8938.4;

		MatrixXcd c(6, 6);

		c(0, 0) = 20.787e10;
		c(0, 1) = 10.906e10;
		c(0, 2) = 9.341e10;
		c(0, 3) = 1.6574e10;
		c(0, 4) = -1.615e10;
		c(0, 5) = -2.3188e10;
		c(1, 1) = 16.77e10;
		c(1, 2) = 13.624e10;
		c(1, 3) = -2.4719e10;
		c(1, 4) = 1.128e10;
		c(1, 5) = 0.86831e10;
		c(2, 2) = 18.591e10;
		c(2, 3) = 0.81453e10;
		c(2, 4) = 0.82076e10;
		c(2, 5) = 1.4505e10;
		c(3, 3) = 10.023e10;
		c(3, 4) = 1.4505e10;
		c(3, 5) = 0.58388e10;
		c(4, 4) = 3.5110e10;
		c(4, 5) = 1.6574e10;
		c(5, 5) = 5.9472e10;

		double vt = sqrt(5.9472e10 / rho);

		int n = 70;

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXcd x = result.x.cast<complex<double>>();

		x = x * h;

		MatrixXcd d1 = pow(h, -1) * result.dm[0].cast<complex<double>>();
		MatrixXcd d2 = pow(h, -2) * result.dm[1].cast<complex<double>>();

		MatrixXcd o = MatrixXcd::Zero(n, n);

		MatrixXcd k = MatrixXcd::Zero(1, 1);

		int m = 0;

		complex<double> i(0, 1);

		MatrixXcd l11 = -pow(k(0, m), 2) * c(4, 4) * MatrixXcd::Identity(n, n) + 2.0 * c(4, 5) * i * k(0, m) * d1 + c(5, 5) * d2;
		MatrixXcd l12 = -pow(k(0, m), 2) * c(3, 4) * MatrixXcd::Identity(n, n) + (c(3, 5) + c(1, 4)) * i * k(0, m) * d1 + c(1, 5) * d2;
		MatrixXcd l13 = -pow(k(0, m), 2) * c(2, 4) * MatrixXcd::Identity(n, n) + (c(2, 5) + c(3, 4)) * i * k(0, m) * d1 + c(3, 5) * d2;

		MatrixXcd l22 = -pow(k(0, m), 2) * c(3, 3) * MatrixXcd::Identity(n, n) + 2.0 * c(1, 3) * i * k(0, m) * d1 + c(1, 1) * d2;
		MatrixXcd l23 = -pow(k(0, m), 2) * c(2, 3) * MatrixXcd::Identity(n, n) + (c(1, 2) + c(3, 3)) * i * k(0, m)* d1 + c(1, 3) * d2;

		MatrixXcd l33 = -pow(k(0, m), 2) * c(2, 2) * MatrixXcd::Identity(n, n) + 2.0 * c(2, 3) * i * k(0, m) * d1 + c(3, 3) * d2;

		MatrixXcd lp;

		lp << l11, l12, l13,
			l12, l22, l23,
			l13, l23, l33;

		MatrixXcd l = lp;

		MatrixXcd s11 = c(1, 5) * d1 + i * k(0, m) * c(1, 4) * MatrixXcd::Identity(n, n);
		MatrixXcd s12 = c(1, 1) * d1 + i * k(0, m) * c(1, 3) * MatrixXcd::Identity(n, n);
		MatrixXcd s13 = c(1, 3) * d1 + i * k(0, m) * c(1, 2) * MatrixXcd::Identity(n, n);

		MatrixXcd s21 = c(3, 5) * d1 + i * k(0, m) * c(3, 4) * MatrixXcd::Identity(n, n);
		MatrixXcd s22 = c(1, 3) * d1 + i * k(0, m) * c(3, 3) * MatrixXcd::Identity(n, n);
		MatrixXcd s23 = c(3, 3) * d1 + i * k(0, m) * c(2, 3) * MatrixXcd::Identity(n, n);

		MatrixXcd s31 = c(5, 5) * d1 + i * k(0, m) * c(4, 5) * MatrixXcd::Identity(n, n);
		MatrixXcd s32 = c(1, 5) * d1 + i * k(0, m) * c(3, 5) * MatrixXcd::Identity(n, n);
		MatrixXcd s33 = c(3, 5) * d1 + i * k(0, m) * c(2, 5) * MatrixXcd::Identity(n, n);

		MatrixXcd s;
		
		s << s31, s32, s33,
			s11, s12, s13,
			s21, s22, s23;

		l.row(0) = s.row(0);
		l.row(n - 1) = s.row(n - 1);

		l.row(n) = s.row(n);
		l.row(2 * n - 1) = s.row(2 * n - 1);

		l.row(2 * n) = s.row(2 * n);
		l.row(3 * n - 1) = s.row(3 * n - 1);

		MatrixXcd mp = -rho*MatrixXcd::Identity(n, n); 
		mp(0, 0) = 0; 
		mp(n - 1, n - 1) = 0;

		MatrixXcd m2;

		m2 << mp, o, o,
			o, mp, o,
			o, o, mp;

		MatrixXcd p, e;

		pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_ceig(l, m2);

		p = eig.first;
		e = eig.second;

		MatrixXcd w = e.diagonal().cwiseSqrt();

	}

}