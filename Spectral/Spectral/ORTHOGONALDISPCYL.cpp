#include "ORTHOGONALDISPCYL.h"
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

namespace OrthogonalDispcyl
{
	void OrthogonalDispcylMethods::orthogonalDispcyl(void)
	{

		double a = 5;
		double b = 5.005;
		double rho = 5300;

		MatrixXcd c;

		c(0, 0) = 23.9e9;
		c(0, 1) = 10.4e9;
		c(0, 2) = 5e9;
		c(1, 1) = 24.7e9;
		c(1, 2) = 5.2e9;
		c(2, 2) = 13.5e9;
		c(3, 3) = 6.5e9;
		c(4, 4) = 6.6e9;
		c(5, 5) = 7.6e9;

		int n = 100;

		complex<double> vt = sqrt(c(5, 5) / rho);

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXcd x = result.x.cast<complex<double>>();

		double h = b - a;

		MatrixXcd a_temp(n, 1);

		a_temp.fill(a);

		MatrixXcd b_temp(n, 1);

		b_temp.fill(b);

		MatrixXcd r = (h * x + a_temp + b_temp) / 2;

		MatrixXcd d1 = (2 / h) * result.dm[0].cast<complex<double>>();
		MatrixXcd d2 = pow((2 / h), 2) * result.dm[1].cast<complex<double>>();

		MatrixXcd o = MatrixXcd::Zero(n, n);

		MatrixXcd k = MatrixXcd::Zero(1, 1);

		int m = 0;

		complex<double> i(0, 1);

		MatrixXcd l11 = c(0, 0) * (d2 + r.array().pow(-1).matrix().diagonal() * d1) - pow(k(0, m), 2) * c(4, 4) * MatrixXcd::Identity(n, n) - c(1, 1) * r.array().pow(-2).matrix().diagonal();
		MatrixXcd l12 = i * k(0, m) * ((c(0, 2) + c(4, 4)) * d1 + (c(0, 2) - c(1, 2))*r.array().pow(-1).matrix().diagonal());

		MatrixXcd l21 = i * k(0, m) * (c(0, 2) * d1 + c(1, 2) * r.array().pow(-1).matrix().diagonal() + c(4, 4)*(d1 + r.array().pow(-1).matrix().diagonal()));
		MatrixXcd l22 = c(4, 4)*(d2 + r.array().pow(-1).matrix().diagonal() * d1) - pow(k(0, m), 2) * c(2, 2)*MatrixXcd::Identity(n, n);

		MatrixXcd lp; 
		
		lp << l11, l12, 
			l21, l22;

		MatrixXcd l = lp;

		MatrixXcd s11 = (c(0, 0) * d1 + c(0, 1) * r.array().pow(-1).matrix().diagonal());
		MatrixXcd s12 = (c(0, 2) * i * k(0, m) * MatrixXcd::Identity(n, n));
		MatrixXcd s21 = (i * k(0, m)*c(4, 4) * MatrixXcd::Identity(n, n));
		MatrixXcd s22 = c(4, 4)*d1;
		
		MatrixXcd s;
			
		s << s11, s12,
		  s21, s22;

		l.row(0) = s.row(0);
		l.row(n - 1) = s.row(n - 1);
		l.row(n) = s.row(n);
		l.row(2 * n - 1) = s.row(2 * n - 1);

		MatrixXcd mp = -rho*MatrixXcd::Identity(n, n); 
		mp(0, 0) = 0; 
		mp(n - 1, n - 1) = 0;

		MatrixXcd m2;

		m2 << mp, o,
			o, mp;

		MatrixXcd p, e;

		pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_ceig(l, m2);

		p = eig.first;
		e = eig.second;

		MatrixXcd w = e.diagonal().cwiseSqrt();
	}

}