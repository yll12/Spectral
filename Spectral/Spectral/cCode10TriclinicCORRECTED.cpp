#include "cCode10TriclinicCORRECTED.h"
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

namespace Code10Triclinic
{
	void Code10TriclinicMethods::code10triclinic(void)
	{

		double a = 5;
		double b = 5.005;

		double rho = 8938.4;

		MatrixXcd c = MatrixXcd::Zero(6, 6);

		c(0, 0) = 182.1e9;
		c(0, 1) = 119.1e9;
		c(0, 2) = 111.8e9;
		c(0, 3) = 26.2e9;
		c(0, 4) = -17.3e9;
		c(0, 5) = 0;
		c(1, 1) = 167.7e9;
		c(1, 2) = 126.2e9;
		c(1, 3) = -26.2e9;
		c(1, 4) = 17.3e9;
		c(1, 5) = 0;
		c(2, 2) = 174.9e9;
		c(2, 3) = 0;
		c(2, 4) = 4.237e9;
		c(2, 5) = 0;
		c(3, 3) = 92.1e9;
		c(3, 4) = 0;
		c(3, 5) = 17.3e9;
		c(4, 4) = 53.5e9;
		c(4, 5) = 26.2e9;
		c(5, 5) = 67.6e9;

		double alpha = M_PI / 9.3;

		double vt = sqrt(67.6e9 / rho);

		int n = 90;

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXcd x = result.x.cast<complex<double>>();

		double h = b - a;

		MatrixXcd r(n, 1);

		r = ((h*x).array() + a + b).matrix() / 2;
		UtilityMethods::eigenToCSV(r.real(), "../../r_code10_cpp_n70.csv");

		MatrixXcd d1 = (2 / h) * result.dm[0].cast<complex<double>>();
		MatrixXcd d2 = pow((2 / h), 2) * result.dm[1].cast<complex<double>>();

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

		complex<double> i(0, 1);

		MatrixXd W = MatrixXd::Zero(2 * n, k.cols());

		//for (int m = 0; m <= k.cols() - 101; m+= 100)
		//for (int m = 0; m <= 10; m++)
		{
			int m = 0;
			MatrixXd diag_r_pow_1 = r.array().pow(-1).matrix().asDiagonal().toDenseMatrix().real();
			MatrixXd diag_r_pow_2 = r.array().pow(-2).matrix().asDiagonal().toDenseMatrix().real();
			MatrixXcd l11 = c(1, 1) * (d2 + diag_r_pow_1 * d1) - c(0, 0)* diag_r_pow_2 * MatrixXcd::Identity(n, n) + i * k(0, m) * c(1, 3) * (2 * d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n)) - pow(k(0, m), 2) * c(3, 3) * MatrixXcd::Identity(n, n);

			//L12 = 1i * k(1, (m + 1))*C(2, 5)*(D1 + diag(r.^-1)*eye(N)) - 1i * k(1, (m + 1))*C(1, 5)*diag(r.^-1)*eye(N) - k(1, (m + 1)) ^ 2 * C(4, 5)*eye(N) + C(2, 6)*D2-C(1, 6)*diag(r.^-1)*(D1 - diag(r.^-1)*eye(N)) + 1i * k(1, (m + 1))*C(4, 6)*(D1 - diag(r.^-1)*eye(N));
			MatrixXcd l12 = i * k(0, m) * c(1, 4) * (d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n)) - i * k(0, m) * c(0, 4) * diag_r_pow_1 *MatrixXcd::Identity(n, n) -pow(k(0, m), 2) * c(3, 4) * MatrixXcd::Identity(n, n) + c(1, 5) * d2;
			l12 += -c(0, 5) * diag_r_pow_1 * (d1 - diag_r_pow_1 * MatrixXcd::Identity(n, n));
			l12 += i * k(0, m) * c(3, 5) * (d1 - diag_r_pow_1 * MatrixXcd::Identity(n, n));

			MatrixXcd l13 = i * k(0, m) * c(1, 2) * (d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n)) - i * k(0, m) * c(0, 2) * diag_r_pow_1 * MatrixXcd::Identity(n, n);
			l13 += -pow(k(0, m), 2) * c(2, 3) * MatrixXcd::Identity(n, n) + c(1, 3) * (d2 + diag_r_pow_1 * d1);
			l13 += -c(0, 3) * diag_r_pow_1 * d1 + i * k(0, m) * c(3, 3) * d1;

			MatrixXcd l21 = i * k(0, m) * c(1, 4) * d1 + c(1, 5) * (d2 + 2 * diag_r_pow_1 * d1) + i * k(0, m) * c(0, 4) * diag_r_pow_1 * MatrixXcd::Identity(n, n) + c(0, 5) * (diag_r_pow_1 * d1 + diag_r_pow_2 * MatrixXcd::Identity(n, n));
			l21 += -pow(k(0, m), 2) * c(3, 4) * MatrixXcd::Identity(n, n) + i * k(0, m) * c(3, 5) * (d1 + 2 * diag_r_pow_1 * MatrixXcd::Identity(n, n));
			
			MatrixXcd l22 = -pow(k(0, m), 2) * c(4, 4) * MatrixXcd::Identity(n, n) + i * k(0, m) * c(4, 5) * (2 * d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n)) + c(5, 5) * (d2 + diag_r_pow_1 * d1 - diag_r_pow_2 * MatrixXcd::Identity(n, n));
			MatrixXcd l23 = -pow(k(0, m), 2) * c(2, 4) * MatrixXcd::Identity(n, n) + i * k(0, m) * c(2, 5) * (d1 + 2 * diag_r_pow_1 * MatrixXcd::Identity(n, n)) + i * k(0, m) * c(3, 4) * d1 + c(3, 5) * (d2 + 2 * diag_r_pow_1 * d1);

			MatrixXcd l31 = i * k(0, m) * c(1, 2) * d1 + c(1, 3) * (d2 + diag_r_pow_1 * d1) + i * k(0, m) * c(0, 2) * diag_r_pow_1 * MatrixXcd::Identity(n, n) + c(0, 3) * (diag_r_pow_1 * d1);
			l31 += -pow(k(0, m), 2) * c(2, 3) * MatrixXcd::Identity(n, n) + i * k(0, m) * c(3, 3) * (d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n));
			MatrixXcd l32 = -pow(k(0, m), 2) * c(2, 4) * MatrixXcd::Identity(n, n) + i * k(0, m) * c(3, 4) * (d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n)) + i * k(0, m) * c(2, 5) * (d1 - diag_r_pow_1 * MatrixXcd::Identity(n, n)) + c(3, 5) * d2;
			MatrixXcd l33 = -pow(k(0, m), 2) * c(2, 2) * MatrixXcd::Identity(n, n) + i * k(0, m) * c(2, 3) * (d1 + diag_r_pow_1 * MatrixXcd::Identity(n, n)) + i * k(0, m) * c(2, 3) * d1 + c(3, 3) * (d2 + diag_r_pow_1 * d1);

			MatrixXcd lp(l11.rows() + l21.rows() + l31.rows(), l11.cols() + l12.cols() + l13.cols());

			lp << l11, l12, l13,
				l21, l22, l23,
				l31, l32, l33;

			MatrixXcd l = lp;
			
			MatrixXcd s11 = c(1, 1) * d1 + c(0, 1) * diag_r_pow_1 * MatrixXcd::Identity(n, n) + i * k(0, m) * c(1, 3) * MatrixXcd::Identity(n, n);
			MatrixXcd s12 = i * k(0, m) * c(1, 4) * MatrixXcd::Identity(n, n) + c(1, 5) * (d1 - diag_r_pow_1 * MatrixXcd::Identity(n, n));
			MatrixXcd s13 = i * k(0, m) * c(1, 2) * MatrixXcd::Identity(n, n) + c(1, 3) * d1;

			MatrixXcd s21 = c(1, 3) * d1 + c(0, 3) * diag_r_pow_1 * MatrixXcd::Identity(n, n) + i * k(0, m) * c(3, 3) * MatrixXcd::Identity(n, n);
			MatrixXcd s22 = i * k(0, m) * c(3, 4) * MatrixXcd::Identity(n, n) + c(3, 5) * (d1 - diag_r_pow_1 * MatrixXcd::Identity(n, n));
			MatrixXcd s23 = i * k(0, m) * c(2, 3) * MatrixXcd::Identity(n, n) + c(3, 3) * d1;

			MatrixXcd s31 = c(1, 5) * d1 + c(0, 5) * diag_r_pow_1 * MatrixXcd::Identity(n, n) + i * k(0, m) * c(3, 5) * MatrixXcd::Identity(n, n);
			MatrixXcd s32 = i * k(0, m) * c(4, 5) * MatrixXcd::Identity(n, n) + c(5, 5) * (d1 - diag_r_pow_1 * MatrixXcd::Identity(n, n));
			MatrixXcd s33 = i * k(0, m) * c(2, 5) * MatrixXcd::Identity(n, n) + c(3, 5) * d1;
			
			MatrixXcd s(s11.rows() + s21.rows() + s31.rows(), s11.cols() + s12.cols() + s13.cols());


			//S = [S11, S12, S13; S21, S22, S23; S31, S32, S33];
			s << s11, s12, s13,
				s21, s22, s23,
				s31, s32, s33;

			l.row(0) = s.row(0);
			l.row(n - 1) = s.row(n - 1);

			l.row(n) = s.row(n);
			l.row(2 * n - 1) = s.row(2 * n - 1);

			l.row(2 * n) = s.row(2 * n);
			l.row(3 * n - 1) = s.row(3 * n - 1);

			MatrixXcd mp = -rho  *MatrixXcd::Identity(n, n);
			mp(0, 0) = 0;
			mp(n - 1, n - 1) = 0;

			MatrixXcd m2(mp.rows() + o.rows() + o.rows(), mp.cols() + o.cols() + o.cols());

			m2 << mp, o, o,
				o, mp, o,
				o, o, mp;

			MatrixXd p, e;

			pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_eig(l, m2);

			p = eig.first.real();
			e = eig.second.real();

			MatrixXd w = e.cwiseSqrt();

			std::sort(w.data(), w.data() + w.size(), UtilityMethods::compare);

			MatrixXd p_(3 * n, steps);
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
			
			//UtilityMethods::eigenToCSV(W.col(m), "../../bigw5cpp_n90.csv");
		
		}


		//UtilityMethods::eigenToCSV(p, "../../p5cpp_n120.csv");
		//UtilityMethods::eigenToCSV(e, "../../e5cpp_n120.csv");
		//UtilityMethods::eigenToCSV(w, "../../w5cpp_n120.csv");
		
			
	}

}