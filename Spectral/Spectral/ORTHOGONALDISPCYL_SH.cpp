#include "ORTHOGONALDISPCYL_SH.h"
#include "Utils.h"
#include <complex>
#include <Eigen\Dense>

//Delete later
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace Utility;

namespace OrthogonalDispcylsh
{
	void OrthogonalDispcylshMethods::orthogonalDispcylsh(void)
	{

		double a = 5;
		double b = 5.005;
		double rho = 5300;

		MatrixXd c;

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

		double vt = sqrt(c(5, 5) / rho);

		std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> result = UtilityMethods::chebdif(n, 2);

		MatrixXd x = result.first;

		double h = b - a;

		MatrixXd a_temp(n, 1);

		a_temp.fill(a);

		MatrixXd b_temp(n, 1);

		b_temp.fill(b);

		MatrixXd r = (h * x + a_temp + b_temp) / 2;

		std::vector<Eigen::MatrixXd> dm = result.second;
		MatrixXd d1 = (2 / h) * dm[0];
		MatrixXd d2 = pow((2 / h), 2) * dm[1];

		MatrixXd o = MatrixXd::Zero(n, n);

		MatrixXd k = MatrixXd::Zero(1, 1);

		int m = 0;

		MatrixXd lp = c(5, 5)*(r.array().pow(2).matrix().diagonal()*d2 + r.diagonal() * d1 - MatrixXd::Identity(n, n)) - c(3, 3) * r.array().pow(2).matrix().diagonal() * pow(k(0, m), 2);

		MatrixXd l = lp;

		MatrixXd s = c(5, 5)*(r.array().pow(2).matrix().diagonal() * d1 - r.diagonal());

		l.row(0) = s.row(0);
		l.row(n - 1) = l.row(n - 1);

		MatrixXd mp = -rho * (r.array().pow(2).matrix().diagonal());
		mp *= MatrixXd::Identity(n, n);
		mp(0, 0) = 0; 
		mp(n - 1, n - 1) = 0;

		MatrixXd m2 = mp;

		MatrixXd p, e;

		pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_eig(l, m2);

		p = eig.first.real();
		e = eig.second.real();

		MatrixXd w = e.diagonal().cwiseSqrt();
	}

}