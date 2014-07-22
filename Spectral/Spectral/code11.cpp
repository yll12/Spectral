#include "code11.h"
#include "spectral.h"
#include "derivativeMatrix.h"
#include "utils.h"
#include <Eigen\Dense>
#include <cmath>
#include <utility>

//Delete later
#include <iostream>
using namespace Utility;

using namespace Eigen;
using namespace Spectral;
using namespace std;

namespace Code11
{
	void Code11Method::code11(void)
	{
		double h = 5e-3;
		double rho = 5300;

		MatrixXd c(6, 6);

		c(0, 0) = 23.9e10;
		c(0, 1) = 10.4e10;
		c(0, 2) = 5e10;
		c(1, 1) = 24.7e10;
		c(1, 2) = 5.2e10;
		c(2, 2) = 13.5e10;
		c(3, 3) = 6.5e10;
		c(4, 4) = 6.6e10;
		c(5, 5) = 7.6e10;

		cout << c << endl;

		double vt = sqrt(c(5, 5) / rho);

		MatrixXd k = MatrixXd::Zero(1,1);
		int n = 70;

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXd x = result.x;

		MatrixXd temp(n, 1);

		temp.fill(1);

		x = (x + temp) * (h / 2);

		MatrixXd d1 = result.dm[0];
		MatrixXd d2 = result.dm[1];

		MatrixXd o = MatrixXd::Zero(n, n);

		// Set up m for loop

		int m = 0;

		MatrixXd lp(n, n);

		lp = -pow(k(0, m), 2) * c(4, 4) * MatrixXd::Identity(n, n) + c(5, 5) * d2;

		MatrixXd l = lp;

		MatrixXd s = c(5, 5) * d1;

		l.row(0) = s.row(0);
		l.row(n - 1) = s.row(n - 1);

		MatrixXd m2(n, n);
		
		m2 = rho * MatrixXd::Identity(n, n);
		m2(0, 0) = 0;
		m2(n - 1, n - 1) = 0;

		MatrixXd p, e;

		std::pair<MatrixXd, MatrixXd> eig = UtilityMethods::matlab_eig(l, m2);

		p = eig.first;
		e = eig.second;

		MatrixXd w = e.diagonal().cwiseSqrt();

		UtilityMethods::eigenToCSV(p, "../../p11cpp.csv");
		UtilityMethods::eigenToCSV(e, "../../e11cpp.csv");
		UtilityMethods::eigenToCSV(w, "../../w11cpp.csv");

	}

}