#include "code11.h"
#include "spectral.h"
#include "derivativeMatrix.h"
#include "utils.h"
#include <Eigen\Dense>
#include <Eigen\Eigenvalues>
#include <cmath>
#include <utility>

//Delete later
#include <iostream>

using namespace Utility;
using namespace Eigen;
using namespace Spectral;
using namespace std;

// Testing
#include <cfloat>
#include <iostream>
#include <sstream>

namespace Code11
{
	
	void Code11Method::code11(void)
	{
		double h = 5e-3;
		double rho = 5300;

		MatrixXd c = MatrixXd::Zero(6, 6);

		c(0, 0) = 23.9e10;
		c(0, 1) = 10.4e10;
		c(0, 2) = 5e10;
		c(1, 1) = 24.7e10;
		c(1, 2) = 5.2e10;
		c(2, 2) = 13.5e10;
		c(3, 3) = 6.5e10;
		c(4, 4) = 6.6e10;
		c(5, 5) = 7.6e10;

		double vt = sqrt(c(5, 5) / rho);

		MatrixXd k = MatrixXd::Zero(1,1);
		int n = 10;

		DerivativeMatrix result = SpectralMethods::chebdif(n, 2);

		MatrixXd x = result.x;

		x = (x.array() + 1).matrix() * (h / 2);

		MatrixXd d1 = pow((h / 2), -1)*result.dm[0];
		MatrixXd d2 = pow((h / 2), -2)*result.dm[1];

		MatrixXd o = MatrixXd::Zero(n, n);

		// Set up m for loop

		int m = 0;

		MatrixXd lp = MatrixXd::Zero(n, n);

		lp = -pow(k(0, m), 2) * c(4, 4) * MatrixXd::Identity(n, n) + c(5, 5) * d2;

		MatrixXd l(n, n);

		l = lp;

		MatrixXd s = c(5, 5) * d1;

		l.row(0) = s.row(0);
		l.row(n - 1) = s.row(n - 1);

		MatrixXd m2 = MatrixXd::Identity(n, n);
		
		m2 *= -rho;
		m2(0, 0) = 0;
		m2(n - 1, n - 1) = 0;

		pair<MatrixXcd, MatrixXcd> eigs = Utility::UtilityMethods::matlab_eig(l, m2);
		MatrixXd p = eigs.first.real();
		MatrixXd e = eigs.second.real();

		cout << "p: \n" << endl;
		cout << p << endl;

		cout << "e: \n" << endl;
		cout << e << endl;

		MatrixXd w = e.cwiseSqrt();

		std::sort(w.data(), w.data() + w.size());

		UtilityMethods::eigenToCSV(p, "../../p11cpp_n10.csv");
		UtilityMethods::eigenToCSV(e, "../../e11cpp_n10.csv");
		UtilityMethods::eigenToCSV(w, "../../w11cpp_n10.csv");

	}

}