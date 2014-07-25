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

// Testing
#include <cfloat>
#include <iostream>
#include <sstream>
#include <cmath>

string FPClass(double x)
{
	int i = _fpclass(x);
	string s;
	switch (i)
	{
	case _FPCLASS_SNAN: s = "Signaling NaN";                break;
	case _FPCLASS_QNAN: s = "Quiet NaN";                    break;
	case _FPCLASS_NINF: s = "Negative infinity (-INF)";     break;
	case _FPCLASS_NN:   s = "Negative normalized non-zero"; break;
	case _FPCLASS_ND:   s = "Negative denormalized";        break;
	case _FPCLASS_NZ:   s = "Negative zero (-0)";           break;
	case _FPCLASS_PZ:   s = "Positive 0 (+0)";              break;
	case _FPCLASS_PD:   s = "Positive denormalized";        break;
	case _FPCLASS_PN:   s = "Positive normalized non-zero"; break;
	case _FPCLASS_PINF: s = "Positive infinity (+INF)";     break;
	}
	return s;
}

string HexDump(double x)
{
	unsigned long* pu;
	pu = (unsigned long*)&x;
	ostringstream os;
	os << hex << pu[0] << " " << pu[1];
	return os.str();
}

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

		MatrixXd temp(n, 1);

		temp.fill(1);

		x = (x + temp) * (h / 2);

		MatrixXd d1 = pow((h / 2), -1)*result.dm[0];
		MatrixXd d2 = pow((h / 2), -2)*result.dm[1];

		UtilityMethods::eigenToCSV(d1, "../../d1code11cpp.csv");
		UtilityMethods::eigenToCSV(d2, "../../d2code11cpp.csv");

		MatrixXd o = MatrixXd::Zero(n, n);

		// Set up m for loop

		int m = 0;

		MatrixXd lp = MatrixXd::Zero(n, n);

		lp = -pow(k(0, m), 2) * c(4, 4) * MatrixXd::Identity(n, n) + c(5, 5) * d2;

		UtilityMethods::eigenToCSV(lp, "../../lp11cpp.csv");

		MatrixXd l(n, n);

		l = lp;

		MatrixXd s = c(5, 5) * d1;

		UtilityMethods::eigenToCSV(s, "../../s11cpp.csv");

		l.row(0) = s.row(0);
		l.row(n - 1) = s.row(n - 1);

		MatrixXd m2 = MatrixXd::Identity(n, n);
		
		m2 *= -rho;
		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n; j++) {
				if (i != j) {
					m2(i, j) = abs(m2(i, j));
				}
			}
		}
		
		m2(0, 0) = 0;
		m2(n - 1, n - 1) = 0;

		UtilityMethods::eigenToCSV(l, "../../l11cpp.csv");
		UtilityMethods::eigenToCSV(m2, "../../m211cpp.csv");

		cout << m2 << endl;

		GeneralizedSelfAdjointEigenSolver<MatrixXd> solver(l, m2);

		MatrixXd p = solver.eigenvectors();
		MatrixXd e = solver.eigenvalues();

		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n; j++) {
				cout << HexDump(m2(i, j)) << " _fpclass(z) = " << FPClass(m2(i, j)) << "\n";
			}
		}

		cout << "p: \n" << endl;

		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < n; j++) {
				cout << HexDump(p(i, j)) << " _fpclass(z) = " << FPClass(p(i, j)) << "\n";
			}
		}

		cout << "e: \n" << endl;

		for (size_t i = 0; i < n; i++)
		{
			cout << HexDump(e(i, 0)) << " _fpclass(z) = " << FPClass(e(i, 0)) << "\n";
		}

		MatrixXd w = e.diagonal().cwiseSqrt();

		UtilityMethods::eigenToCSV(p, "../../p11cpp.csv");
		UtilityMethods::eigenToCSV(e, "../../e11cpp.csv");
		UtilityMethods::eigenToCSV(w, "../../w11cpp.csv");

	}

}