#define _USE_MATH_DEFINES
#include <cmath>

#include "spectral.h"

#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

namespace Spectral
{
	double SpectralMethods::chebdif(int n, int m)
	{

		MatrixXd i = MatrixXd::Identity(n, n);

		double n1 = floor(n * (1.0 / 2.0));

		double n2 = ceil(n * (1.0 / 2.0));

		MatrixXd k(n, 1);

		for (int i = 0; i < n; i++)
		{
			k(i, 0) = i;
		}

		MatrixXd th(n, 1);

		th = k * M_PI / (n - 1);

		MatrixXd x(n - 1, 1);

		int index = 0;
		for (int i = n - 1; i > 1-n; i-=2)
		{
			x(index, 0) = sin(M_PI * i / (2 * (n - 1)));
			index++;
		}

		MatrixXd t(n, n);
		MatrixXd th_temp(n, 1);
		th_temp = th / 2;
		t = th_temp.replicate(1, n);

		MatrixXd t_transpose(n, n);
		t_transpose = t.transpose();

		MatrixXd dx(n, n);

		dx = (t_transpose + t);

		cout << "Here's the matrix dx: \n" << dx << "\n";

		dx = dx.sin();

		//dx = 2 * (t_transpose + t).sin().cwiseProduct((t_transpose - t).sin());

		cout << "Here's the matrix dx: \n" << dx << "\n";

		return n + m;
	}

}