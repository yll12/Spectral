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

		dx = (2 * (t_transpose + t).unaryExpr(ptr_fun(sin))).cwiseProduct((t_transpose - t).unaryExpr(ptr_fun(sin)));

		dx += MatrixXd::Identity(n, n);

		MatrixXd h(n, 1);

		h.fill(-1);

		for (int i = 0; i < n; i++)
		{
			h(i, 0) = pow(h(i, 0), i);
		}
		
		MatrixXd c(n, n);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				c(i, j) = h(abs(i - j));
			}

		}

		c.row(0) *= 2;
		c.row(n-1) *= 2;
		c.col(0) /= 2;
		c.col(n-1) /= 2;

		MatrixXd z(n, n);

		z = dx.cwiseInverse();

		z -= MatrixXd::Identity(n, n);

		MatrixXd d = MatrixXd::Identity(n, n);

		MatrixXd *dm = new MatrixXd[m];
		for (int i = 0; i < m; i++){
			dm[i] = MatrixXd::Zero(n, n);
		}

		for (int i = 0; i < m; i++)
		{ 
			MatrixXd j = c.cwiseProduct(d.diagonal().replicate(1, n)) - d;
			d = (i + 1) * z.cwiseProduct(j);
			MatrixXd x = (-(d.transpose().colwise().sum()));
			for (int j = 0; j < n; j++)
			{
				d(j, j) = x(0, j);
			}
			dm[i] = d;
		}

		for (int i = 0; i < m; i++)
		{
			cout << "Here's the matrix dm[" << i << "]: \n" << dm[i] << "\n";
		}

		delete[] dm;

		return n + m;

	}

}