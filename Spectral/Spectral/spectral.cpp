#define _USE_MATH_DEFINES
#include <cmath>

#include "spectral.h"
#include "utils.h"
#include <stdexcept>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

// Delete later
#include <fstream>
#include <iomanip>

using namespace std;
using namespace Eigen;
using namespace Utility;

namespace Spectral
{
	DerivativeMatrix SpectralMethods::chebdif(int n, int m)
	{
	
		int n1 = (int) floor(n / 2.0); // Indices used for flipping trick.	
		int n2 = (int) ceil(n / 2.0);
	
		MatrixXd th(n, 1); // Compute theta vector.

		for (int j = 0; j < n; j++)
		{
		    th(j, 0) = j * M_PI / (n - 1);
		}

		MatrixXd x(n, 1); // Compute Chebyshev points.

		for (int index = 0, j = n - 1; j >= 1 - n; index++, j -= 2)
		{
			x(index, 0) = sin(M_PI * j / (2 * (n - 1)));
		}


		MatrixXd t(n, n);
		t = (th / 2.0).replicate(1, n);
		
		MatrixXd t_transpose(n, n);
		t_transpose = t.transpose();
		
		MatrixXd dx(n, n);

		// Trigonometric identity.
		// dx = (2 * (t_transpose + t).unaryExpr(ptr_fun(sin))).cwiseProduct((t_transpose - t).unaryExpr(ptr_fun(sin)));
		dx = 2 * (t_transpose + t).array().sin().cwiseProduct((t_transpose - t).array().sin());
		
		// Flipping trick
		// DX = [DX(1:n1,:); -flipud(fliplr(DX(1:n2,:)))];
		MatrixXd top = dx.topRows(n2);
		dx << dx.topRows(n1),
			-UtilityMethods::matlab_flipud(UtilityMethods::matlab_fliplr(top));
		
		// Put 1's on the main diagonal of dx.
		dx += MatrixXd::Identity(n, n);

		MatrixXd h(n, 1);
		h.fill(-1);

		for (int i = 0; i < n; i += 2)
		{
			h(i, 0) = -h(i, 0);
		}
		
		MatrixXd c(n, n); // C is the matrix with entries c(k) / c(j)
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

		MatrixXd z(n, n); // Z contains entries 1 / (x(k) - x(j)) with zeros on the diagonal.
		z = dx.array().inverse();

		z -= MatrixXd::Identity(n, n);

		MatrixXd d = MatrixXd::Identity(n, n); // D contains diff. matrices.
		vector<MatrixXd> _dm;
		_dm.resize(m);
		for (int i = 0; i < m; i++){
			_dm[i] = MatrixXd::Zero(n, n);
		}

		for (int i = 0; i < m; i++)
		{ 
			MatrixXd j = c.cwiseProduct(d.diagonal().replicate(1, n)) - d;
			d = (i + 1) * z.cwiseProduct(j); // Off-diagonals.
			MatrixXd b = (-(d.transpose().colwise().sum()));
			for (int q = 0; q < n; q++)
			{
				d(q, q) = b(0, q); // Correct main diagonal of D
			}
			_dm[i] = d; // Store current D in DM
		} 

		DerivativeMatrix result; // Contains both matrix for the differentiation matrices and chebyshev points

		result.dm = _dm;
		result.x = x;

		return result;

	}

}