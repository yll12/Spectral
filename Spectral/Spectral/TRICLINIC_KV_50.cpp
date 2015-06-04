#include "TRICLINIC_KV_50.h"
#include "Utils.h"
#include <complex>
#include <Eigen\Dense>

using namespace Eigen;
using namespace std;
using namespace Utility;

namespace TRICLINIC_KV_50
{
	void TRICLINIC_KV_50Methods::TRICLINIC_KV_50(void)
	{
		/*
		double h = 50e-3; //plate thickness in m
		double rho = 1500;  //kg / m3;
		double fn = 2e6;   //in Hz

		MatrixXcd e;

		e(0, 0) = 2.18e8;
		e(0, 1) = 7.65e7;
		e(0, 2) = 1.64e7;
		e(0, 3) = -3.60e6;
		e(0, 4) = 6.88e5;
		e(0, 5) = 1.16e8;

		e(1, 1) = 7.11e7;
		e(1, 2) = 1.92e7;
		e(1, 3) = -7.71e5;
		e(1, 4) = 2.15e6;
		e(1, 5) = 5.00e7;

		e(2, 2) = 4.22e7;
		e(2, 3) = -9.644e5;
		e(2, 4) = 6.27e5;
		e(2, 5) = -3.07e6;

		e(3, 3) = 1.11e7;
		e(3, 4) = 2.89e6;
		e(3, 5) = -1.15e6;

		e(4, 4) = 1.36e7;
		e(4, 5) = 1.48e6;

		e(5, 5) = 9.35e7;

		//%%%%%%%%%%%%%%%%% VELOCITY TO RESCALE THE PLOT%%%%%%%%%%%%%%%%%
		double Vt = sqrt(28.29e9 / rho);

		int n = 100; //checked for 80 and 120

		std::pair<Eigen::MatrixXd, std::vector<Eigen::MatrixXd>> result = UtilityMethods::chebdif(n, 2);
		MatrixXcd x = result.first.cast<complex<double>>();
		x = ((x.array() + 1) / 2) * h; //y
		x = x / h; //y hat
		
		std::vector<Eigen::MatrixXd> dm = result.second;
		MatrixXcd d1 = pow( 1 / 2 , -1) * dm[0].cast<complex<double>>();
		MatrixXcd d2 = pow( 1 / 2 , -2) * dm[1].cast<complex<double>>();

		double wmin = 1e3 * 2 * M_PI; //fmin 1kHz
		double wmax = 0.1e6 * 2 * M_PI; //fmax = 0.1MHz
		double step = 1e2 * 2 * M_PI; //step = 1kHz
		
		const int steps = (int) ( (wmax - wmin) / step ) + 1;

		//w = wmin: step : wmax; %%CROSSING POINT AT 25
		MatrixXd w = MatrixXd::Zero(1, steps);

		{
			for (int j = 0; j <= (steps - 1); wmin += step, j++)
			{
				w(0, j) = wmin;
			}

		}

		// detail OF CROSSING POINT

		//for m = 400:25 : (size(w, 2) - 5900); %from 5000 onwards something is wrong
		//Set up m FOR loop
		// for m=201:1:310;
		//for m = 25; %%CROSSING POINT AT 25
		int m = 201;

		MatrixXcd c = MatrixXcd::Zero(6, 6);

		complex<double> i(0, 1);

		c(1, 1) = 25.69e9 - i * (w(0, (m- 1)) / (2 * M_PI*2e6))*e(1, 1);
		c(1, 2) = 5.65e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(1, 2);
		c(1, 3) = 9.28e7 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(1, 3);
		c(1, 4) = -8.01e7 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(1, 4);
		c(1, 5) = 17.52e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(1, 5);

		c(2, 2) = 12.11e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(2, 2);
		c(2, 3) = 1.33e7 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(2, 3);
		c(2, 4) = -0.86e7 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(2, 4);
		c(2, 5) = 0.22e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(2, 5);

		c(3, 3) = 4.18e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(3, 3);
		c(3, 4) = 1.31e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(3, 4);
		c(3, 5) = 9.49e7 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(3, 5);

		c(4, 4) = 5.35e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(4, 4);
		c(4, 5) = -7.05e7 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(4, 5);

		c(5, 5) = 28.29e9 - i * (w(0, (m - 1)) / (2 * M_PI*2e6))*e(5, 5);

		complex<double> k1 = c(1, 1) / c(5, 5); //kappa1
		complex<double> k2 = c(1, 2) / c(5, 5); //kappa2
		complex<double> k3 = c(1, 3) / c(5, 5); //kappa3
		complex<double> k4 = c(1, 4) / c(5, 5); //kappa4
		complex<double> k5 = c(1, 5) / c(5, 5); //kappa5

		complex<double> k6 = c(2, 2) / c(5, 5); //kappa6
		complex<double> k7 = c(2, 3) / c(5, 5); //kappa7
		complex<double> k8 = c(2, 4) / c(5, 5); //kappa8
		complex<double> k9 = c(2, 5) / c(5, 5); //kappa9

		complex<double> k10 = c(3, 3) / c(5, 5); //kappa10
		complex<double> k11 = c(3, 4) / c(5, 5); //kappa11
		complex<double> k12 = c(3, 5) / c(5, 5); //kappa12

		complex<double> k13 = c(4, 4) / c(5, 5); //kappa13
		complex<double> k14 = c(4, 5) / c(5, 5); //kappa14

		complex<double> k15 = c(5, 5) / c(5, 5); //kappa15

		complex<double> what = w(0, (m - 1)) * (h / (sqrt(c(5, 5) / rho))); //what
			
		//Q matrices

		MatrixXd Z = MatrixXd::Zero(n, n);

		//2nd degree in K_hat
		
		MatrixXcd Q2(3*n, 3*n);

		Q2 << -MatrixXcd::Identity(n, n)*k13, -MatrixXcd::Identity(n, n)*k11, -MatrixXcd::Identity(n, n)*k8,
			-MatrixXcd::Identity(n, n)*k11, -MatrixXcd::Identity(n, n)*k10, -MatrixXcd::Identity(n, n)*k7,
			-MatrixXcd::Identity(n, n)*k8, -MatrixXcd::Identity(n, n)*k7, -MatrixXcd::Identity(n, n)*k6;

		MatrixXcd S2 = MatrixXcd::Zero(3 * n);

		Q2.row(0) = S2.row(0);
		Q2.row(n - 1) = S2.row(n - 1);
		Q2.row(n) = S2.row(n);
		Q2.row(2 * n - 1) = S2.row(2 * n - 1);
		Q2.row(2 * n) = S2.row(2 * n);
		Q2.row(3 * n - 1) = S2.row(3 * n - 1);

		//1st degree in K_hat

		MatrixXcd Q1(3 * n, 3 * n);
		Q1 << i * (k14 + k14) * d1, i * (k12 + k4) * d1, i * (k9 + k11)*d1,
			i * (k12 + k4)*d1, i * (k3 + k3)*d1, i * (k2 + k10)*d1,
			i * (k9 + k11)*d1, i * (k2 + k10)*d1, i * (k7 + k7)*d1;

		MatrixXcd S1(3 * n, 3 * n);
		S1 << i * k4 * MatrixXcd::Identity(n, n), i * k3*MatrixXcd::Identity(n, n), i * k2*MatrixXcd::Identity(n, n),
			i * k11*MatrixXcd::Identity(n, n), i * k10*MatrixXcd::Identity(n, n), i * k7*MatrixXcd::Identity(n, n),
			i * k14*MatrixXcd::Identity(n, n), i * k12*MatrixXcd::Identity(n, n), i * k9*MatrixXcd::Identity(n, n);

		Q1.row(0) = S1.row(0);
		Q1.row(n - 1) = S1.row(n - 1);
		Q1.row(n) = S1.row(n);
		Q1.row(2 * n - 1) = S1.row(2 * n - 1);
		Q1.row(2 * n) = S1.row(2 * n);
		Q1.row(3 * n - 1) = S1.row(3 * n - 1);

		//0th degree in K_hat
		MatrixXcd L0u = k15*d2 + pow(what, 2) * MatrixXcd::Identity(n, n);
		MatrixXcd L0v = k1*d2 + pow(what, 2) * MatrixXcd::Identity(n, n);
		MatrixXcd L0w = k10*d2 + pow(what, 2) * MatrixXcd::Identity(n, n);

		MatrixXcd Q0(3 * n, 3 * n);
		Q0 << L0u, k5*d2, k12*d2,
		k5*d2, L0v, k3*d2,
		k12*d2, k3*d2, L0w;

		MatrixXcd S0(3 * n, 3 * n);
		S0 << k5*d1, k1*d1, k3*d1,
		k12*d1, k3*d1, k10*d1,
		d1, k5*d1, k12*d1;

		Q0.row(0) = S0.row(0);
		Q0.row(n - 1) = S0.row(n - 1);
		Q0.row(n) = S0.row(n);
		Q0.row(2 * n - 1) = S0.row(2 * n - 1);
		Q0.row(2 * n) = S0.row(2 * n);
		Q0.row(3 * n - 1) = S0.row(3 * n - 1);

		//companion matrices
		MatrixXcd M1(6 * n, 6 * n);
		MatrixXcd M2(6 * n, 6 * n);
		M1 << -Q1, -Q0,
			MatrixXcd::Identity(3 * n, 3 * n), MatrixXcd::Zero(3 * n, 3 * n);
		M2 << Q2, MatrixXcd::Zero(3 * n, 3 * n),
			MatrixXcd::Zero(3 * n, 3 * n), MatrixXcd::Identity(3 * n, 3 * n);

		MatrixXcd p, W;
		pair<MatrixXcd, MatrixXcd> eig = UtilityMethods::matlab_eig(M1, M2);
		p = eig.first;
		W = eig.second;

		// khat = diag(k);
		// kh = sort(khat);
		// kh = kh(find(isfinite(real(kh)) == 1));
		// kh = kh(find(real(kh) >= 0));
		// kh = kh(find(imag(kh)>-1e-9));
	*/
	}
}