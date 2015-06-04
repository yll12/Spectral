#include "TRICLINIC_KV_50.h"

#define _CRTDBG_MAP_ALLOC
#include <string>
#include <Utils.h>

std::string numToStr(int n) {
	std::string str(1, 'a' + n % 26);
	n = n / 26;
	while (n != 0) {
		str = (char)('a' + (n - 1) % 26) + str;
		n = (n - 1) / 26;
	}
	return str;
}

int main()
{

	// Check Memory leak

	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(2, 2);
	c << 1, 2,
		3, 4;
	int G = 2;
	Eigen::MatrixXcd Q12X(G, G);
	Eigen::MatrixXcd Q12Y(G, G);
	Eigen::MatrixXcd Q12Z(G, G);
	Q12X << 1, 2,
		3, 4;
	Q12Y << 5, 6,
		7, 8;
	Q12Z << 9, 10,
		11, 12;
	std::vector<Eigen::MatrixXcd> Q1(6);
	Q1[1] = Eigen::MatrixXcd::Zero(G, 3*G);
	Q1[1] << Q12X, Q12Y, Q12Z;
	std::cout << Q1[1]<< "\n";
	std::vector<Eigen::MatrixXcd> T(6);
	double h = 2;
	double k = 3;
	T[1] = (h * k) * Q1[1];
	std::cout << T[1] << "\n";
	//Code11::Code11Method::code11();
	//TRICLINIC_KV_50::TRICLINIC_KV_50Methods::TRICLINIC_KV_50();
	return 0;
}
