#include <iostream>
#include <Eigen/Dense>
#include "spectral.h" 

using namespace std;
using namespace Eigen;

int main()
{
	int n = 10;
	int m = 2;

	cout << "chebdif(" << n << ", " << m  << ")= " << 
		Spectral::SpectralMethods::chebdif(n, m) << endl;
	
	return 0;
}