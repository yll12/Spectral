#include <iostream>

#include "spectral.h" 

using namespace std;

int main()
{
	double n = 10;
	double m = 2;

	cout << "chebdif(" << n << ", " << m  << ")= " <<
		Spectral::SpectralMethods::chebdif(n, m) << endl;
	
	return 0;
}