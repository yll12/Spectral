#include <iostream>

#include "spectral.h" 

using namespace std;

int main()
{
	double a = 7.4;
	int b = 99;

	cout << "a + b = " <<
		Spectral::SpectralMethods::chebdif(a, b) << endl;
	
	return 0;
}