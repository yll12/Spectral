#include <iostream>
#include "Code11.h"
#include "code12.h"
#include "code5Triclinic.h"
#include "utils.h"
#include "ORTHOGONALDISPCYL.h"
#include "spectral.h" 
#include <Eigen/Dense>
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <iomanip>
#include <ctime> 
#include <chrono>
#include <utility>

using namespace std;
using namespace Eigen;
using namespace Spectral;
using namespace Utility;


string numToStr(int n) {
	string str(1, 'a' + n % 26);
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

	/*
	int n = 20;
	int m = 2;
	
	auto start = chrono::system_clock::now();
	DerivativeMatrix result = SpectralMethods::chebdif(n, m);
	auto end = chrono::system_clock::now();

	for (int i = 0; i < m; i++)
	{
		cout << "Here's the matrix dm[" << i << "]: \n" << 
			setprecision(15) << result.dm[i] << "\n";
		cout << "Here's the matrix x:\n" << result.x << "\n";
		ostringstream buffer;
		buffer << "../../D" << (i + 1) << "cpp.csv";
		string filename = buffer.str();
		UtilityMethods::eigenToCSV(result.dm[i], filename);
	}
	UtilityMethods::eigenToCSV(result.x, "../../xcpp.csv");
	auto elapsed = chrono::duration_cast<std::chrono::milliseconds>(end - start);

	cout << "Total time taken is: " << setprecision(15) << elapsed.count() / double(CLOCKS_PER_SEC)
		<< " seconds" << endl;

	/*
	ofstream a_file("../../example.txt"); //this creates it. 
	if (!a_file) {
		cout << "Cannot open file.\n";
		return 1;
	}
	cout << "Reached here" << endl;
	a_file << "Hello World!" << endl;
	a_file.close();
	*/




	/*
	cout << "chebdif(" << n << ", " << m  << ")= " << 
		Spectral::SpectralMethods::chebdif(n, m) << endl;
	*/
	//Spectral::SpectralMethods::chebdif(90, 2);
	//Code11::Code11Method::code11();
	Code12::Code12Method::code12();
	//cout << numToStr(1300) <<endl;
	//Code5Triclinic::Code5TriclinicMethods::code5Triclinic();

	return 0;
}
