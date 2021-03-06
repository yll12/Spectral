#include "utility.h"

#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace Eigen;

namespace Utility
{
	/*!
	*  Write Eigen Matrices to files
	*/
	void UtilityMethods::EigenToCSV(MatrixXd& x, string filename)
	{
		const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
		ofstream of(filename);
		if (!of) {
			cout << "Cannot open file.\n";
		}
		cout.setf(ios::fixed);
		of << setprecision(numeric_limits<double>::digits10 + 2) << x.format(CSVFormat) << endl;
		of.close();
	}

	Eigen::MatrixXd UtilityMethods::Matlab_fliplr(const MatrixXd& x) 
	{
		return x.rowwise().reverse().eval();
	}

	Eigen::MatrixXd UtilityMethods::Matlab_flipud(const MatrixXd& x) 
	{
		return x.colwise().reverse().eval();
	}
	

}
