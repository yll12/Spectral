#define _CRTDBG_MAP_ALLOC
#include <Utils.h>

int main()
{

	// Check Memory leak

	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	Eigen::MatrixXd c = Eigen::MatrixXd::Zero(6, 6);
	c << 235.54e9, 88.64e9, 88.78e9, 19.05e9, -2.97e9, -6.54e9,
		0, 215.85e9, 108.43e9, -14.14e9, 2.35e9, 22.16e9,
		0, 0, 215.77e9, -4.90e9, 0.61e9, -15.57e9,
		0, 0, 0, 58.24e9, -17.68e9, -7.02e9,
		0, 0, 0, 0, 39.69e9, 11.66e9,
		0, 0, 0, 0, 0, 44.02e9;
	return 0;
}
