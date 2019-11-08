#include <cmath>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

VectorXd RHS(double t, VectorXd Y, inputClass &input) {
	const int n = input.q0.size();
	double pi = M_PI;
	VectorXd dY(2*n);
	dY << sin(2*pi*t), cos(2*pi/2*t), sin(2*pi*t), cos(2*pi/2*t);

	return dY;
}
