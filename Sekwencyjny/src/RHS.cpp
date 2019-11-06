#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;
extern const double PI;

VectorXd RHS(double t, VectorXd Y, inputClass &input) {
	const int n = input.q0.size();
	VectorXd dY(2*n);
	dY << sin(2*PI*t), cos(2*PI/2*t), sin(2*PI*t), cos(2*PI/2*t);

	return dY;
}
