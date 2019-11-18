#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

VectorXd RHS(const double t, const VectorXd &Y, const inputClass &input) {
	const int n = Y.size() / 2;
	const VectorXd p = Y.head(n);
	const VectorXd q = Y.tail(n);

	data_set datas = data_set(n);
	datas.set_S(q, input);

	double pi = M_PI;
	VectorXd dY(2*n);
	dY << sin(2*pi*t), cos(2*pi/2*t), sin(2*pi*t), cos(2*pi/2*t);

	return dY;
}

struct ksi {

};
