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

	//	Wyznacz wspolczynniki dla wszystkich czlonow
	std::vector<ksi_coef> ksi;
	ksi.reserve(input.Nbodies);
	for (int i = 0; i < input.Nbodies; i++) {
		ksi.emplace_back(p, datas, i);
	}

	//	Assembly phase


	double pi = M_PI;
	VectorXd dY(2*n);
	dY << sin(2*pi*t), cos(2*pi/2*t), sin(2*pi*t), cos(2*pi/2*t);

	return dY;
}

ksi_coef::ksi_coef(const VectorXd &p, const data_set &data, int i)
 	 : _i(i) {
	int Nbodies = p.size();
	i11 = data.tab[i].M1().inverse();
	i22 = data.tab[i].M2().inverse();
	i12 = data.tab[i].M1().llt().solve( data.tab[i].S12() );
	i21 = data.tab[i].M2().llt().solve( data.tab[i].S21() );

	Vector3d rhs10, rhs20;
	rhs10 = data.tab[i].H()*p[i];
	rhs20 = data.tab[i].S21()*data.tab[i].H()*p[i];
	if (i < Nbodies-2) {
		rhs10 -= data.tab[i+1].S12()*data.tab[i+1].H()*p[i+1];
		rhs20 -= data.tab[i+1].H()*p[i+1];
	}
	i10 = data.tab[i].M1().llt().solve(rhs10);
	i20 = data.tab[i].M2().llt().solve(rhs20);
}

