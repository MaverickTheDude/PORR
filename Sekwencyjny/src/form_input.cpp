#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;


body::body(int e_id, double e_L, double e_m)
	: id(e_id), L(e_L), m(e_m), sC1(-L/2,0), sC2(L/2,0), s12(L,0) {
	DiagonalMatrix<double, 3> M_tmp(m, m, 1.0/12.0*m*L*L);
	_M = M_tmp;
}

void body::print() const {
	cout << "czlon " << id << ": " << endl;
	cout << "dlugosc: \t" << L << endl;
	cout << "massa: \t" << m << endl;
	cout << "interfejsy: \n" << "sC1: " << sC1.transpose() << ",\tsC2: " <<
			 sC2.transpose() << ",\ts12: " << s12.transpose() << endl;
	cout << "Macierz masowa w CM: \t: " << endl << _M << endl;
}

inputClass::inputClass(const int &e_Nbodies, VectorXd e_p0, VectorXd e_q0)
 : Nbodies(e_Nbodies), p0(e_p0), q0(e_q0) {
	//bodies = new body[e_Nbodies];
	if (e_q0.rows() != e_Nbodies || e_p0.rows() != e_Nbodies)
			stop = true;
	for (int i = 0; i < e_Nbodies; i++) {
		const int ind = i + 1;
		bodies.push_back(body(ind, L, m));
	}
}

body inputClass::getBody(int ith) const {
	return bodies[ith];
}

inputClass::~inputClass() {
	//delete [] bodies;
}

double inputClass::getTk() const {
	return Tk;
}

double inputClass::getdt() const {
	return dt;
}

void inputClass::print() {
	cout << "Liczba czlonow: \t" << Nbodies << endl;
	cout << "Dziedzina czasowa: \t [0, " << dt << ", " << Tk << ", ]" << endl;
	cout << "Warunki poczatkowe: \n p0 = " << p0.transpose() << "\n q0 = " << q0.transpose() << endl;
}

solution::solution(VectorXd &e_T, MatrixXd &e_pTab, MatrixXd &e_qTab)
	: T(e_T), pTab(e_pTab), qTab(e_qTab) {
	error_messege = "-- nothing to report -- \n";
	stopped = false;
}

solution::solution(std::string e_error_msg)
	: error_messege(e_error_msg) {
	stopped = true;
}


