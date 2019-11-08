#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;


body::body(int e_id, double e_L, double e_m)
: id(e_id), L(e_L), m(e_m), sC1(-L/2,0),
  sC2(L/2,0), s12(L,0), M(m, m, 1/12*m*L*L)
{}

void body::print() {
	cout << "czlon " << id << ": " << endl;
	cout << "dlugosc: \t" << L << endl;
	cout << "massa: \t" << m << endl;
	cout << "interfejsy: \n" << "sC1: " << sC1.transpose() << ",\tsC2: " <<
			 sC2.transpose() << ",\ts12: " << s12.transpose() << endl;
//	cout << "Macierz masowa w CM: \t: " << M << endl; <- to do
	cout << "Macierz masowa w CM: \t: "  << endl;
}


inputClass::inputClass(const int e_Nbodies, VectorXd e_p0, VectorXd e_q0)
 : Nbodies(e_Nbodies), p0(e_p0), q0(e_q0) {
	if (e_q0.rows() != e_Nbodies || e_p0.rows() != e_Nbodies)
			stop = true;
}

double inputClass::getTk() {
	return Tk;
}

double inputClass::getdt() {
	return dt;
}

void inputClass::print() {
	cout << "Liczba czlonow: \t" << Nbodies << endl;
	body1.print();
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


