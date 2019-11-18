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

inputClass::inputClass(const int &e_Nbodies, VectorXd e_q0, VectorXd e_v0)
 : Nbodies(e_Nbodies), q0(e_q0), v0(e_v0) {
	//bodies = new body[e_Nbodies];
	if (e_q0.rows() != e_Nbodies || e_v0.rows() != e_Nbodies)
			stop = true;
	for (int i = 0; i < e_Nbodies; i++) {
		const int ind = i + 1;
		bodies.push_back(body(ind, L, m));
	}
	v0_to_p0();
}

void inputClass::v0_to_p0() {
	int N = Nbodies;
	inputClass &input = *this;
	VectorXd p0(N);
	MatrixXd sig0(2, N);

	data_set datas(N);
	datas.set_S(q0, input);
	datas.set_dS(q0, v0, input);
	MatrixXd V1 = getVelocity(v0, datas);

	// Baza indukcyjna dla ostatniego czlonu:
	sig0.rightCols(1) = datas.tab[N-1].D().transpose() * datas.tab[N-1].M1() * V1.rightCols(1);
	p0.tail(1) = 		datas.tab[N-1].H().transpose() * datas.tab[N-1].M1() * V1.rightCols(1);
	// Rekursja do korzenia:
	for (int i = N-2; i >= 0; i--) {
		Vector3d P1Bart = P1(p0[i+1], sig0.col(i+1), datas.tab[i+1]);
		Vector3d P1A = datas.tab[i].M1() * V1.col(i);
		p0[i] =   	  datas.tab[i].H().transpose() * (P1A + datas.tab[i].S12() * P1Bart);
		sig0.col(i) = datas.tab[i].D().transpose() * (P1A + datas.tab[i].S12() * P1Bart);
	}
	this->_p0 = p0;
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
	cout << "Warunki poczatkowe: \n p0 = " << v0.transpose() << "\n q0 = " << q0.transpose() << endl;
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


