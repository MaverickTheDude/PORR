#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

VectorXd RHS(const double t, const VectorXd &Y, const inputClass &input) {
	const unsigned int n = Y.size() / 2;
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

ksi_coef::ksi_coef(const VectorXd &p, const data_set &data, int i) {
	int Nbodies = p.size();
	i11 = data.tab[i].M1().inverse();
	i22 = data.tab[i].M2().inverse();
	i12 = data.tab[i].M1().ldlt().solve( data.tab[i].S12() );
	i21 = data.tab[i].M2().ldlt().solve( data.tab[i].S21() );

	Vector3d rhs10, rhs20;
	rhs10 = data.tab[i].H()*p[i];
	rhs20 = data.tab[i].S21()*data.tab[i].H()*p[i];
	if (i < Nbodies-1) {
		rhs10 -= data.tab[i+1].S12()*data.tab[i+1].H()*p[i+1];
		rhs20 -= data.tab[i+1].H()*p[i+1];
	}
	i10 = data.tab[i].M1().ldlt().solve(rhs10);
	i20 = data.tab[i].M2().ldlt().solve(rhs20);
}

bool ksi_coef::check_if_ok() const {
	Matrix3d tmp = i21 - i12.transpose();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (abs(tmp(i,j)) > 1e-10){
				std::cout << "Sprawdzamy czy ok: \n" << tmp << std::endl; return false; }
		}
	}
	return true;
}

ksi_coef::ksi_coef() {/*Note: trzeba zdefiniowac, mimo ze jest 100% defaultowy*/}

ksi_coef::ksi_coef(const ksi_coef &_ksi) {
	i11 = _ksi.i11;
	i22 = _ksi.i22;
	i12 = _ksi.i12;
	i21 = _ksi.i21;
	i10 = _ksi.i10;
	i20 = _ksi.i20;
}

Matrix<double, 3, Dynamic> set_forces_at_H1(data_set &datas, inputClass &input) {
	const double g = 9.80665;
	MatrixXd Q1(3, input.Nbodies);
	Vector3d F = Vector3d(0.0, 0.0, 0.0);
	for (int i = 0; i < input.Nbodies; i++) {
		double m = input.bodies[i].m();
		F(1) = - m * g;
		Q1.col(i) = datas.tab[i].S1c() * F;
	}
	return Q1;
}

Assembly::Assembly(const ksi_coef &_ksi, const Matrix3d &_S12)
: ksi(_ksi), S12(_S12)  {
	AssA = nullptr;
	AssB = nullptr;
	D << 1, 0, 0, 1, 0, 0;
}

Assembly::Assembly(const Assembly &A, const Assembly &B) {
	D << 1, 0, 0, 1, 0, 0;
	S12 = A.S12 * B.S12;
	Matrix2d C  = - D.transpose() * (B.ksi.i11 + A.ksi.i22) * D;
	Vector2d b  =   D.transpose() * (B.ksi.i10 - A.ksi.i20);
	Matrix3d W  =   D * C.ldlt().solve(D.transpose());
	Vector3d beta = D * C.ldlt().solve(b);

	ksi.i11 =  A.ksi.i11 + A.ksi.i12 * W * A.ksi.i21;
	ksi.i22 =  B.ksi.i22 + B.ksi.i21 * W * B.ksi.i12;
	ksi.i12 = -A.ksi.i12 * W * B.ksi.i12;
	ksi.i21 = -B.ksi.i21 * W * A.ksi.i21;
	ksi.i10 =  A.ksi.i10 - A.ksi.i12 * beta;
	ksi.i20 =  B.ksi.i20 + B.ksi.i21 * beta;

	AssA = &A;
	AssB = &B;
}

// pierwszy arg. by val, bo error (czemu?)
acc_force::acc_force(Vector3d _Q1, const Matrix3d &_S12)
	: Q1(_Q1), S12(_S12) {
	AssA = nullptr;
	AssB = nullptr;
}

acc_force::acc_force(const acc_force &A, const acc_force &B) {
	S12 = A.S12 * B.S12;
	Q1 = A.Q1 + A.S12 * B.Q1;
	AssA = &A;
	AssB = &B;
}

void Assembly::connect_base_body() {
	Matrix2d c = - D.transpose() * ksi.i11 * D;
	T1 = D * c.ldlt().solve(D.transpose()) * ksi.i10;
	T2 << 0.0, 0.0, 0.0;
}

void acc_force::connect_base_body() {
	Q1art = Q1;
	Q2art << 0.0, 0.0, 0.0;
}

void Assembly::disassemble() {
Matrix2d C = -D.transpose() * (AssB->ksi.i11 + AssA->ksi.i22) * D;
Matrix3d W =  D * C.ldlt().solve(D.transpose());
Vector2d b =  D.transpose() * (AssB->ksi.i10 - AssA->ksi.i20);
Vector3d beta = D * C.ldlt().solve(b);
AssB->T1 = W * AssB->ksi.i12 * T2 - W * AssA->ksi.i21 * T1 + beta;
AssA->T2 = (-1) * AssB->T1;
AssA->T1 = T1;
AssB->T2 = T2;
}

void acc_force::disassemble() {

}
