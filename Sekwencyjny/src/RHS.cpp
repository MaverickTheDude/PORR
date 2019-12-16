#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

VectorXd RHS(const double t, const VectorXd &Y, const inputClass &input) {
	const VectorXd p = Y.head(input.Nbodies);
	const VectorXd q = Y.tail(input.Nbodies);
	VectorXd dq(input.Nbodies), dp(input.Nbodies);
	Vector3d H = Vector3d(0.0, 0.0, 1.0);

	data_set datas = data_set(input.Nbodies);
	datas.set_S(q, input);

	//	Wspolczynniki calkowego rownania ruchu
	std::vector<ksi_coef> ksi;
	ksi.reserve(input.Nbodies);
	for (int i = 0; i < input.Nbodies; i++) {
		ksi.emplace_back(p, datas, i);
	}

	//	Assembly - disassemlby phase
	MatrixXd Q1(3, input.Nbodies);
	Q1 = set_forces_at_H1(datas, input);

	std::vector<Assembly> base_assembly;
	std::vector<acc_force> base_Qacc;
	base_assembly.reserve(input.Nbodies);
	base_Qacc.reserve(input.Nbodies);
	for (int i = 0; i < input.tiers_info[0]; i++) {
		base_assembly.emplace_back(ksi[i], datas.tab[i].S12());
		base_Qacc.emplace_back(Q1.col(i), datas.tab[i].S12());
	}

	Assembly AssemblyC = Assembly(base_assembly[0], base_assembly[1]);
	acc_force Qacc_C = acc_force(base_Qacc[0], base_Qacc[1]);
	Assembly AssemblyD = Assembly(base_assembly[2], base_assembly[3]);
	acc_force Qacc_D = acc_force(base_Qacc[2], base_Qacc[3]);
	Assembly AssemblyS = Assembly(AssemblyC, AssemblyD);
	acc_force Qacc_S = acc_force(Qacc_C, Qacc_D);

	AssemblyS.connect_base_body();
	Qacc_S.connect_base_body();

	AssemblyS.disassemble();
	Qacc_S.disassemble();

	AssemblyC.disassemble();
	AssemblyD.disassemble();
	Qacc_C.disassemble();
	Qacc_D.disassemble();

	// velocity calculation
	MatrixXd P1art(3,input.Nbodies);
	dq(0) = H.transpose() * base_assembly[0].calculate_V1();
	P1art.col(0) = base_assembly[0].T1() + H*p(0);

	for (int i = 1; i < input.tiers_info[0]; i++) {
		Vector3d V1B = base_assembly[i].calculate_V1();
		Vector3d V2A = base_assembly[i-1].calculate_V2();
		dq(i) = H.transpose() * (V1B - V2A);
		P1art.col(i) = base_assembly[i].T1() + H*p(i);
	}
	datas.set_dS(dq);

	// Momenta calculation
	MatrixXd des(3,input.Nbodies);
	des.rightCols(1) << 0, 0, 0;
	for (int i = input.Nbodies - 2; i >= 0; i--) {
		des.col(i) =  (datas.tab[i].dSc2() - datas.tab[i+1].dSc1()) * P1art.col(i+1)
					  + des.col(i+1);
	}

	for (int i = 0; i < input.Nbodies; i++) {
		dp(i) = H.transpose() * (des.col(i) +
				base_Qacc[i].Q1art() - datas.tab[i].dSc1()*P1art.col(i) );
	}

	VectorXd dY(2*input.Nbodies);
	dY << dp, dq;
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

MatrixXd set_forces_at_H1(const data_set &datas, const inputClass &input) {
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
: ksi(_ksi), S12(_S12), AssA(nullptr), AssB(nullptr) {
	D << 1, 0, 0, 1, 0, 0;
}

// Czemu const ref generuje blad fpermissive przy deklarowaniu wskaznika?
Assembly::Assembly(/*const*/ Assembly &A, /*const*/ Assembly &B)
	: AssA(&A), AssB(&B) {
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
}

// pierwszy arg. by val, bo error (czemu?)
acc_force::acc_force(Vector3d _Q1, const Matrix3d &_S12)
	: Q1(_Q1), S12(_S12), AssA(nullptr), AssB(nullptr) { }

acc_force::acc_force( acc_force &A, /*const*/ acc_force &B)
	: AssA(&A), AssB(&B) {
	S12 = A.S12 * B.S12;
	Q1 = A.Q1 + A.S12 * B.Q1;
}

void Assembly::connect_base_body() {
	Matrix2d c = - D.transpose() * ksi.i11 * D;
	_T1 = D * c.ldlt().solve(D.transpose()) * ksi.i10;
	_T2 << 0.0, 0.0, 0.0;
}

void acc_force::connect_base_body() {
	_Q1art = Q1;
	_Q2art << 0.0, 0.0, 0.0;
}

void Assembly::disassemble() {
	Matrix2d C = -D.transpose() * (AssB->ksi.i11 + AssA->ksi.i22) * D;
	Matrix3d W =  D * C.ldlt().solve(D.transpose());
	Vector2d b =  D.transpose() * (AssB->ksi.i10 - AssA->ksi.i20);
	Vector3d beta = D * C.ldlt().solve(b);
	AssB->_T1 = W * AssB->ksi.i12 * _T2 - W * AssA->ksi.i21 * _T1 + beta;
	AssA->_T2 = (-1) * AssB->_T1;
	AssA->_T1 = _T1;
	AssB->_T2 = _T2;
}

Vector3d Assembly::calculate_V1() {
	Vector3d V1 = ksi.i11*_T1 + ksi.i12*_T2 + ksi.i10;
	return V1;
}

Vector3d Assembly::calculate_V2() {
	Vector3d V2 = ksi.i21*_T1 + ksi.i22*_T2 + ksi.i20;
	return V2;
}

void acc_force::disassemble() {
	AssB->_Q1art =  AssB->Q1 - AssB->S12 * _Q2art;
	AssA->_Q2art = -AssB->_Q1art;
	AssA->_Q1art =  _Q1art;
	AssB->_Q2art =  _Q2art;
}
