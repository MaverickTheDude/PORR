#include "PORR.h"
using namespace Eigen;

Matrix2d Rot(double fi) {
	Matrix2d R;
	R << cos(fi), -sin(fi), sin(fi), cos(fi);
	return R;
}

VectorXd get_abs_angles(const VectorXd &q) {
	VectorXd angles(q.size());
	double tmpAngle = 0.0;
	for (int i = 0; i < angles.size(); i++) {
		tmpAngle += q(i);
		angles[i] = tmpAngle;
	}
	return angles;
}

data::data() {
	_D << 1, 0, 0, 1, 0, 0;
	_S12 = Matrix3d::Identity();	_dS12 = Matrix3d::Zero();
	_S21 = Matrix3d::Identity();	_dS21 = Matrix3d::Zero();
	_S1c = Matrix3d::Identity();	_dS1c = Matrix3d::Zero();
	_S2c = Matrix3d::Identity();	_dSc1 = Matrix3d::Zero();
	_Sc1 = Matrix3d::Identity();	_dS2c = Matrix3d::Zero();
	_Sc2 = Matrix3d::Identity();	_dSc2 = Matrix3d::Zero();
}

void data::set_S(const int ith, const double &fi, const inputClass &input) {
	const Matrix3d & M = input.getBody(ith).M();
	const Vector2d & s12_loc = input.getBody(ith).s12;
	const Vector2d & sC1_loc = input.getBody(ith).sC1;
	const Vector2d & sC2_loc = input.getBody(ith).sC2;

	_s12 = Rot(fi)*s12_loc;   	_s21 = -Rot(fi)*s12_loc;
    _sC1 = Rot(fi)*sC1_loc;   	_s1C = -Rot(fi)*sC1_loc;
    _sC2 = Rot(fi)*sC2_loc;   	_s2C = -Rot(fi)*sC2_loc;

    _S12.block<1,2>(2,0) = (Om*_s12).transpose();
    _Sc1.block<1,2>(2,0) = (Om*_sC1).transpose();
    _Sc2.block<1,2>(2,0) = (Om*_sC2).transpose();
    _S21.block<1,2>(2,0) = (Om*_s21).transpose();
    // Dodatkowo inna skladania:
    _S1c.bottomLeftCorner<1,2>() = (Om*_s1C).transpose();
    _S2c.bottomLeftCorner<1,2>() = (Om*_s2C).transpose();

    _M1 = _S1c * M * _S1c.transpose();
    _M2 = _S2c * M * _S2c.transpose();
}

void data::set_dS(const double &om) {
    _dS12.block<1,2>(2,0) = -om*s12().transpose();
    _dSc1.block<1,2>(2,0) = -om*sC1().transpose();
    _dSc2.block<1,2>(2,0) = -om*sC2().transpose();
    _dS21.block<1,2>(2,0) = -om*s21().transpose();
    _dS1c.block<1,2>(2,0) = -om*s1C().transpose();
    _dS2c.block<1,2>(2,0) = -om*s2C().transpose();
}

data_set::data_set(int e_n) : n(e_n) {
	tab = new data[e_n];
	for (int i = 0; i < e_n; i++) {
		tab[i] = data();
	}
}
data_set::~data_set() {
	delete [] tab;
}

void data_set::set_S(const VectorXd &q, const inputClass &input) {
	const VectorXd fi = get_abs_angles(q);
	for (int i = 0; i < n; i++) {
		tab[i].set_S(i, fi[i], input);
	}
}

void data_set::set_dS(const VectorXd &q, const VectorXd &dq, const inputClass &input) {
	const VectorXd omega = get_abs_angles(dq);
	for (int i = 0; i < n; i++) {
		tab[i].set_dS(omega[i]);
	}
}

MatrixXd getVelocity(const VectorXd &dq, const data_set &data) {
	signed int Nbodies = dq.size();
	MatrixXd V1(3, Nbodies);
	V1.col(0) = data.tab[0].H() * dq(0);
	for (int i = 1; i < Nbodies; i++) {
		V1.col(i) = data.tab[i-1].S12().transpose() * V1.col(i-1) +
					data.tab[i].H() * dq(i);
	}
	return V1;
}

Vector3d P1(const double &p, const Ref<VectorXd> &sig, const data &ith_data) {
	return ith_data.H()*p + ith_data.D()*sig;
}


void ksi_coef::print() const {
	std::cout << "ksi_11 \n" << i11 << std::endl << std::endl;
	std::cout << "ksi_12 \n" << i12 << std::endl << std::endl;
	std::cout << "ksi_21 \n" << i21 << std::endl << std::endl;
	std::cout << "ksi_22 \n" << i22 << std::endl << std::endl;
	std::cout << "ksi_10 \n" << i10 << std::endl << std::endl;
	std::cout << "ksi_20 \n" << i20 << std::endl << std::endl;
}

