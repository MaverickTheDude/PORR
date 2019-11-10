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
	tab = new data*[e_n];
	for (int i = 0; i < e_n; i++) {
		tab[i] = new data();
	}
}
data_set::~data_set() {
	for(int i = 0; i < n; i++) {
		delete tab[i];
	}
	delete [] tab;
}

void data_set::set_S(const VectorXd &q, const inputClass &input) {
	const VectorXd fi = get_abs_angles(q);
	for (int i = 0; i < n; i++) {
		tab[i]->set_S(i, fi[i], input);
	}
}

void data_set::set_dS(const VectorXd &q, const VectorXd &dq, const inputClass &input) {
	const VectorXd omega = get_abs_angles(dq);
	for (int i = 0; i < n; i++) {
		tab[i]->set_dS(omega[i]);
	}
}

/*class datas {
std::vector<data> *tab;
std::vector<data> datas;

public:
datas(int e_n) {
	tab = new std::vector<data>;
	for (int i = 0; i < e_n; i++) {
		data *ith_data = new data();
		tab[i].push_back(ith_data);
	}
}
};*/

void done() {
	// Declaring argument for time()
	 time_t tt;
	 // Applying time()
	 time (&tt);

	 // Using localtime()
	 tm * ti = localtime(&tt);

	 std::cout << "\nDone\n";
	 std::cout << "Compilation time: \t" << __TIME__ << std::endl;
//	 std::cout << "Current Day, Date and Time is = " << asctime(ti) << std::endl;
}
