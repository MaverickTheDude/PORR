// # question 1: czemu musze deklarowac friendship pomiedzy obiema
//   klasami do zadeklarowania body wewnatrz input?
//	 Czy da sie to zrobic inaczej, niz przez friendship?
// # Jaka jest roznica miedzy {} a () w liscie inicjalizacyjnej?
// # obiekt klasy data_set lepiej dac wewnatrz RHS(), czy przed petla RK_ode()?
//   Lepiej zrobic dynamiczna tablice dynamicznych obiektow data, czy wektor?
//   Wybor na chwile obecna: lokalnie utworzony obiekt data_set z dynamiczna tablica
//   (bo nie-const rozmiar) statycznie tworzonych obiektow 'data' (dynamiczna opcja
//    z "new" zwraca blad Egiena - czemu?)
// # acc_force::pierwszy arg. by val, bo error (czemu?)
// # Assembly::Czemu const ref generuje blad fpermissive przy deklarowaniu wskaznika?

#ifndef PORR_H_
#define PORR_H_

#undef __STRICT_ANSI__
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <ctime>
#include <time.h>
using namespace Eigen;

Matrix2d Rot(double fi);
const Matrix3d I = Matrix3d::Identity();
const Matrix2d Om = Rot(M_PI / 2.0);

class body {
	friend class inputClass;
	const int id;
	const double L, _m;
	Matrix3d _M;

public:
	body(int e_id, double e_L, double e_m);
	const Vector2d sC1, sC2, s12;
	double m() const {return _m;}
	const Matrix3d & M() const {return _M;}
	void print() const;
};

class inputClass {
	friend class body;
	// parametry do ustawienia:
	const double L = 0.4;
	const double m = 0.5;
	const double dt = 0.02;
	const double Tk = 1.0;
	VectorXd _p0;
	void v0_to_p0();

public:
	const int Nbodies;
	const VectorXd q0;
	const VectorXd v0;
	std::vector<body> bodies;
	const unsigned int tiers;
	int* tiers_info;	// pytanie: jak zrobic aby tiers_info bylo const i mialo const elementy?
	const int N = static_cast<int>( round(Tk / dt) ) + 1;
	bool stop = false;
	const VectorXd &p0() const {return _p0;}
	inputClass(const int &e_Nbodies, VectorXd e_q0, VectorXd e_p0);
	~inputClass();
	double getTk() const;
	double getdt() const;
	body getBody(int ith) const;
	void print();
};

class data {
	Matrix3d _M1, _M2, _S12, _S21, _S1c, _S2c, _Sc1, _Sc2,
			 _dS12, _dS21, _dS1c, _dSc1, _dS2c, _dSc2;
	Vector2d _s12, _s21, _s1C, _s2C, _sC1, _sC2;
	Vector3d _H = Vector3d(0.0, 0.0, 1.0);
	MatrixXd _D = MatrixXd(3,2);

public:
	// zgodnie z case 1:
	// http://eigen.tuxfamily.org/dox-devel/group__TopicUnalignedArrayAssert.html
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	data();

	const Matrix3d & M1() const {return _M1;}
	const Matrix3d & M2() const {return _M2;}
	const Matrix3d & S12() const {return _S12;}
	const Matrix3d & S21() const {return _S21;}
	const Matrix3d & S1c() const {return _S1c;}
	const Matrix3d & Sc1() const {return _Sc1;}
	const Vector2d & s12() const {return _s12;}
	const Vector2d & s1C() const {return _s1C;}
	const Vector2d & sC1() const {return _sC1;}

	const Vector2d & s21() const {return _s21;}
	const Vector2d & s2C() const {return _s2C;}
	const Vector2d & sC2() const {return _sC2;}

	const Vector3d & H() const {return _H;}
	const MatrixXd & D() const {return _D;}

	const Matrix3d & dS12() const {return _dS12;}
	const Matrix3d & dS1c() const {return _dS1c;}
	const Matrix3d & dSc1() const {return _dSc1;}
	const Matrix3d & dSc2() const {return _dSc2;}

	void set_S(const int ith, const double &fi, const inputClass &input);
	void set_dS(const double &om);
};

class data_set {
	int n;

public:
	data* tab;
	data_set(int e_n);
	~data_set();
	void set_S(const VectorXd &q, const inputClass &input);
	void set_dS(const VectorXd &dq);
};

struct solution {
	VectorXd T;
	MatrixXd pTab;
	MatrixXd qTab;
	std::string error_messege;
	bool stopped;

	solution(VectorXd &e_T, MatrixXd &e_pTab, MatrixXd &e_qTab);
	solution(std::string e_error_msg);
};

solution RK_solver(inputClass &input);
VectorXd RHS(const double t, const VectorXd &Y, const inputClass &input);
void done();
VectorXd get_abs_angles(const VectorXd &q0);

MatrixXd getVelocity(const VectorXd &dq, const data_set &data);
Vector3d P1(const double &p, const Ref<VectorXd> &sig, const data &ith_data);
VectorXd v0_to_p0(const VectorXd &q0, const VectorXd &v0, const inputClass &input);
MatrixXd set_forces_at_H1(const data_set &datas, const inputClass &input);

struct ksi_coef {
	Matrix3d i11, i12, i21, i22;
	Vector3d i10, i20;
	ksi_coef(const VectorXd &p, const data_set &data, int i);
	ksi_coef(const ksi_coef &_ksi);
	ksi_coef();
	bool check_if_ok() const;
	Matrix3d g11() const {return i11;}
	Matrix3d g12() const {return i12;}
	Matrix3d g21() const {return i21;}
	Matrix3d g22() const {return i22;}
	Vector3d g10() const {return i10;}
	Vector3d g20() const {return i20;}
	void print() const;
};

class Assembly {
public:
	Assembly(const ksi_coef &_ksi, const Matrix3d &_S12);
	Assembly( Assembly &A, /*const*/ Assembly &B);
	void connect_base_body();
	void disassemble();
	Vector3d calculate_V1();
	Vector3d calculate_V2();
	Vector3d T1() const {return _T1;}
	Vector3d T2() const {return _T2;}

private:
	ksi_coef ksi;
	Matrix3d S12;
	Vector3d _T1, _T2;
	Assembly *const AssA, *const AssB;

	// wersja nie-ogolna (rozwiazanie wygodniejsze)
	const Vector3d H = Vector3d(0.0, 0.0, 1.0);
	Matrix<double, 3, 2> D;
};

class acc_force {
public:
	acc_force(Vector3d _Q1, const Matrix3d &_S12); // pierwszy arg. by val, bo error (czemu?)
	acc_force(/*const*/ acc_force &A, /*const*/ acc_force &B);
	void connect_base_body();
	void disassemble();
	Vector3d Q1art() const {return _Q1art;}
	Vector3d Q2art() const {return _Q2art;}

private:
	Vector3d Q1; // accumulated force at H1
	Matrix3d S12;
	Vector3d _Q1art, _Q2art; // acrticulated forces
	acc_force *const AssA, *const AssB;
};

#endif /* PORR_H_ */
