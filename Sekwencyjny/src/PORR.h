// # question 1: czemu musze deklarowac friendship pomiedzy obiema
//   klasami do zadeklarowania body wewnatrz input?
//	 Czy da sie to zrobic inaczej, niz przez friendship?
// # Jaka jest roznica miedzy {} a () w liscie inicjalizacyjnej?
// # obiekt klasy data_set lepiej dac wewnatrz RHS(), czy przed petla RK_ode()?
//   Lepiej zrobic dynamiczna tablice dynamicznych obiektow data, czy wektor?
//   Wybor na chwile obecna: lokalnie utworzony obiekt data_set z dynamiczna tablica
//   (bo nie-const rozmiar) statycznie tworzonych obiektow 'data' (dynamiczna opcja
//    z "new" zwraca blad Egiena - czemu?)
// # Czemu do funkcji RHS() nie moge wyslac t i Y "by ref" ?

#ifndef PORR_H_
#define PORR_H_

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <ctime>
using namespace Eigen;

Matrix2d Rot(double fi);
const Matrix3d I = Matrix3d::Identity();
const Matrix2d Om = Rot(M_PI / 2.0);

class body {
	friend class inputClass;
	const int id;
	double L, m;
	Matrix3d _M;

	body(int e_id, double e_L, double e_m);
public:
	const Vector2d sC1, sC2, s12;
	const Matrix3d & M() const {return _M;}
	void print() const;
};

class inputClass {
	friend class body;
	const int Nbodies;
	// parametry do ustawienia:
	const double L = 1.0;
	const double m = 0.4;
	const double dt = 0.02;
	const double Tk = 1.0;
	//body body1 = body(1, L, m);
	//body *bodies;
	std::vector<body> bodies;

public:
	const int N = static_cast<int>( std::round(Tk / dt) ) + 1;
	const VectorXd p0;
	const VectorXd q0;
	bool stop = false;

	inputClass(const int &e_Nbodies, VectorXd e_p0, VectorXd e_q0);
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
	data();

	const Matrix3d & M1() const {return _M1;}
	const Matrix3d & M2() const {return _M2;}
	const Matrix3d & S12() const {return _S12;}
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
	void set_dS(const VectorXd &q, const VectorXd &dq, const inputClass &input);
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
VectorXd RHS(const double t, VectorXd Y, inputClass &input);
void done();
VectorXd get_abs_angles(const VectorXd &q0);

#endif /* PORR_H_ */
