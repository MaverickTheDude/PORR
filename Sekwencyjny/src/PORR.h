// # question 1: czemu musze deklarowac friendship pomiedzy obiema
//   klasami do zadeklarowania body wewnatrz input?
//	 Czy da sie to zrobic inaczej, niz przez friendship?

#ifndef PORR_H_
#define PORR_H_

#include <iostream>
#include <Eigen/Dense>
using namespace Eigen;
const double PI = 4*atan(1);

class body {
	friend class inputClass;
	const int id;
	const double L, m;
	const Vector2d sC1, sC2, s12;
	const DiagonalMatrix<double, 3> M;

	body(int e_id, double e_L, double e_m);
	void print();
};

class inputClass {
	friend class body;
	const int Nbodies;
	// parametry do ustawienia:
	const double L = 1.0;
	const double m = 0.4;
	const double dt = 0.02;
	const double Tk = 1.0;
	body body1 = body(1, L, m);

public:
	const VectorXd p0;
	const VectorXd q0;
	bool stop = false;

	inputClass(const int e_Nbodies, VectorXd e_p0, VectorXd e_q0);
	double getTk();
	double getdt();
	void print();
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

VectorXd RHS(double t, VectorXd Y, inputClass &input);


#endif /* PORR_H_ */
