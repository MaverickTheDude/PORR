//============================================================================
// Name        : Sekwencyjny.cpp
// Author      : PM
// Version     :
// Copyright   : Homework and training
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <Eigen/Dense>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

int main() {
	int Nbodies = 1;
	VectorXd q0(Nbodies), p0(Nbodies);
	q0 << 0;
	p0 << 1;
	inputClass input(Nbodies, p0, q0);
	//input.print();
	cout << endl << input.stop << endl;
	solution sol = RK_solver(input);

	if (sol.stopped) {
		cout << " PROBLEM!" << endl << sol.error_messege << endl;
		return 1;
	}

	cout << sol.T.transpose() << endl;
	cout << sol.pTab.transpose() << endl;
	cout << sol.qTab.transpose() << endl;
	return 0;


	//	double m = 0.4, L = 1;
	//	DiagonalMatrix<double, 3> M(m, m, m*L*L/12);
	//	Matrix3d S;
	//	S << 1, 0, 0, 0, 1, 0, m, L, 1;
	//	Matrix3d M1;
	//	M1 = S.transpose() * M * S;
	//	cout << M1 << endl;

	Vector3d v(1,2,3);
	Vector3d w(0,1,2);

	cout << "Dot product: " << v.dot(w) << endl;
	double dp = v.adjoint()*w; // automatic conversion of the inner product to a scalar
	cout << "Dot product via a matrix product: " << dp << endl;
	cout << "Cross product:\n" << v.cross(w) << endl;

	cout << "Here is the matrix A:\n" << v.adjoint() << endl;

	std::cout << "\n\n" << std::endl;



	std::cout << "\n\n" << std::endl;

	Matrix2f A, b;
	A << 2, -1, -1, 3;
	b << 1, 2, 3, 1;
	cout << "Here is the matrix A:\n" << A << endl;
	cout << "Here is the right hand side b:\n" << b << endl;
	Matrix2f x = A.ldlt().solve(b);
	cout << "The solution is:\n" << x << endl;
	return 0;
}
