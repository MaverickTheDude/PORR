//============================================================================
// Name        : Sekwencyjny.cpp
// Author      : PM
// Version     :
// Copyright   : Homework and training
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <fstream>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

IOFormat exportFmt(FullPrecision, 0, " ", "\n", "", "", "", "");

int main() {

	int Nbodies = 2;
	VectorXd q0(Nbodies), p0(Nbodies);
	p0 << 1, 2;
	q0 << M_PI/2.0, M_PI/2.0;
	inputClass input(Nbodies, p0, q0);

	solution sol = RK_solver(input);

	if (sol.stopped) {
		cout << " PROBLEM!" << endl << sol.error_messege << endl;
		return 1;
	}

	data_set datas(Nbodies);
	int b = 1;
	const VectorXd fi = get_abs_angles(q0);
	const VectorXd omega = get_abs_angles(p0);

	datas.set_S(q0, input);
	datas.set_dS(q0, p0, input);
	cout << datas.tab[b].S1c() << endl << datas.tab[b].s1C() << endl;
	done();
	return 0;

	for (int i = 0; i < Nbodies; i++) {
		input.getBody(i).print();
	}


	std::ofstream outFile;
	outFile.open("results.txt");

	if (outFile.fail() ) {
		std::cerr << "nie udalo sie otworzyc pliku.";
		return 2;
	}

	outFile << sol.T.transpose() << endl;
	outFile << sol.pTab.format(exportFmt) << '\n' << sol.qTab.format(exportFmt) << endl;

	cout << get_abs_angles(input.q0) << endl;
	done();
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

/*	Matrix2f A, b;
	A << 2, -1, -1, 3;
	b << 1, 2, 3, 1;
	cout << "Here is the matrix A:\n" << A << endl;
	cout << "Here is the right hand side b:\n" << b << endl;
	Matrix2f x = A.ldlt().solve(b);
	cout << "The solution is:\n" << x << endl;*/
	return 0;
}
