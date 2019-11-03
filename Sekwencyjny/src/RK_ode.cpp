#include <iostream>
#include <Eigen>
#include <cmath>
#include "PORR.h"
using namespace Eigen;

solution RK_solver(inputClass &input) {
	const double dt = input.getdt();
	const double Tk = input.getTk();
	const int N = static_cast<int>( std::round(Tk/dt) ) + 1;

	if (input.stop)
		return solution("niewlasciwa liczba argumentow wejsciowych");

	// Inicjalizacja dziedziny czasowej
	VectorXd T(N);
	for (int i = 0; i < T.size(); i++) {
		T[i] = static_cast<double>(i) * dt;
	}

	// TO DO

	MatrixXd pTab(2,2);
	MatrixXd qTab(2,2);
	solution soln(T, pTab, qTab);
	return soln;
};
