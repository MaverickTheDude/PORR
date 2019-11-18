#include <fstream>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

IOFormat exportFmt(FullPrecision, 0, " ", "\n", "", "", "", "");

int main() {
	// STREFA WARUNKOW POCZATKOWYCH
	int Nbodies = 2;
	VectorXd q0(Nbodies), v0(Nbodies);
//	q0 << 0, -M_PI/4, 0, M_PI/4;
//	v0 << -1.0, 2, 0, 0.5;
	q0 << 0, -M_PI/4;
	v0 << -1.0, 2;
	inputClass input(Nbodies, q0, v0);

	cout << input.p0() << endl;

	// STREFA TESTOW



	// OBLICZENIA

	solution sol = RK_solver(input);

	if (sol.stopped) {
		cout << " PROBLEM!" << endl << sol.error_messege << endl;
		return 1;
	}

	// DRUKOWANIE WYNIKOW

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

/*	Matrix2f A, b;
	A << 2, -1, -1, 3;
	b << 1, 2, 3, 1;
	cout << "Here is the matrix A:\n" << A << endl;
	cout << "Here is the right hand side b:\n" << b << endl;
	Matrix2f x = A.ldlt().solve(b);
	cout << "The solution is:\n" << x << endl;*/
}


void done() {
/*
	 time_t tt;					// Declaring argument for time()
	 time (&tt);				// Applying time()
	 tm * ti = localtime(&tt);	// Using localtime()
	 std::cout << "Current Day, Date and Time is = " << asctime(ti) << std::endl;
*/
	 std::cout << "\nDone\n";
	 std::cout << "Compilation time: \t" << __TIME__ << std::endl;
}
