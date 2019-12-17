#include <fstream>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

IOFormat exportFmt(FullPrecision, 0, " ", "\n", "", "", "", "");
double* times = nullptr;
int main() {
	// STREFA WARUNKOW POCZATKOWYCH
	const int Nbodies = 4;
	VectorXd q0(Nbodies), v0(Nbodies);
	q0 << 0.0, -M_PI/4, 0.0, M_PI/4;
	v0 << -1.0, 2.0, 0.0, 0.5;
	inputClass input(Nbodies, q0, v0);

	// OBLICZENIA
	times = new double[input.N];
	solution sol = RK_solver(input);

	double total_time = 0.0;
	for (int i = 0; i < input.N; i++) {
		total_time += times[i];
		std::cout << "time[" << i << "]: " << times[i] << endl;
	}

	std::cout << "\n\nTotal time: " << total_time << std::endl;

	delete [] times;
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
