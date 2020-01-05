#include <fstream>
#include <stdlib.h>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

IOFormat exportFmt(FullPrecision, 0, " ", "\n", "", "", "", "");
double* times = nullptr;

int main(int argc, char* argv[]) {

	// Ustawienie maks liczby watkow (argument aplikacji)
	int Nthr = atoi(argv[1]);
	omp_set_num_threads(Nthr);

	// STREFA WARUNKOW POCZATKOWYCH
// 	const int Nbodies = 4;
	const int Nbodies = atoi(argv[2]);
	VectorXd q0(Nbodies), v0(Nbodies);
 	if (Nbodies == 4) {
	// Warunki poczatkowe dla czterech czlonow
 		q0 << 0.0, -M_PI/4, 0.0, M_PI/4;
 		v0 << -1.0, 2.0, 0.0, 0.5;
 	}
 	else {
 		q0 = VectorXd::Zero(Nbodies);
		v0 = VectorXd::Zero(Nbodies);
		q0(0) = M_PI/6;
 	}

	// OBLICZENIA
	inputClass input(Nbodies, q0, v0);
	times = new double[input.N];

	double t = omp_get_wtime(); //tic
	solution sol = RK_solver(input);
	t =  omp_get_wtime() - t; //toc
	std::cout << "calkowity czas: " << t << std::endl << std::endl;

/*	double total_time = 0.0;
	for (int i = 0; i < input.N; i++) {
		total_time += times[i];
		std::cout << "time[" << i << "]: " << times[i] << endl;
	std::cout << "\n\nTotal time: " << total_time << std::endl;
	}	*/

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
	//	std::cout << "Max threads =" << omp_get_max_threads() << endl;
	std::cout << "Liczba watkow: " << Nthr << std::endl;
	std::cout << "Liczba czlonow: " << input.Nbodies << std::endl;
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
