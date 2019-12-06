#include <fstream>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

IOFormat exportFmt(FullPrecision, 0, " ", "\n", "", "", "", "");

int main() {
	// STREFA WARUNKOW POCZATKOWYCH
	const int Nbodies = 4;
	VectorXd q0(Nbodies), v0(Nbodies);
	q0 << 0.0, -M_PI/4, 0.0, M_PI/4;
	v0 << -1.0, 2.0, 0.0, 0.5;
	inputClass input(Nbodies, q0, v0);

	// STREFA TESTOW
	data_set datas = data_set(Nbodies);
	datas.set_S(q0, input);

	std::vector<ksi_coef> ksi;
	ksi.reserve(Nbodies);
	for (int i = 0; i < input.Nbodies; i++) {
		ksi.emplace_back(input.p0(), datas, i);
	}

	std::vector<Assembly> base_assembly;
	base_assembly.reserve(input.Nbodies);
	for (int i = 0; i < input.tiers_info[0]; i++) {
		base_assembly.emplace_back(ksi[i], datas.tab[i].S12());
	}

	// ok
/*	for (int i = 0; i < input.tiers_info[0]; i++) {
		std::cout << "assembly[" << i << "] \n\n";
		base_assembly[i].ksi.print();
	}*/
//	TO DO:
//
//	ksi_coef::print (zrobic ogolna wersje z kopiowaniem obiektu *this (utilities)
//	assembly phase - debug dlaczego ass[2] jest wyjatkowo zle ...
//		konstruktor ksi_coef::ksi_coef(const VectorXd &p, const data_set &data, int i)
//		note: jako jedyne p0(2) ma 0.0
//	>> dissasembly phase

	base_assembly[2].ksi.print();

	return 0;

//	Assembly AssemblyC = Assembly(base_assembly[0], base_assembly[1]);
	Assembly AssemblyD = Assembly(base_assembly[2], base_assembly[3]);
//	Assembly AssemblyS = Assembly(AssemblyC, AssemblyD);

/*	std::cout << "ksiC : " << std::endl;
	AssemblyC.ksi.print();

	std::cout << "ksiD : " << std::endl;
	AssemblyD.ksi.print();

	std::cout << "ksiS : " << std::endl;
	AssemblyS.ksi.print();*/

	done();
	return 0;

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
