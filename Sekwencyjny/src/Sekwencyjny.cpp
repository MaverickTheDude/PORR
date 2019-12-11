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
	Vector3d H = Vector3d(0.0, 0.0, 1.0);
	MatrixXd D = MatrixXd(3,2);
	D << 1, 0, 0, 1, 0, 0;

	// STREFA TESTOW
	data_set datas = data_set(Nbodies);
	datas.set_S(q0, input);

	std::vector<ksi_coef> ksi;
	ksi.reserve(Nbodies);
	for (int i = 0; i < input.Nbodies; i++) {
		ksi.emplace_back(input.p0(), datas, i);
	}

	Matrix<double, 3, Dynamic> Q1;
	Q1 = set_forces_at_H1(datas, input);
	std::cout << Q1 << endl;

	std::vector<Assembly> base_assembly;
	std::vector<acc_force> base_Qacc;
	base_assembly.reserve(input.Nbodies);
	base_Qacc.reserve(input.Nbodies);
	for (int i = 0; i < input.tiers_info[0]; i++) {
		base_assembly.emplace_back(ksi[i], datas.tab[i].S12());
		base_Qacc.emplace_back(Q1.col(i), datas.tab[i].S12());
	}

	Assembly AssemblyC = Assembly(base_assembly[0], base_assembly[1]);
	acc_force Qacc_C = acc_force(base_Qacc[0], base_Qacc[1]);
	Assembly AssemblyD = Assembly(base_assembly[2], base_assembly[3]);
	acc_force Qacc_D = acc_force(base_Qacc[2], base_Qacc[3]);
	Assembly AssemblyS = Assembly(AssemblyC, AssemblyD);
	acc_force Qacc_S = acc_force(Qacc_C, Qacc_D);

	//base body connection
	Matrix2d c = - D.transpose() * AssemblyS.ksi.i11 * D;
	Vector3d T1S = D * c.ldlt().solve(D.transpose()) * AssemblyS.ksi.i10;

	std::cout << Qacc_S.Q1 << std::endl << Qacc_S.S12 << std::endl;
	std::cout << "-----------" << endl;
	std::cout << c << std::endl << T1S << std::endl;

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
