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

	AssemblyS.connect_base_body();
	Qacc_S.connect_base_body();

	AssemblyS.disassemble();
	Qacc_S.disassemble();

	AssemblyC.disassemble();
	AssemblyD.disassemble();
	Qacc_C.disassemble();
	Qacc_D.disassemble();

/*	std::cout << base_assembly[0].T1 << std::endl << base_assembly[1].T2 << std::endl;
	std::cout << base_assembly[3].T1 << std::endl << base_assembly[2].T2 << std::endl;
	std::cout << Qacc_C.Q1art << std::endl << Qacc_C.Q2art << std::endl;
	std::cout << Qacc_D.Q1art << std::endl << Qacc_D.Q2art << std::endl;
	std::cout << "---------------" << std::endl;
	std::cout << std::endl << base_Qacc[0].Q1art << std::endl;
	std::cout << std::endl << base_Qacc[1].Q1art << std::endl;
	std::cout << std::endl << base_Qacc[2].Q1art << std::endl;
	std::cout << std::endl << base_Qacc[3].Q1art << std::endl;*/

	MatrixXd V1(3,input.Nbodies), V2(3,input.Nbodies);
	for (int i = 0; i < input.tiers_info[0]; i++) {
		V1.col(i) = base_assembly[i].calculate_V1();
		V2.col(i) = base_assembly[i].calculate_V2();
	}
	Vector4d dq;
	dq(0) = H.transpose() * V1.col(0);
	for (int i = 1; i < input.tiers_info[0]; i++)
		dq(i) = H.transpose() * (V1.col(i) - V2.col(i-1));

	datas.set_dS(dq);

	MatrixXd P1art(3,input.Nbodies);
	VectorXd p = input.p0();
	for (int i = 0; i < input.tiers_info[0]; i++)
		P1art.col(i) = base_assembly[i].T1() + H*p(i);

	for (int i = 0; i < input.tiers_info[0]; i++) {
		std::cout << P1art.col(i) << std::endl << std::endl;
	}

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
