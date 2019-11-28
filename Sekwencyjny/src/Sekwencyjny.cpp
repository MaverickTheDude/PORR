#include <fstream>
#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

IOFormat exportFmt(FullPrecision, 0, " ", "\n", "", "", "", "");

class Assembly {
public:
	ksi_coef ksi;
	Matrix3d S12;
	ksi_coef *ksiA, *ksiB;

	Assembly(const ksi_coef &_ksi, const Matrix3d &_S12)
	: ksi(_ksi), S12(_S12)  {
		ksiA = nullptr;
		ksiB = nullptr;
	}

	Assembly(const Assembly &A, const Assembly &B) {

		//To do: zrealizowac konstruktor defaultowy (albo taki, ktory zadziala) dla
		// 			ksi_coef
		// To do: zrozumiec, czemu nie moge zdefiniowac wskaznikow prosto na obiekty A i B
		Vector3d H = Vector3d(0.0, 0.0, 1.0);
		MatrixXd D = MatrixXd(3,2);
		D << 1, 0, 0, 1, 0, 0;
		S12 = A.S12 * B.S12;
		Matrix2d C  = - D.transpose() * (B.ksi.i11 + A.ksi.i22) * D;
		Vector2d b  =   D.transpose() * (B.ksi.i10 - A.ksi.i20);
		Matrix3d W  =   D * C.llt().solve(D.transpose());
		Vector3d beta = D * C.llt().solve(b);

		ksi.i11 =  A.ksi.i11 + A.ksi.i12 * W * A.ksi.i21;
		ksi.i22 =  B.ksi.i22 + B.ksi.i21 * W * B.ksi.i12;
		ksi.i12 = -A.ksi.i12 * W * B.ksi.i12;
		ksi.i21 = -B.ksi.i21 * W * A.ksi.i21;
		ksi.i10 =  A.ksi.i10 - A.ksi.i12 * beta;
		ksi.i20 =  B.ksi.i20 + B.ksi.i21 * beta;

		ksiA = *A.ksi;
		ksiB = *B.ksi;

//		C = - D.' * (ksiB.i11 + ksiA.i22) * D; % no inverse here (!)
//		W = D * (C \ D.');
//		b = D.' * (ksiB.i10 - ksiA.i20);
//		beta = D * (C \ b);
//		ksi11 = ksiA.i11 + ksiA.i12 * W * ksiA.i21;
//		ksi22 = ksiB.i22 + ksiB.i21 * W * ksiB.i12;
//		ksi12 = -ksiA.i12 * W * ksiB.i12;
//		ksi21 = -ksiB.i21 * W * ksiA.i21;
//		ksi10 = ksiA.i10 - ksiA.i12 * beta;
//		ksi20 = ksiB.i20 + ksiB.i21 * beta;
	}
};



int main() {
	// STREFA WARUNKOW POCZATKOWYCH
	int Nbodies = 4;
	VectorXd q0(Nbodies), v0(Nbodies);
	q0 << 0, -M_PI/4, 0, M_PI/4;
	v0 << -1.0, 2, 0, 0.5;
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

	Assembly AssemblyC = Assembly(base_assembly[0], base_assembly[1]);
	Assembly AssemblyD = Assembly(base_assembly[2], base_assembly[3]);
	Assembly AssemblyS = Assembly(AssemblyC, AssemblyD);


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
