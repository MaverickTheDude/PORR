#include "PORR.h"
using namespace Eigen;
using std::cout;
using std::endl;

VectorXd RHS(const double t, const VectorXd &Y, const inputClass &input) {
	const VectorXd p = Y.head(input.Nbodies);
	const VectorXd q = Y.tail(input.Nbodies);
	VectorXd dq(input.Nbodies), dp(input.Nbodies);
	Vector3d H = Vector3d(0.0, 0.0, 1.0);

	data_set datas = data_set(input.Nbodies);
	datas.set_S(q, input);

	//--- Wspolczynniki calkowego rownania ruchu -----------------------------------------

	//	https://stackoverflow.com/questions/18669296/
	//	c-openmp-parallel-for-loop-alternatives-to-stdvector
	std::vector<ksi_coef> ksi;
	ksi.reserve(input.Nbodies);
	size_t *prefix;

#pragma omp parallel
{
	int ithread  = omp_get_thread_num();
	int nthreads = omp_get_num_threads();
#pragma omp single
	{
		prefix = new size_t[nthreads+1];
		prefix[0] = 0;
	}
	std::vector<ksi_coef> ksi_private;

#pragma omp for schedule(static) nowait
	for (int i = 0; i < input.Nbodies; i++)
		ksi_private.emplace_back(p, datas, i);

	prefix[ithread+1] = ksi_private.size();
	#pragma omp barrier
	#pragma omp single
	{
		for(int i=1; i<(nthreads+1); i++)
			prefix[i] += prefix[i-1];
		//note: tutaj nie potrzebna funkcja resize()
	}
	std::copy(ksi_private.begin(), ksi_private.end(), ksi.begin() + prefix[ithread]);
}
	delete [] prefix;

	//Assembly - disassemlby phase

	MatrixXd Q1(3, input.Nbodies);
	Q1 = set_forces_at_H1(datas, input);

	const unsigned int tiers = input.tiers-1;	// ostatnie pietro definiujemy recznie
	std::vector< std::vector<Assembly> > tree;
	std::vector< std::vector<acc_force> > Q_tree;
	tree.reserve(tiers);
	Q_tree.reserve(tiers);
	tree.resize(tiers);
	Q_tree.resize(tiers);

	//--- Pierwsza galaz drzewa binarnego -----------------------------------------

	std::vector<Assembly> & base_assembly =   tree[0];
	std::vector<acc_force>& base_Qacc     = Q_tree[0];
	base_assembly.reserve(input.Nbodies);
	base_Qacc.reserve(input.Nbodies);

#pragma omp parallel
{
	int ithread  = omp_get_thread_num();
	int nthreads = omp_get_num_threads();
#pragma omp single
	{
		prefix = new size_t[nthreads+1];
		prefix[0] = 0;
	}
	std::vector<Assembly> base_assembly_p;
	std::vector<acc_force> base_Qacc_p;
	unsigned long long int chunk = divide(input.tiers_info[0], nthreads);
	base_assembly_p.reserve(chunk);
	base_Qacc_p.reserve(chunk);

#pragma omp for schedule(static, chunk) nowait
	for (int i = 0; i < input.Nbodies; i++) {
		base_assembly_p.emplace_back(ksi[i], datas.tab[i].S12());
		base_Qacc_p.emplace_back(Q1.col(i), datas.tab[i].S12());
	}
	prefix[ithread+1] = base_assembly_p.size();
	#pragma omp barrier
	#pragma omp single
	{
		for(int i=1; i<(nthreads+1); i++)
			prefix[i] += prefix[i-1];
		// note: funkcja resize() jest niezbedna
		base_assembly.resize(base_assembly.size() + prefix[nthreads]);
		base_Qacc.resize(base_Qacc.size() + prefix[nthreads]);
	}
	std::copy(base_assembly_p.begin(), base_assembly_p.end(), base_assembly.begin() + prefix[ithread]);
	std::copy(base_Qacc_p.begin(),     base_Qacc_p.end(),         base_Qacc.begin() + prefix[ithread]);
}
	delete [] prefix;


	//--- Rekursja od lisci do korzenia -----------------------------------------

	for (unsigned int i = 1; i < tiers; i++) {
		std::vector<Assembly> & branch   =   tree[i];
		std::vector<acc_force>& Q_branch = Q_tree[i];
		branch.reserve(input.tiers_info[i]);
		Q_branch.reserve(input.tiers_info[i]);
		const int end_of_branch = input.tiers_info[i-1] - 1;

		size_t *prefix;
#pragma omp parallel
{
		int ithread  = omp_get_thread_num();
		int nthreads = omp_get_num_threads();
#pragma omp single
		{
			prefix = new size_t[nthreads+1];
			prefix[0] = 0;
		}
		std::vector<Assembly> 	 branch_p;
		std::vector<acc_force> Q_branch_p;
		branch_p.reserve( divide(input.tiers_info[i], nthreads) );
		Q_branch_p.reserve( divide(input.tiers_info[i], nthreads) );

#pragma omp for schedule(static) nowait
		for (int j = 0; j < end_of_branch; j+=2) {
			branch_p.emplace_back(tree[i-1][j], tree[i-1][j+1]);
			Q_branch_p.emplace_back(Q_tree[i-1][j], Q_tree[i-1][j+1]);
		}
	    prefix[ithread+1] = branch_p.size();
	    #pragma omp barrier
	    #pragma omp single
	    {
	        for(int i=1; i<(nthreads+1); i++)
				prefix[i] += prefix[i-1];
	        branch.resize(branch.size() + prefix[nthreads]);
	        Q_branch.resize(Q_branch.size() + prefix[nthreads]);
	    }

	    std::copy(  branch_p.begin(),   branch_p.end(),   branch.begin() + prefix[ithread]);
	    std::copy(Q_branch_p.begin(), Q_branch_p.end(), Q_branch.begin() + prefix[ithread]);
}
		delete [] prefix;

		// Nieparzysta liczba czlonow w galezi powyzej
		if (input.tiers_info[i-1] % 2 == 1) {
			branch.emplace_back(tree[i-1].back());
			Q_branch.emplace_back(Q_tree[i-1].back());
		}
	}

	//--- Base body connection -----------------------------------------

	Assembly AssemblyS = Assembly(tree[tiers-1][0], tree[tiers-1][1]);
	acc_force Q_AssemblyS = acc_force(Q_tree[tiers-1][0], Q_tree[tiers-1][1]);
	AssemblyS.connect_base_body();
	AssemblyS.disassemble();
	Q_AssemblyS.connect_base_body();
	Q_AssemblyS.disassemble();

	//--- Disassembly -----------------------------------------

	for (int i = tiers-1; i > 0; i--) {
		const int end_of_branch = input.tiers_info[i-1] % 2 == 0 ?
				input.tiers_info[i] : input.tiers_info[i]-1;

#pragma omp parallel for schedule(static)
		for (int j = 0; j < end_of_branch; j++) {
			tree[i].at(j).disassemble();
			Q_tree[i][j].disassemble();
		}
		if (input.tiers_info[i-1] % 2 == 1) {
			tree[i-1].back().set(tree[i].back());
			Q_tree[i-1].back().set(Q_tree[i].back());
//			cout << "disassembly: galaz powyzej zawiera nieparzysta liczbe elementow" << endl << endl;
		}
	}

	//--- velocity calculation -----------------------------------------

	MatrixXd P1art(3,input.Nbodies);
	dq(0) = H.transpose() * tree[0][0].calculate_V1();
	P1art.col(0) = tree[0][0].T1() + H*p(0);

#pragma omp parallel for schedule(static)
	for (int i = 1; i < input.tiers_info[0]; i++) {
		Vector3d V1B = tree[0][i].calculate_V1();
		Vector3d V2A = tree[0][i-1].calculate_V2();
		dq(i) = H.transpose() * (V1B - V2A);
		P1art.col(i) = tree[0][i].T1() + H*p(i);
	}
	datas.set_dS(dq);

	//--- descendants term calculation (reverse cummulative sum) -------
	MatrixXd des = MatrixXd::Zero(3, input.Nbodies);
	Vector3d sum;
	MatrixXd sumThr(3, omp_get_max_threads()+1);
	sumThr.leftCols(1).setZero();

	// Oblicz cumsum w podprzedzialach i zapisz ostatnia (najwieksza) wartosc w tablicy sumThr
# pragma omp parallel for schedule(static) private(sum)
	for (int i = input.Nbodies - 2; i >= 0; i--) {
		sum += (datas.tab[i].dSc2() - datas.tab[i+1].dSc1()) * P1art.col(i+1);
		des.col(i) = sum;
		sumThr.col(omp_get_thread_num()+1) = sum;
	}

	// Oblicz poprawke dla poszczegolnych przedzialow. Tylko Nthr (lub Nthr-1) cykli
	for (int i=1; i < omp_get_max_threads(); i++)
		sumThr.col(i) += sumThr.col(i-1);

	// Zastosuj poprawke
# pragma omp parallel for schedule(static)
	for (int i = input.Nbodies - 2; i >= 0; i--) {
		des.col(i) += sumThr.col(omp_get_thread_num());
	}

#pragma omp parallel for schedule(static)
	for (int i = 0; i < input.Nbodies; i++) {
		dp(i) = H.transpose() * (des.col(i) +
				Q_tree[0][i].Q1art() - datas.tab[i].dSc1()*P1art.col(i) );
	}

	VectorXd dY(2*input.Nbodies);
	dY << dp, dq;
	return dY;
}

ksi_coef::ksi_coef(const VectorXd &p, const data_set &data, int i) {
	int Nbodies = p.size();
	i11 = data.tab[i].M1().inverse();
	i22 = data.tab[i].M2().inverse();
	i12 = data.tab[i].M1().ldlt().solve( data.tab[i].S12() );
	i21 = data.tab[i].M2().ldlt().solve( data.tab[i].S21() );

	Vector3d rhs10, rhs20;
	rhs10 = data.tab[i].H()*p[i];
	rhs20 = data.tab[i].S21()*data.tab[i].H()*p[i];
	if (i < Nbodies-1) {
		rhs10 -= data.tab[i+1].S12()*data.tab[i+1].H()*p[i+1];
		rhs20 -= data.tab[i+1].H()*p[i+1];
	}
	i10 = data.tab[i].M1().ldlt().solve(rhs10);
	i20 = data.tab[i].M2().ldlt().solve(rhs20);
}

bool ksi_coef::check_if_ok() const {
	Matrix3d tmp = i21 - i12.transpose();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (abs(tmp(i,j)) > 1e-10){
				std::cout << "Sprawdzamy czy ok: \n" << tmp << std::endl; return false; }
		}
	}
	return true;
}

ksi_coef::ksi_coef() { }

Assembly::Assembly() : AssA(nullptr), AssB(nullptr) { }

acc_force::acc_force() : AssA(nullptr), AssB(nullptr) { }

ksi_coef::ksi_coef(const ksi_coef &_ksi) {
	i11 = _ksi.i11;
	i22 = _ksi.i22;
	i12 = _ksi.i12;
	i21 = _ksi.i21;
	i10 = _ksi.i10;
	i20 = _ksi.i20;
}

MatrixXd set_forces_at_H1(const data_set &datas, const inputClass &input) {
	const double g = 9.80665;
	MatrixXd Q1(3, input.Nbodies);
	Vector3d F = Vector3d(0.0, 0.0, 0.0);

// w tym przypadku dziala taki prosty one-liner. Eigen chyba zawiera wlasna obsluge openMP
#pragma omp parallel for schedule(static)
	for (int i = 0; i < input.Nbodies; i++) {
		double m = input.bodies[i].m();
		F(1) = - m * g;
		Q1.col(i) = datas.tab[i].S1c() * F;
	}

	return Q1;
}

Assembly::Assembly(const ksi_coef &_ksi, const Matrix3d &_S12)
: ksi(_ksi), S12(_S12), AssA(nullptr), AssB(nullptr) {
	D << 1, 0, 0, 1, 0, 0;
}

// const ref generuje blad fpermissive przy deklarowaniu wskaznika (mozliwosc pozniejszej zmiany)
Assembly::Assembly(Assembly &A, Assembly &B)
	: AssA(&A), AssB(&B) {
	D << 1, 0, 0, 1, 0, 0;
	S12 = A.S12 * B.S12;
	Matrix2d C  = - D.transpose() * (B.ksi.i11 + A.ksi.i22) * D;
	Vector2d b  =   D.transpose() * (B.ksi.i10 - A.ksi.i20);
	Matrix3d W  =   D * C.ldlt().solve(D.transpose());
	Vector3d beta = D * C.ldlt().solve(b);

	ksi.i11 =  A.ksi.i11 + A.ksi.i12 * W * A.ksi.i21;
	ksi.i22 =  B.ksi.i22 + B.ksi.i21 * W * B.ksi.i12;
	ksi.i12 = -A.ksi.i12 * W * B.ksi.i12;
	ksi.i21 = -B.ksi.i21 * W * A.ksi.i21;
	ksi.i10 =  A.ksi.i10 - A.ksi.i12 * beta;
	ksi.i20 =  B.ksi.i20 + B.ksi.i21 * beta;
}

// wskazniki przy konstruktorze kopiujacym nie sa wykorzystywane
// Zadanie przejete przez metode set
// Zadeklarowano konstruktor kopiujacy (dokladny, tj. const Assembly&), poniewaz
// stl nie przyjmuje wersji bez const przy zapelnianiu vectora
Assembly::Assembly(const Assembly &A)
	: ksi(A.ksi), S12(A.S12), AssA(A.AssA), AssB(A.AssB), D(A.D) { }

acc_force::acc_force(const acc_force &A)
	: Q1(A.Q1), S12(A.S12), AssA(A.AssA), AssB(A.AssB) { }

// pierwszy arg. by val, bo error (czemu?)
acc_force::acc_force(Vector3d _Q1, const Matrix3d &_S12)
	: Q1(_Q1), S12(_S12), AssA(nullptr), AssB(nullptr) { }

acc_force::acc_force(acc_force &A, acc_force &B)
	: AssA(&A), AssB(&B) {
	S12 = A.S12 * B.S12;
	Q1 = A.Q1 + A.S12 * B.Q1;
}

void Assembly::connect_base_body() {
	Matrix2d c = - D.transpose() * ksi.i11 * D;
	_T1 = D * c.ldlt().solve(D.transpose()) * ksi.i10;
	_T2 << 0.0, 0.0, 0.0;
}

void acc_force::connect_base_body() {
	_Q1art = Q1;
	_Q2art << 0.0, 0.0, 0.0;
}

void Assembly::disassemble() {
	Matrix2d C = -D.transpose() * (AssB->ksi.i11 + AssA->ksi.i22) * D;
	Matrix3d W =  D * C.ldlt().solve(D.transpose());
	Vector2d b =  D.transpose() * (AssB->ksi.i10 - AssA->ksi.i20);
	Vector3d beta = D * C.ldlt().solve(b);

	AssB->_T1 = W * AssB->ksi.i12 * _T2 - W * AssA->ksi.i21 * _T1 + beta;
	AssA->_T2 = (-1) * AssB->_T1;
	AssA->_T1 = _T1;
	AssB->_T2 = _T2;
}

Vector3d Assembly::calculate_V1() {
	Vector3d V1 = ksi.i11*_T1 + ksi.i12*_T2 + ksi.i10;
	return V1;
}

Vector3d Assembly::calculate_V2() {
	Vector3d V2 = ksi.i21*_T1 + ksi.i22*_T2 + ksi.i20;
	return V2;
}

void acc_force::disassemble() {
	AssB->_Q1art =  AssB->Q1 - AssB->S12 * _Q2art;
	AssA->_Q2art = -AssB->_Q1art;
	AssA->_Q1art =  _Q1art;
	AssB->_Q2art =  _Q2art;
}

void Assembly::set(const Assembly &A) {
	_T1 = A._T1;
	_T2 = A._T2;
}

void acc_force::set(const acc_force &A) {
	_Q1art = A._Q1art;
	_Q2art = A._Q2art;
}


// Bez tych operatorow nie dziala std::copy()
Assembly& Assembly::operator=(const Assembly& A) {
	this->ksi = A.ksi;
	this->S12 = A.S12;
	this->_T1 = A._T1;
	this->_T2 = A._T2;
	this->AssA = A.AssA;
	this->AssB = A.AssB;
	this->D = A.D;
	this->H = A.H;

	return *this;
}
acc_force& acc_force::operator=(const acc_force& A) {
	this->Q1 = A.Q1;
	this->S12 = A.S12;
	this->_Q1art = A._Q1art;
	this->_Q2art = A._Q2art;
	this->AssA = A.AssA;
	this->AssB = A.AssB;

	return *this;
}
