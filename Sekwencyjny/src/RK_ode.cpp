#include "PORR.h"
using namespace Eigen;
extern double* times;

solution RK_solver(inputClass &input) {
	const double dt = input.getdt();
	const double Tk = input.getTk();
	const int N = input.N;
	const unsigned int n = input.q0.size();
	clock_t start, end;

	if (input.stop)
		return solution("niewlasciwa liczba argumentow wejsciowych");

	VectorXd T = VectorXd::LinSpaced(N, 0, Tk);
	VectorXd y_m1(2*n);
	MatrixXd pTab(n, N);	pTab.col(0) = input.p0();
	MatrixXd qTab(n ,N);	qTab.col(0) = input.q0;
	y_m1.head(n) = input.p0();
	y_m1.tail(n) = input.q0;

	int j = 0;
//	double t = omp_get_wtime(); //tic
	for (int i = 1; i < T.size(); i++) {
		start = clock();

		VectorXd k1 = RHS(T[i-1], y_m1, input);
		VectorXd tmp = dt/2.0 * k1;
		VectorXd k2 = RHS(T[i-1] + dt/2.0, y_m1 + dt/2.0*k1, input);
		VectorXd k3 = RHS(T[i-1] + dt/2.0, y_m1 + dt/2.0*k2, input);
		VectorXd k4 = RHS(T[i-1] + dt,   y_m1 + dt*k3, input);

		VectorXd y = y_m1 +  dt/6 * (k1 + 2*k2 + 2*k3 + k4);

		pTab.col(i) = y.head(n);
		qTab.col(i) = y.tail(n);
		y_m1 = y;

		end = clock();
		double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
		times[j] = time_taken;
		j++;
	}
//	t =  omp_get_wtime() - t; //toc
//	std::cout << "calkowity czas: " << t << std::endl << std::endl;

	solution soln(T, pTab, qTab);
	return soln;
};
