//============================================================================
// Name        : Sekwencyjny.cpp
// Author      : PM
// Version     :
// Copyright   : Homework and training
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <Eigen>
using namespace std;
using namespace Eigen;

int main() {

	  Vector3d v(1,2,3);
	  Vector3d w(0,1,2);

	  cout << "Dot product: " << v.dot(w) << endl;
	  double dp = v.adjoint()*w; // automatic conversion of the inner product to a scalar
	  cout << "Dot product via a matrix product: " << dp << endl;
	  cout << "Cross product:\n" << v.cross(w) << endl;

	  cout << "Here is the matrix A:\n" << v.adjoint() << endl;

	  std::cout << "\n\n" << std::endl;



	  std::cout << "\n\n" << std::endl;

	   Matrix2f A, b;
	   A << 2, -1, -1, 3;
	   b << 1, 2, 3, 1;
	   cout << "Here is the matrix A:\n" << A << endl;
	   cout << "Here is the right hand side b:\n" << b << endl;
	   Matrix2f x = A.ldlt().solve(b);
	   cout << "The solution is:\n" << x << endl;
	   return 0;
}
