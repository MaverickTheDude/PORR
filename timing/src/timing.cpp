#include <iostream>
#include <string>
#include <cstdlib>		// getenv()
#include <stdlib.h>
#include <time.h>
#include <fstream>

using std::cout;
using std::endl;
const std::string getJobId(int, int);
double timeTask(int, int, int);



int main() {
	int threads[] = {1, 2, 4};      int Nth = 3;
	int bodies[] = {512, 1024, 1500, 2048};  int Nbodies = 4;
	int mean = 4;

	double** tab = new double*[Nth];
	for(int j = 0; j < Nth; j++)
		tab[j] = new double[Nbodies];

	for (int i = 0; i < Nth; i++) {
		for (int j = 0; j < Nbodies; j++)
			tab[i][j] = timeTask(threads[i], bodies[j], mean);
	}

	std::ofstream outFile;
	outFile.open("times.txt");
	if (outFile.fail() ) {
		std::cerr << "nie udalo sie otworzyc pliku.";
		return 2;
	}
	for (int i = 0; i < Nth; i++) {
		for (int j = 0; j < Nbodies; j++)
			outFile << tab[i][j] << '\t';
		outFile << endl;
	}

	outFile.close();
	return 0;
}


double timeTask(int threads, int bodies, int mean) {
	double* tab = new double[mean];
	clock_t start, end;
	double t = 0.0;
	const std::string job = getJobId(threads, bodies);

	for (int i = 0; i < mean; i++) {
		start = clock();
		system(job.c_str());
		end = clock();
		tab[i] = double(end - start) / double(CLOCKS_PER_SEC);
	}

	cout << "czasy dla zadania: Nthreads = " << threads << " Nbodies = " << bodies << endl;
	for (int i = 0; i < mean; i++)
		cout << tab[i] << '\t';
	cout << endl << endl;

	// oblicz œredni¹ wszystkich czasów:
	for (int i = 0; i < mean; i++)
		t += tab[i];
	t = t / (double) mean;

	delete [] tab;
	return t;
}

const std::string getJobId(int threads, int bodies) {
	const std::string path = "D:\\code\\PORR\\Sekwencyjny\\Debug\\Sekwencyjny.exe";
	std::string threads_c = std::to_string(threads);
	std::string bodies_c = std::to_string(bodies);
	return path + " " + threads_c + " " + bodies_c;
}

// Niezwiazane z tematem, ale moze kiedys przydatne
//int main()
//{
//    char *homePath(getenv("PWD"));
//    if (homePath == NULL)
//    {
//        std::cout << "$HOME is not set!" << std::endl;
//    }
//    else
//    {
//        std::cout << "$HOME is set to '" << homePath << "'" << std::endl;
//    }
//
//    return homePath == NULL;
//}
