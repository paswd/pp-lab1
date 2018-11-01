#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int main(void) {
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	/*cout << "Hello server!" << endl;
	cout << "Started core count = " << size << endl;
	cout << "My name is core #" << rank << endl << endl;*/

	return 0;
}

