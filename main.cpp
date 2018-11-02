#include "mpi.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdint>

using namespace std;

const size_t DEFAULT_STEP = 1;
const size_t INIT_PARAMS_CNT = 5;

typedef enum {
	PARENT_CHILD_INIT,
	CHILD_PARENT_INIT,
	CHILD_CHILD,
	PARENT_CHILD_RESULT,
	CHILD_PARENT_RESULT
} Tags;

class Puasson3dSolver {
private:
	void getDataFromFile(string filename) {
		ifstream in(filename.c_str());
		in >> SizeX >> SizeY >> SizeZ >>
			IterCnt >> Coeff;
		in.close();
	}
protected:
	//Solution variables
	size_t SizeX;
	size_t SizeY;
	size_t SizeZ;

	double StepSqX = 1.;
	double StepSqY = 1.;
	double StepSqZ = 1.;

	size_t IterCnt;
	double Coeff;

	double ***Func;
	double ***NewFunc;
	double ***Ro;

	//Cluster variables
	size_t MachinesCnt;
	size_t MachineId;
	size_t LayersOnMachine;

	char StatusBuffer = 0;
	double InitBuffer[INIT_PARAMS_CNT];
	double *ResultBuffer;

	MPI_Status Status;
	

	virtual bool isBorder(size_t i, size_t j, size_t k) {
		if (i == 0 || j == 0 || k == 0) {
			return true;
		}
		if (i == SizeX - 1 || j == SizeY - 1 || k == SizeZ - 1) {
			return true;
		}
		return false;
	}

	virtual void arrInit() {
		Func = new double**[SizeX];

		for (size_t i = 0; i < SizeX; i++) {
			Func[i] = new double*[SizeY];

			for (size_t j = 0; j < SizeY; j++) {
				Func[i][j] = new double[SizeZ];
			}
		}

	}

	void arrClear() {
		for (size_t i = 0; i < SizeX; i++) {
			for (size_t j = 0; j < SizeY; j++) {
				delete [] Func[i][j];
			}
			delete [] Func[i];
		}
		delete [] Func;
	}

	double getElementNewValue(size_t i, size_t j, size_t k) {
		return (
				(Func[i + 1][j][k] + Func[i - 1][j][k]) / (double) StepSqX +
				(Func[i][j + 1][k] + Func[i][j - 1][k]) / (double) StepSqY +
				(Func[i][j][k + 1] + Func[i][j][k - 1]) / (double) StepSqZ -
				Ro[i][j][k]
			) / (
				2. / (double) StepSqX +
				2. / (double) StepSqY +
				2. / (double) StepSqZ +
				Coeff
			);
	}
	void setClusterParams(size_t machinesCnt, size_t machineId) {
		MachinesCnt = machinesCnt;
		MachineId = machineId;
		LayersOnMachine = SizeZ / machinesCnt;
		if (SizeZ % machinesCnt > 0) {
			LayersOnMachine++;
		}
	}

	void putInitDataIntoBuffer() {
		InitBuffer[0] = (double) SizeX;
		InitBuffer[1] = (double) SizeY;
		InitBuffer[2] = (double) SizeZ;
		InitBuffer[3] = (double) IterCnt;
		InitBuffer[4] = Coeff;
	}
	void getInitDataFromBuffer() {
		SizeX = (size_t) InitBuffer[0];
		SizeY = (size_t) InitBuffer[1];
		SizeZ = (size_t) InitBuffer[2];
		IterCnt = (size_t) InitBuffer[3];
		Coeff = InitBuffer[4];
	}
	size_t getResultBufferSize() {
		return SizeX * SizeY * LayersOnMachine;
	}

	virtual void allocBuffer() {
		/*InitBuffer = malloc(getInitBufferSize());
		LayerBuffer = malloc(getLayerBufferSize());*/
		//InitBuffer = new double[5];
		ResultBuffer = new double[getResultBufferSize()];
		LayerBuffer = new double[SizeX * SizeY];
	}
	virtual void freeBuffer() {
		//free(InitBuffer);
		//free(LayerBuffer);
		delete [] LayerBuffer;
		delete [] ResultBuffer;
	}
private:
	void sendDataToChildren() {
		putInitDataIntoBuffer();
		for (size_t machine = 1; machine < MachinesCnt; machine++) {
			MPI_SENDRECV(InitBuffer, INIT_PARAMS_CNT, MPI_DOUBLE, (int) machine, Tags.PARENT_CHILD_INIT,
				&StatusBuffer, 1, MPI_CHAR, (int) machine, Tags.CHILD_PARENT_INIT, MPI_COMM_WORLD, &Status);
		}
	}

public:
	Puasson3dSolver(string filename, size_t machinesCnt) {
		allocBuffer();
		getDataFromFile(filename);
		setClusterParams(machinesCnt, 0);
		sentDataToChildren();
		arrInit();
	}
	~Puasson3dSolver() {
		arrClear();
		freeBuffer();
	}

	virtual void solve() {
		/*size_t untillX = SizeX - 1;
		size_t untillY = SizeY - 1;
		size_t untillZ = SizeZ - 1;

		for (size_t iteration = 0; iteration < IterCnt; iteration++) {
			for (size_t i = 1; i < untillX; i++) {
				for (size_t j = 1; j < untillY; j++) {
					for (size_t k = 1; k < untillZ; k++) {
						NewFunc[i][j][k] = getElementNewValue(i, j, k);
						
						double ***tmp = Func;
						Func = NewFunc;
						NewFunc = tmp;
					}
				}
			}
		}*/

		for (size_t machine = 1; machine < MachinesCnt; machine++) {
			MPI_SENDRECV(&StatusBuffer, 1, MPI_CHAR, (int) machine, Tags.PARENT_CHILD_RESULT,
				ResultBuffer, getResultBufferSize(), MPI_DOUBLE, (int) machine, Tags.CHILD_PARENT_RESULT,
				MPI_COMM_WORLD, &Status);

			size_t shift = machine - 1;
			for (size_t i = 0; i < SizeX; i++) {
				for (size_t j = 0; j < SizeY; j++) {
					for (size_t k = 0; k < LayersOnMachine; k++) {
						Func[i][j][shift * LayersOnMachine + k] = ResultBuffer[SizeX * SizeY * i + SizeY * j + k];
					}
				}
			}
		}
	}

	void writeResultToFile(string filename) {
		ofstream of(filename.c_str());
		of << "Solution for " << SizeX << "x" << SizeY << "x" << SizeZ << endl;
		of << "Alpha = " << Coeff << endl;
		of << "Iterations count: " << IterCnt << endl << endl;
		for (size_t k = 0; k < SizeZ; k++) {
			of << "k = " << k << endl;

			for (size_t j = 0; j < SizeY; j++) {
				for (size_t i = 0; i < SizeX; i++) {
					if (i > 0) {
						of << " ";
					}
					of << round(Func[i][j][k] * 100.) / 100.;
				}
				of << endl;
			}
			of << endl;
		}
		of.close();
	}
};

class Puasson3dClusterElement : public Puasson3dSolver {
protected:
	double *LayerBuffer = NULL;

	size_t getLayerBufferSize() {
		return sizeof(double) * SizeX * SizeY;
	}

	void allocBuffer() {
		/*InitBuffer = malloc(getInitBufferSize());
		LayerBuffer = malloc(getLayerBufferSize());*/
		//InitBuffer = new double[5];
		ResultBuffer = new double[getResultBufferSize()];
		LayerBuffer = new double[SizeX * SizeY];
	}
	void freeBuffer() {
		//free(InitBuffer);
		//free(LayerBuffer);
		delete [] LayerBuffer;
		delete [] ResultBuffer;
	}

	void getDataFromParent() {
		MPI_SENDRECV(&StatusBuffer, 1, MPI_CHAR, 0, Tags.CHILD_PARENT_INIT,
			InitBuffer, INIT_PARAMS_CNT, MPI_DOUBLE, 0, Tags.PARENT_CHILD_INIT, MPI_COMM_WORLD, &Status);

		getInitDataFromBuffer();
	}

	bool isBorder(size_t i, size_t j, size_t k) {
		if (i == 0 || j == 0 || (k == 0 && MachineId == 1)) {
			return true;
		}
		if (i == SizeX - 1 || j == SizeY - 1 || (k == SizeZ - 1 && MachineId == MachinesCnt - 1)) {
			return true;
		}
		return false;
	}

	void arrInit() {
		Func = new double**[SizeX];
		NewFunc = new double**[SizeX];
		Ro = new double**[SizeX];

		for (size_t i = 0; i < SizeX; i++) {
			Func[i] = new double*[SizeY];
			NewFunc[i] = new double*[SizeY];
			Ro[i] = new double*[SizeY];

			for (size_t j = 0; j < SizeY; j++) {
				Func[i][j] = new double[LayersOnMachine];
				NewFunc[i][j] = new double[LayersOnMachine];
				Ro[i][j] = new double[LayersOnMachine];

				for (size_t k = 0; k < LayersOnMachine; k++) {
					if (isBorder(i, j, k)) {
						Func[i][j][k] = 1.;
						NewFunc[i][j][k] = 1.;
					} else {
						Func[i][j][k] = 0.;
						NewFunc[i][j][k] = 0.;
					}
					Ro[i][j][k] = 0.;
				}
			}
		}
	}
	void putDataToResultBuffer() {
		for (size_t i = 0; i < SizeX; i++) {
			for (size_t j = 0; j < SizeY; j++) {
				for (size_t k = 0; k < LayersOnMachine; k++) {
					ResultBuffer[SizeX * SizeY * i + SizeY * j + k] = Func[i][j][k];
				}
			}
		}
	}
public:
	Puasson3dClusterElement(size_t machinesCnt, size_t machineId) {
		getDataFromParent();
		setClusterParams(machinesCnt, machineId);
		allocBuffer();
		arrInit();
	}
	~Puasson3dClusterElement() {
		arrClear();
		freeBuffer();
	}

	void solve() {
		//Make change between machines
		size_t untillX = SizeX - 1;
		size_t untillY = SizeY - 1;
		size_t untillZ = LayersOnMachine - 1;

		for (size_t iteration = 0; iteration < IterCnt; iteration++) {
			for (size_t i = 1; i < untillX; i++) {
				for (size_t j = 1; j < untillY; j++) {
					for (size_t k = 1; k < untillZ; k++) {
						NewFunc[i][j][k] = getElementNewValue(i, j, k);
						
						double ***tmp = Func;
						Func = NewFunc;
						NewFunc = tmp;
					}
				}
			}
		}
	}

	void returnResult() {
		putDataToResultBuffer();
		MPI_SENDRECV(ResultBuffer, getResultBufferSize, MPI_DOUBLE, 0, Tags.CHILD_PARENT_RESULT,
			&StatusBuffer, 1, MPI_CHAR, 0, Tags.PARENT_CHILD_RESULT, MPI_COMM_WORLD, &Status);
	}
};


int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	/*cout << "Hello server!" << endl;
	cout << "Started core count = " << size << endl;
	cout << "My name is core #" << rank << endl << endl;*/

	if (rank == 0) {
		Puasson3dSolver solver("input.txt", (size_t) size - 1);
		solver.solve();
		solver.writeResultToFile("output.txt");
	} else {
		Puasson3dClusterElement element((size_t) size - 1, (size_t) rank);
		element.solve();
		element.returnResult();
	}

	return 0;
}

