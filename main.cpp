#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

const size_t DEFAULT_STEP = 1;

typedef enum {
	LAYER
} Tags;

class Puasson3dSolver {
private:
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

	double *SendBuffer = NULL;
	double *GetBuffer = NULL;
	size_t BufferSize;
	MPI_Status Status;

	bool isFirstMachine() {
		return MachineId == 0;
	}
	bool isLastMachine() {
		return MachineId == (MachinesCnt - 1);
	}
	size_t getLayersOnMachine(size_t machineId) {
		size_t res = SizeZ / MachinesCnt;
		if (machineId < SizeZ % MachinesCnt) {
			res++;
		}
		return res;
	}
	size_t getOriginalZPos(size_t posInLayer) {
		size_t prevLayersCnt = 0;
		for (size_t i = 0; i < MachineId; i++) {
			prevLayersCnt += getLayersOnMachine(i);
		}
		size_t res = prevLayersCnt + posInLayer;
		if (!isFirstMachine()) {
			res--;
		}
		return res;
	}

	bool isBorder(size_t i, size_t j, size_t k) {
		size_t origK = getOriginalZPos(k);
		if (i == 0 || j == 0 || origK == 0) {
			return true;
		}
		if (i == SizeX - 1 || j == SizeY - 1 || origK == SizeZ - 1) {
			return true;
		}
		return false;
	}


	void setClusterParams(size_t machinesCnt, size_t machineId) {
		MachinesCnt = machinesCnt;
		MachineId = machineId;
		/*LayersOnMachine = SizeZ / MachinesCnt;

		if (MachineId < SizeZ % MachinesCnt) {
			LayersOnMachine++;
		}*/
		LayersOnMachine = getLayersOnMachine(MachineId);
		if (isFirstMachine() || isLastMachine()) {
			LayersOnMachine++;
		} else {
			LayersOnMachine += 2;
		}
	}

	void getDataFromFile(string filename) {
		ifstream in(filename.c_str());
		in >> SizeX >> SizeY >> SizeZ >>
			IterCnt >> Coeff;
		in.close();
		BufferSize = SizeX * SizeY;
	}
	void copyLayerToBuffer(double *buffer, size_t layerId) {
		for (size_t i = 0; i < SizeX; i++) {
			for (size_t j = 0; j < SizeY; j++) {
				buffer[(i * SizeX) + j] = Func[i][j][layerId];
			}
		}
	}
	void getLayerFromBuffer(double *buffer, size_t layerId) {
		for (size_t i = 0; i < SizeX; i++) {
			for (size_t j = 0; j < SizeY; j++) {
				Func[i][j][layerId] = buffer[(i * SizeX) + j];
			}
		}
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

		SendBuffer = new double[BufferSize];
		GetBuffer = new double[BufferSize];

	}

	void arrClear() {
		for (size_t i = 0; i < SizeX; i++) {
			for (size_t j = 0; j < SizeY; j++) {
				delete [] Func[i][j];
				delete [] NewFunc[i][j];
				delete [] Ro[i][j];
			}
			delete [] Func[i];
			delete [] NewFunc[i];
			delete [] Ro[i];
		}
		delete [] Func;
		delete [] NewFunc;
		delete [] Ro;

		delete [] SendBuffer;
		delete [] GetBuffer;
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

public:
	Puasson3dSolver(string filename, size_t machinesCnt, size_t machineId) {
		getDataFromFile(filename);
		setClusterParams(min(machinesCnt, SizeZ), machineId);
		arrInit();
	}
	~Puasson3dSolver() {
		arrClear();
	}

	void solve() {
		size_t untillX = SizeX - 1;
		size_t untillY = SizeY - 1;
		size_t untillZ = LayersOnMachine - 1;

		for (size_t iteration = 0; iteration < IterCnt; iteration++) {
			//Exchange data with previous layer
			if (MachineId > 0) {
				copyLayerToBuffer(SendBuffer, 1);
				MPI_Sendrecv(SendBuffer, BufferSize, MPI_DOUBLE, MachineId - 1, Tags.LAYER,
					GetBuffer, BufferSize, MPI_DOUBLE, MachineId - 1, Tags.LAYER, MPI_COMM_WORLD, &Status);
				getLayerFromBuffer(GetBuffer, 0);
			}
			if (MachineId < MachinesCnt - 1) {
				copyLayerToBuffer(SendBuffer, MachinesCnt - 2);
				MPI_Sendrecv(SendBuffer, BufferSize, MPI_DOUBLE, MachineId + 1, Tags.LAYER,
					GetBuffer, BufferSize, MPI_DOUBLE, MachineId + 1, Tags.LAYER, MPI_COMM_WORLD, &Status);
				getLayerFromBuffer(GetBuffer, MachinesCnt - 1);
			}
			//Exchange data with next layer
			for (size_t i = 1; i < untillX; i++) {
				for (size_t j = 1; j < untillY; j++) {
					for (size_t k = 1; k < untillZ; k++) {
						NewFunc[i][j][k] = getElementNewValue(i, j, k);
					}
				}
			}
			double ***tmp = Func;
			Func = NewFunc;
			NewFunc = tmp;
		}
	}

	void writeResultToFile(string filename) {
		for (size_t i = 0; i < MachinesCnt; i++) {
			if (i == MachineId) {
				ofstream of(filename.c_str());
				if (MachineId == 0) {
					of << "Solution for " << SizeX << "x" << SizeY << "x" << SizeZ << endl;
					of << "Alpha = " << Coeff << endl;
					of << "Iterations count: " << IterCnt << endl << endl;
				}
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
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	bool isUnclaimed() {
		return MachineId >= MachinesCnt;
	}
};


int main(void) {
	//int argc, char * argv[]
	Puasson3dSolver solver("input.txt");
	if (!isUnclaimed) {
		solver.solve();
		solver.writeResultToFile("output.txt");
	}

	return 0;
}
