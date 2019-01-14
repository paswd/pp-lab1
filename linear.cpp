#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

const size_t DEFAULT_STEP = 1;

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

	bool isBorder(size_t i, size_t j, size_t k) {
		if (i == 0 || j == 0 || k == 0) {
			return true;
		}
		if (i == SizeX - 1 || j == SizeY - 1 || k == SizeZ - 1) {
			return true;
		}
		return false;
	}

	void getDataFromFile(string filename) {
		ifstream in(filename.c_str());
		in >> SizeX >> SizeY >> SizeZ >>
			IterCnt >> Coeff;
		in.close();
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

				Func[i][j] = new double[SizeZ];
				NewFunc[i][j] = new double[SizeZ];
				Ro[i][j] = new double[SizeZ];

				for (size_t k = 0; k < SizeZ; k++) {
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
	Puasson3dSolver(string filename) {
		getDataFromFile(filename);
		arrInit();
	}
	~Puasson3dSolver() {
		arrClear();
	}

	void solve() {
		size_t untillX = SizeX - 1;
		size_t untillY = SizeY - 1;
		size_t untillZ = SizeZ - 1;

		for (size_t iteration = 0; iteration < IterCnt; iteration++) {
			//cout << "iteration = " << iteration << endl;
			//cout << "IterCnt = " << IterCnt << endl << endl;
			for (size_t i = 1; i < untillX; i++) {
				//cout << "i = " << i << endl;
				for (size_t j = 1; j < untillY; j++) {
					//cout << "j = " << i << endl;
					for (size_t k = 1; k < untillZ; k++) {
						//cout << "k = " << i << endl;
						NewFunc[i][j][k] = getElementNewValue(i, j, k);
						
						double ***tmp = Func;
						Func = NewFunc;
						NewFunc = tmp;
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



int main(void) {
	//int argc, char * argv[]
	Puasson3dSolver solver("input.txt");
	solver.solve();
	solver.writeResultToFile("output.txt");

	return 0;
}
