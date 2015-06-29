#include <cmath>

// read element of matrix
double ReadMat(double*XX, int *pNN, int *pDD, int i, int d){
	double out=0.0;
	out=XX[(d-1)*(pNN[0])+i-1];
	return out;
}

// write element of matrix
void WriteMat(double*XX, int *pNN, int *pDD, int i, int d, double input){
	XX[(d-1)*(pNN[0])+i-1]=input;
}

// get row of matrix
template <typename RealVector, typename RealMatrix>
inline RealVector matrixRow(const RealMatrix& X, const unsigned rowIdx) {
	const unsigned colNum = X.ncol();
	const unsigned rowNum = X.nrow();
	RealVector rowVector(colNum);
	for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
		rowVector[colIdx] = X[rowIdx + colIdx * rowNum];
	}
	return rowVector;
}

// oneKernel
template <typename RealVector1, typename RealVector2, typename RealMatrix>
inline double oneKernel(const RealVector1& point, const RealMatrix& X, const double h, const RealVector2& weight) {
	const unsigned dimension = X.ncol();
	const unsigned sampleNum = X.nrow();
	double sum, tmp;
	double oneKernelValue = 0.0;

	if (weight.size() == 1) {
		for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
			sum = 0.0;
			for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
				tmp = point[dimIdx] - X[sampleIdx + dimIdx * sampleNum];
				sum += tmp * tmp;
			}
			oneKernelValue += exp(-sum / (2 * h * h));
		}
		return (oneKernelValue / sampleNum);

	} else {
		for (unsigned sampleIdx = 0; sampleIdx < sampleNum; ++sampleIdx) {
			sum = 0.0;
			for (unsigned dimIdx = 0; dimIdx < dimension; ++dimIdx) {
				tmp = point[dimIdx] - X[sampleIdx + dimIdx * sampleNum];
				sum += tmp * tmp;
			}
			oneKernelValue += exp(-sum / (2 * h * h)) * weight[sampleIdx];
		}
		return (oneKernelValue / std::accumulate(weight.begin(), weight.end(), 0.0));
	}
}

// print frame of progress
template <typename Print>
inline void printProgressFrame(Print print) {
	print("0   10   20   30   40   50   60   70   80   90   100");
	print("\n");
	print("|----|----|----|----|----|----|----|----|----|----|\n");
	print("*");
}

// print progress amount
template <typename Print>
inline void printProgressAmount(Print print, int& counter, const int totalCount, int& percentageFloor) {
	int progressAmount = std::floor((100 * (++counter) / totalCount - percentageFloor) / 2);
	if (progressAmount > 0) {
		for (int progressIdx = 1; progressIdx <= progressAmount; ++progressIdx) {
			print("*");
			percentageFloor += 2;
		}
	}
}

