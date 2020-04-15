#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
using namespace std;
typedef vector<vector<double> > Matrix;

// LU DECOMPOSITION OF A MATRIX AND GAUSS PROCEDURE //


	// The j-ith step of TRIANGULARISING a MATRIX A.
	
void Gauss(Matrix& A, int j){
	int n = A.size();
	for(int i= j+1; i < n; ++i){
		A[i][j] /= A[j][j];
		for(int k=j+1; k < n; ++k){
			A[i][k] -= A[i][j]*A[j][k];
		}
	}
}

	// Scaled Partial Pivoting (SPP).
	
int pivotatge(Matrix& A, int j){
	int n = A.size();
	int resul = -1;
	double maximum = -1.0;
	
	for(int i = j; i < n; ++i){
		double maxrow = abs(A[i][j]);
		for(int k = j+1; k < n; ++k){
			maxrow = max(maxrow, abs(A[i][k]));
		}
		double p = abs(A[i][j]/maxrow);
		if (p > maximum){
			maximum = p;
			resul = i;
		}
	}
	return resul;
}

	// Returns 1 if even # of permutations, -1 if odd and 0 if the matrix is singular.
	
int lu(Matrix& A, vector<int>& permut){
	int n = A.size();
	int permutacions = 0;
	for(int i=0; i < n; ++i) permut[i] = i;
	
	for(int k = 0; k < n-1; ++k){
		int pivot = pivotatge(A, k);
		if (pivot != k){
			swap(A[k], A[pivot]);
			swap(permut[k], permut[pivot]);
			++permutacions;
		}
		if (abs(A[k][k]) < 0.000000000001) return 0;
		Gauss(A,k);
	}
	if (abs(A[n-1][n-1]) < 0.000000000001) return 0;
	return 2*((permutacions+1)%2)-1;
}







