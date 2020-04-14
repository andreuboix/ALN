#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
typedef vector<vector <double> > Matrix;  
int lu(Matrix& A, vector<int>& permut);
				// GIVEN A & b, SOLVES Ax = b & RETURNS x. //

vector<double> resol(const Matrix& A, const vector<double>& b){ // 
	int n = A.size();
	vector<double> x(n,0);
	vector<int> perm (n,0);
	Matrix LU = A;
	lu(LU,perm);
	
	// B = Pb //
	vector<double> B(n,0);
	for(int i=0; i < n; ++i){
		B[i] = b[perm[i]];
	}
	// Ly = B //
	for(int i=0; i < n; ++i){
		x[i] = B[i];
		for(int j=0; j < i; ++j){
			x[i] -= x[j]*LU[i][j];
		}
	}
	// Ux = y // 
	for(int i= n-1; i >= 0; --i){
		for(int j = n-1; j > i; --j){
			x[i] -= x[j]*LU[i][j];
		}
		x[i] /= LU[i][i];
	}
	return x;
}


		// This procedure gives us the MULTIPLICATORS MATRIX, DOOLITTLE METHOD //
Matrix LUing(Matrix& A){
	int n = A.size();
	Matrix L(n, vector<double> (n,0.0));
	
	for(int i=0; i < n; ++i){
		for(int j=0; j < n; ++j){
			if(j < i) {
				L[i][j] = A[i][j];
				A[i][j] = 0.0;
			}
			if( i == j) L[i][j] = 1.0;
		}
	}
	return L;
}

	// This procedure gives us the RESIDUAL VECTOR //
vector<double> residual(Matrix& A, vector<double>& x, vector<double>& b){
	int n = A.size();
	
	vector<double> r(n);
	
	for(int i=0; i < n; ++i){
		r[i] = -b[i];
		for(int j=0; j < n; ++j){
			r[i] += A[i][j]*x[j];
		}
	}
	return r;
}

	// This procedure gives us the PA matrix //
Matrix PAing (Matrix& A, vector<int>& permut){
	int n = A.size();
	int m = A[0].size();

	Matrix aux(n, vector<double>(m));
	for(int i=0; i < n; ++i){
		aux[i] = A[permut[i]];
	}
	return aux;
}

		// OPERATIONS WITH MATRICES & VECTORS - NORMS & DETERMINANT //
			
			
			
			// WORK IN PROGRESS ........// 
			
			// INVERSE OF A MATRIX //
				// (Observation : Ax_{j} = e_{j}, then LUx_{j} = e_{j}, for all j = 1 _ n.) //
/*Matrix inverse (Matrix& A){
	int n = A.size();
	
	
	
}

			
			// TRANSPOSED OF A MATRIX // 
Matrix transposed (Matrix& A){
	int n = A.size();
	Matrix aux(n);
	for(int i=0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			aux[j][i] = A[i][j];
		}
	}
	return aux;
}			
			// CALCULATING THE S = A^{T}A //
Matrix symmetric (Matrix& A){
	Matrix aux = multiplicacioMM(A, transposed(A));
	return aux;
}

			// CALCULATE THE EIGENVALUES OF A MATRIX // 
							// ... //
							
							
			// CALCULATE THE MATRICIAL 2-NORM //
				// (Observation: It is the maximum in absolute value of the eigenvalues.) //
			*/
			
			// MULTIPLICATION MATRIX - MATRIX //
			
Matrix multiplicacioMM(const Matrix & X, const Matrix & A){
	int n = X.size();
	int m = A.size();
	Matrix mul(n, vector<double>(m,0));
	
	for(int i=0;i<n;i++){    
		for(int j=0;j<m;j++){    
			for(int k=0;k<m;k++){    
				mul[i][j]+=X[i][k]*A[k][j];    
			}    
		}    
	}    	
	return mul;
}

			// MULTIPLICATION MATRIX - VECTOR //
			
vector<double> multiplicacioMV(const Matrix & X, const vector<double> & Y){ // Ep, X ha de ser U.
	int n = X.size();
	int m = Y.size();
	vector<double> mul(n);
	for(int i=0;i<n;i++){     
		for(int k=0;k<m;k++){
			mul[i] += X[i][k]*Y[k];    
		}    
	}        	
	return mul;
}
				
				// SUBSTRACTING VECTOR - VECTOR //
				
vector<double> restaVV(const vector<double>& X, const vector<double>& Y){
	int n = X.size(); // Observe how X.size() == Y.size() in order to be substracted
	vector<double> aux(n,0);
	for(int i=0; i< n; ++i){
		aux[i] = X[i] - Y[i];
	}
	return aux;
}
	

				// SUBSTRACTING MATRIX - MATRIX //

Matrix restaMM (const Matrix& X, const Matrix& Y){
	int n = Y[0].size();
	int m = Y.size();
	Matrix aux = X;
	for(int i=0; i< m; ++i){
		for(int j=0; j < n; ++j){
			aux[i][j] -= Y[i][j];
		}
	}
	return aux;
}


				// DETERMINANT of a TRIANGULAR MATRIX is the TRACE //
				
double det (const Matrix& A){
	int n = A.size();
	double suma = 0;
	for(int k=0; k < n; ++k){
		suma += A[k][k];
	}
	return suma;
}

			// Norm_{1} of a MATRIX//
			
double norma1M(const Matrix& A){
	int n = A.size();
	int m = A[0].size();
	double resul = 0.0;
	
	for(int j=0; j < n; ++j){
		double act = 0.0;
		for(int i=0; i < m; ++i){
			act += abs(A[i][j]);
		}
		resul = max(resul, act);
	}
	return resul;
}


			// Norm_{/infty} of a MATRIX//
			
double normainfM(const Matrix& A){
	int n = A.size();
	int m = A[0].size();
	double resul = 0;
	for(int i=0; i < n; ++i){
		double act = 0;
		for(int j=0; j < m; ++j){
			act += abs(A[i][j]);
			resul = max(resul,act);
		}
	}
	return resul;
}
			// Norm_{1} of a VECTOR //
			
double norma1V(const vector<double>& A){
	double suma = 0;
	int n= A.size();
	for(int i=0; i < n; ++i){	
		suma += abs(A[i]);
	}
	return suma;
}

			// Norm_{2} of a VECTOR //
			
double norma2V(const vector<double>& A){
	int n = A.size();
	double suma = 0;
	for(int i=0; i < n; ++i){
		suma += A[i]*A[i];
	}
	return sqrt(suma);
}

			// Norm_{/infty} of a VECTOR//
			
double normainfV(const vector<double>& A){
	double max = 0;
	double aux = 0;
	bool first = true;
	int n = A.size();
	for(int i =0; i < n; ++i){
		if (first){
			max = abs(A[i]);
			first = false;
		}
		else{
			aux = abs(A[i]);
			if(max < aux) max = aux;
		}
	}
	return max;
}
