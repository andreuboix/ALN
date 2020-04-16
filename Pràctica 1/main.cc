#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
using namespace std;
typedef vector<vector <double> > Matrix;

			// LU PART //
void Gauss (Matrix& A, int j);
int pivotatge(Matrix& A, int j);
int lu(Matrix& A, vector<int>& permut);
		// [......................] //

			// SOLVE PART //
vector<double> resol(const Matrix& A, const vector<double>& b);
Matrix LUing(Matrix& A);
vector<double> residual(const Matrix& A,const vector<double>& x,const vector<double>& b);
Matrix PAing (const Matrix& A,const vector<int>& permut);
void Buid(Matrix& A);

	// ops //
Matrix Id(int n);
Matrix inverse (const Matrix& A);
Matrix symmetric (const Matrix& A);
Matrix transposed(const Matrix& A);
Matrix multiplicacioMM(const Matrix & X, const Matrix & A);
vector<double> multiplicacioMV(const Matrix& X, const vector<double>& Y);
vector<double> restaVV (const vector<double>& X, const vector<double>& Y);
Matrix restaMM (const Matrix& X,const  Matrix& Y);

 // determinant and norms of Matrices and Vectors//
double trace (const Matrix& A);
double det(const Matrix& A);
double norma1M(const Matrix& A);
double normainfM(const Matrix& A);
double norma1V(const vector<double>& A);
double norma2V(const vector<double>& A);
double normainfV(const vector<double>& A);
		// [......................] //


int main() {
				// READING THE INPUT ... //
    int n;
    cin >> n;
    Matrix A(n, vector<double>(n,0));
    vector<double> permut(n);
    for(int i=0; i < n; ++i){
		permut[i] = i;
	}
    int v;
    cin >> v;
    for(int c=0; c < v; ++c){
        int i;
        int j;
        cin >> i >> j >> A[i][j];
    }
    int k;
    cin >> k;
    vector <double> b(k,0);
    for(int p=0; p < k; ++p){
        int i;
        cin >> i >> b[i];
    }
				// THIS IS OUR ORIGINAL A & the PERMUTATION vector, USED TO YIELD THE P MATRIX//
				
    Matrix COPY = A;
    Matrix COPY2 = A; 
    Matrix COPY3 = A;
    Matrix COPY4 = A;
    vector<int> perm(n,0);
    
				// HERE BEGINS THE ACTUAL PROGRAM //
	
	
	 // NOT A SINGULAR MATRIX? C'MON... //
	if (lu(COPY,perm) == 0){
		 cout << "This is a SINGULAR MATRIX, I cannot really help you. " << endl;
		 return 0;
	 }
	 
	 // U & L TRIANGULAR MATRICES //
	 Matrix U = COPY;
	 Matrix L = LUing(COPY);
	 Buid(U);
	 // SOLUTIONS //
	 vector<double> x(n,0);
	 x = resol(A,b);
	 
	 
	 // RESIDUAL VECTOR //
	 vector<double> r(n,0);
	 r = residual(A,x,b);
	 
	 
	 // WE PREPARE THE OUTPUT.TXT
	 Matrix PA = PAing(COPY3, perm);
	 Matrix LU = multiplicacioMM(L,U);
	 Matrix restat = restaMM(PA, LU);
	 double deta = det(A);
	 ofstream output;
	 output.open("output.txt");
	 output << "Out of the MATXX.DAT or MSINGULARXX.DAT... " << "\n";
	 output << "This .txt provides you the next information:" << "\n";
	 output << "               -The dimension." << "\n";
	 output << "               -The determinant." << "\n";
	 output  <<"               -The permutation vector." << "\n";
	 output  <<"               -The parity of the permutations." << "\n";
	 output  << "               -The error estimation in PA = LU : 1 and infinity norms." << "\n";
	 output  <<  "               -The error estimation in Ax = b : 1,2 and infinity norms." << "\n";
	 output   << "               -x, the solution vector." << "\n";
	 output  << "               -The PA matrix." << "\n";
	 output   << "               -The L triangular matrix." << "\n";
	 output  << "               -The U triangular matrix." << "\n";
	 output   << "               -The transposed of the A matrix." << "\n";
	 output  << "               -The inverse of the A matrix." << "\n";
	 output << "               -The symmetric matrix created by the multiplication of Atransposed x A." << "\n" << "\n";
	 output << scientific << setprecision(15);
	 output << "The system dimension is " << n << "\n" << "\n";
	 output << "The determinant of the system is: ";
	 output << deta << "\n" << "\n";
	 output << "The permutation vector is: ( ";
	 for(int i=0; i < n; ++i) output << perm[i] << ' ';
	 output << ")" << "\n" << "\n";
	 
	 output << "The number of permutations is: ";
	 vector<int> perm4 (n,0);
	 if (lu(COPY4,perm4) == 1 ) output << "EVEN." << endl << endl;
	 else output << "ODD." << endl << endl;
	 output << "This is the estimation of the error in the PA = LU with the 1-norm and Infinity-Norm:" << "\n";
	 output << "||PA - LU ||_{1} = " << norma1M(restaMM(PAing(COPY3,perm),multiplicacioMM(L,U))) << "\n";
	 output << "||PA - LU ||_{infty} = " << normainfM(restaMM(PAing(COPY3,perm),multiplicacioMM(L,U))) << "\n";
	 output << endl << "This is the estimation of the error in the solution with the 1,2,Infinity-Norms:" << "\n";
	 output << "||Ax - b||_{1} = " << norma1V(r) << "\n";
	 output << "||Ax - b||_{2} = " << norma2V(r) << "\n";
	 output << "||Ax - b||_{infty} = " << normainfV(r) << "\n"; 
	 output << "\n \n \n" << "And this is the solution vector given: x = ( ";
	 for(int i=0; i < n; ++i){
		 output << x[i] << ' ';
	 }
	 output << ")" << "\n";
	 output << "\n" << "\n";
	 output << "The PA matrix has the next form:" << "\n";
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			 output << setw(20)<< PA[i][j] << " " ;
		 }
		 output << "\n" << "\n";
	 }
	 output << "\n";
	 output << "And the U and L matrices are these: " << "\n";
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			output << setw(20)<< U[i][j] << " " ; 
		 }
		 output << "\n" << "\n";
	 }
	  output << "\n" << "\n" << "\n";
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			output << setw(20)<< L[i][j] << " "; 
		 }
		 output << "\n" << "\n";
	 }
	 output << "\n" << "\n";
	 output << "And the transposed is the next:" << "\n";
	 
	 Matrix T = transposed(COPY2);
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			 output << setw(20) << T[i][j]<< " "; 
		 }
		 output << "\n \n";
	 }
	 output << "\n" << "\n" << "\n";
	 output << "And the inverse is this one:" << "\n";
	 vector<int> perm2(n);
	 lu(COPY2, perm2);
	 Matrix I = inverse(COPY2);
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			 output << setw(20) << I[i][j] << " " ; 
		 }
		 output << "\n \n";
	 }
	 output << "And, finally, the S (always symmetric) matrix is this one:" << "\n";
	 Matrix COPY5 = A;
	 Matrix S = symmetric(A);
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			 output << setw(20) << S[i][j] << " "; 
		 }
		 output << "\n \n";
	 }
	 output.close();
	 return 0;
}
