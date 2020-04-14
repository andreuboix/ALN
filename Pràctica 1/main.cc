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
	 
	 // SOLUTIONS //
	 vector<double> x(n,0);
	 x = resol(A,b);
	 
	 
	 // RESIDUAL VECTOR //
	 vector<double> r(n,0);
	 r = residual(A,x,b);
	 
	 
	 // WE PREPARE THE OUTPUT.TXT
	 Matrix PA = PAing(A, perm);
	 Matrix LU = multiplicacioMM(L,U);
	 Matrix restat = restaMM(PA, LU);
	 double deta = det(A);
	 ofstream output;
	 output.open("output.txt");
	 output << scientific << setprecision(15);
	 output << "The system dimension is " << n << "\n" << "\n";
	 output << "The determinant of the system is: ";
	 output << deta << "\n" << "\n";
	 output << "The permutation vector is: ( ";
	 for(int i=0; i < n; ++i) output << perm[i] << ' ';
	 output << ")" << "\n" << "\n";
	 
	 /*cout << "The number of permutations is: ";
	 if (lu(COPY,perm) == 1 ) cout << "EVEN." << endl << endl;
	 else cout << "ODD." << endl << endl;
	 * */
	 output << "This is the estimation of the error in the PA = LU with the 1-norm and Infinity-Norm:" << "\n";
	 output << "||PA - LU ||_{1} = " << norma1M(restaMM(PAing(A,perm),multiplicacioMM(L,U))) << "\n";
	 output << "||PA - LU ||_{infty} = " << normainfM(restaMM(PAing(A,perm),multiplicacioMM(L,U))) << "\n";
	 output << endl << "This is the estimation of the error in the solution with the 1,2,Infinity-Norms:" << "\n";
	 output << "||Ax - b||_{1} = " << norma1V(r) << "\n";
	 output << "||Ax - b||_{2} = " << norma2V(r) << "\n";
	 output << "||Ax - b||_{infty} = " << normainfV(r) << "\n"; 
	 output << "\n \n \n" << "And this is the solution vector given: x = ( ";
	 for(int i=0; i < n; ++i){
		 output << x[i] << ' ';
	 }
	 output << ")" << "\n \n";
	 output << "And the L and U matrices are these: " << "\n";
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			output << i << " " << j << " " << L[i][j] << "\n"; 
		 }
	 }
	 output << "\n" << "\n" << "\n";
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			output << i << " " << j << " " << U[i][j] << "\n"; 
		 }
	 }
	 output << "\n" << "\n";
	 output << "And the transposed is the next:" << "\n";
	 
	 Matrix T = transposed(COPY2);
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			 output << i << " " << j << " " << T[i][j] << "\n"; 
		 }
	 }
	 output << "\n" << "\n" << "\n";
	 output << "And, finally, the inverse is this one:" << "\n";
	 vector<int> perm2(n);
	 lu(COPY2, perm2);
	 Matrix I = inverse(COPY2);
	 for(int i=0; i < n; ++i){
		 for(int j=0; j < n; ++j){
			 output << i << " " << j << " " << I[i][j] << "\n"; 
		 }
	 }
	 output.close();
	 return 0;
	 
	 
	 // IGNORE THIS... //
			/*cout << "El nombre de permutacions és: ";
			if(c == 1) cout << "PARELL." << endl;
			else cout << "SENAR." << endl;
			vector<double> pb = B(permut, b);
			vector<double> y = Y(L(A),pb,n);
			vector<double> x = X(U(A),y,n);
			Matrix p = P(permut);
			vector<double> Ax = multiplicacioMV( Copia, x );
			vector<double> R = restaVV(Ax, pb);
			Matrix PA = multiplicacioMM(p, Copia);
			Matrix LU = multiplicacioMM(L(A), U(A));
			Matrix E = restaMM ( PA, LU );
			cout << "Aquesta és la dimensió del sistema : " << n << endl;
			cout << "Aquesta és l'estimació de l'error en la descomposició PA = LU en normes 1 i infinit:" << endl;
			cout << norma1M(E) << " " << normainfM(E) << endl;
			cout << "Aquest és el vector de permutació resultant: ( ";
			for(int i=0; i < n; ++i){
				cout << permut[i] << " ";
			} 
			cout << ")" << endl;
			cout << "Aquest és el determinant de la matriu " << det(A) << endl;
			cout << "I aquestes són les estimacions de Ax - b (normes 1,2 i inf respect.): ";
			cout << norma1V(R) << endl;
			cout << norma2V(R) << endl;
			cout <<  normainfV(R) << endl;
			* */
    // Bonus:  càlcul de la inversa de A i dels seus vaps (+ norma2MM(Matrix& A).
    /*//vector<double> s(Matrix& A, int j);
//int swapfiles (Matrix& A, vector<double>& b, int j, vector<double>& permut, int& permutacions, vector<double>& p);
//int parcialesglaonadisim(Matrix& A, vector<double>& b, int n,vector<double>& permut, int& permutacions, vector<double>& p);
//vector<double> X (const Matrix& U, const vector<double>& Y, int n);
//vector<double> Y (const Matrix& L, const vector<double>& b, int n);
//vector<double> inv (const vector<double>& permut);//
* // g++ main.cc lu.cc resol.cc
// ./a.out M18.DAT output.txt
		// Ha de llegir el sistema lineal i escriure la solució a output.txt.
		// es llegeix des de l'arxiu de comandes.
// main.cc : ha de guardar una còpia de A. Aquí és on es calcula l'Ax-b.
// lu.cc :  retorna la LU amb parcial esglaonat, tornant el vector de 
//          permutacio. 1 si #permut parell. -1 si #permut senar. 0 si <tol.
//          
// resol.cc : retorna la solució del sistema lineal. Error de ||PA-LU||. 
//            Agafa el vector de permutació.

//bool singular(const Matrix& A,const vector<double>& b, int n);
//int lu(const Matrix& A, const vector<double>& b, int n, int& permutacions);
* //Matrix L (const Matrix& A);
//Matrix U (const Matrix& A);
//Matrix P(const vector<double>& perm);
* */
}
