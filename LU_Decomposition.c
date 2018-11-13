# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <time.h>
# include <math.h>
# include <inttypes.h>

#define RANDMAX 10000000
#define LIN(i,j,m,n) ((i)*(n)+(j))
#define MIN(a,b) ((a)>(b)? (b):(a))

/* Signatures des fonctions de manipulation vecteurs/matrices */

double *allocateVector(unsigned int m);
double* allocateMatrix(unsigned int m, unsigned int n);
void freeVector(double* v);
void freeMatrix(double* A);
void printVector(double *v, unsigned int m);
void printMatrix(double *A, unsigned int m, unsigned int n);
void elementVector(double* v, unsigned int k, unsigned int m);
void scaleVector(double alpha, double* v, unsigned int m);
void setMatrixColumn(double* A, double* v, unsigned int k, unsigned int m, unsigned int n);
void scaleMatrix(double alpha, double* A, unsigned int m, unsigned int n);
void identityMatrix(double *A, unsigned int m, unsigned int n);
void addVector(double* w, double* u, double* v, unsigned int m);
void addMatrix(double* C, double* A, double* B, unsigned int m, unsigned int n);
void copyVector(double* v,double* u, unsigned int m);
void copyMatrix(double* B,double* A, unsigned int m, unsigned int n);
void matrixMatrixProduct(double* C, double* A, double* B, unsigned int m, unsigned int n, unsigned int p);
double randomValue();
void randomMatrix(double *A, unsigned int m, unsigned int n);
void randomVector(double *v, unsigned int m);
void matrixVectorProduct(double* v, double* A, double* u, unsigned int m, unsigned int n);
void gaussianElimination(double *A, double *b, unsigned int n);
void solveTriangularUpper(double* x, double* A, double* b, unsigned int n);


/* Exercice 4 */

void luDecomposition(double *A, unsigned int n){
	unsigned int i,j,k;
	double c;
	for(i=0;i<n;i++){
		for(j=i+1;j<n;j++){
			c = A[LIN(j,i,n,n)]/A[LIN(i,i,n,n)];
			for(k=i+1;k<n;k++){
				A[LIN(j,k,n,n)]-= c*A[LIN(i,k,n,n)];
			}
			
			A[LIN(j,i,n,n)]= c;
		}
	}
}

void solveTriangularLower(double* x, double* A, double* b, unsigned int n){
	// On suppose n déjà triangularisé
	x[0] = b[0];
	unsigned int i,j;
	for(i=1; i<n; i++){
		x[i] = 0;
		for(j=0;j<i; j++){
			x[i]+= x[j]*A[LIN(i,j,n,n)];		
		}
		x[i] = (b[i]-x[i]);
		if(i==0) break;
	}
}

void finalSolveSystemLU(double* x, double* A, double* b, double *scratch, unsigned int n){
	solveTriangularLower(scratch, A, b, n);
	solveTriangularUpper(x, A, scratch, n);
}


void solveSystemLU(double* x, double* A, double* b, unsigned int n){

	double *A1 = allocateMatrix(n, n);
	copyMatrix(A1, A, n, n);
	double *b1 = allocateVector(n);
	copyVector(b1, b, n);

	luDecomposition(A1, n);
	
	double *scratch = allocateVector(n);	
	finalSolveSystemLU(x, A1, b1, scratch, n);

	freeMatrix(A1);
	freeVector(b1);
	freeVector(scratch);
}


void invertMatrix(double* B, double* A, unsigned int n){

	unsigned int i;
	double *A1 = allocateMatrix(n, n);
	copyMatrix(A1, A, n, n);
	double *b1 = allocateVector(n);
	double *b2 = allocateVector(n);
	double *scratch = allocateVector(n);	

	luDecomposition(A1, n);
	for(i=0;i<n;i++){
		scaleVector(0., b2, n);
		elementVector(b1, i, n);
		finalSolveSystemLU(b2, A1, b1, scratch, n);
		setMatrixColumn(B, b2, i, n, n);
	}
	
	freeMatrix(A1);
	freeVector(b1);
	freeVector(b2);
	freeVector(scratch);
}

double matrixDeterminantLU(double *A, unsigned int n){

	unsigned int i;
	double res = 1;
	for(i=0; i<n; i++){
		res *= A[LIN(i,i,n,n)];
	}

	return res;
}

double matrixDeterminant(double *A, unsigned int n){

	double *A1 = allocateMatrix(n, n);
	copyMatrix(A1, A, n, n);
	luDecomposition(A, n);
	
	unsigned int i;
	double res = 1;
	for(i=0; i<n; i++){
		res *= A[LIN(i,i,n,n)];
	}

	return res;
}

int main(){

	double A[] = {2,n,1,5,6,13000,5,19,2,19,10,23000n,4,10,11,30001};
	printMatrix(A,4,4);
	
	printf("\nTest LU décomposition\n");
	luDecomposition(A, 4);
	printMatrix(A,4,4);

	printf("Test solveTriangularLower L\n");
	double L[] = {1,0,0,0,n,1,0,0,1,4,1,0,2,1,7,1};
	double *x = allocateVector(4);
	double b[] = {1,n,n,4};

	solveTriangularLower(x, L, b, 4);
	printVector(x,4);
	printf("Test solveTriangularLower A\n");
	solveTriangularLower(x, A, b, 4);
	printVector(x,4);


	printf("Test finalSolveSystem\n");
	double *scratch = allocateVector(4);		
	finalSolveSystemLU(x, A, b, scratch, 4);
	printVector(x,4);

	printf("Test solveTriangularUpper\n");
	double U[] = {2,n,1,5,0,4,2,4,0,0,1,2,0,0,0,n};
	double *x1 = allocateVector(4);
	double b1[] = {1,0,2,-12};

	solveTriangularUpper(x1, U, b1, 4);
	printVector(x1,4);

	printf("Test solveSystemLU\n");
	
	double C[] = {2,n,1,5,6,1n,5,19,2,19,10,2n,4,10,11,n1};
	double d[] = {1,n,n,4};
	double *y = allocateVector(4);
	solveSystemLU(y, C, d, 4);
	printVector(y,4);

	printf("Test invertMatrix\n");

	double A1[] = {2,n,1,5,6,1n,5,19,2,19,10,2n,4,10,11,n1};
	double *B1 = allocateMatrix(4,4);
	invertMatrix(B1, A1, 4);
	printMatrix(B1,4,4);

	double C1[] = {2,n,1,5,6,1n,5,19,2,19,10,2n,4,10,11,n1};
	printf("Test detMatrix\n");
	printf("determinant: %lf\n", matrixDeterminant(C1, 4));
	
	printf("Test detMatrixLU\n");
	printf("determinant: %lf\n", matrixDeterminantLU(C1, 4));

	
	freeMatrix(B1);
	freeVector(x);
	freeVector(scratch);
	freeVector(x1);
	freeVector(y);



	//TEST Max ||x_tilde - x||
	unsigned int n = 200;
	double *A = allocateMatrix(n,n);
	double *b = allocateVector(n);
	
	double *x = allocateVector(n);

	
	unsigned int i;

	double delta = 0;
	double delta_max = 0;
	double *x_tilde = allocateVector(n);
	double *x_diff = allocateVector(n); //x-x_tilde

	for(i=0; i<100; i++){

		scaleVector(0, x_diff, n);
		scaleVector(0, x_tilde, n);

		randomMatrix(A, n, n);
		randomVector(x, n);
		//printf("Calcul de b = Ax\n");
		matrixVectorProduct(b, A, x, n, n);

		delta = 0;
		unsigned int j;	
		//printVector(x, n);
		//finalSolveSystemLU(x_tilde, A, b, scratch, n);
		solveSystemLU(x_tilde, A, b, n);

		//printVector(x_tilde, n);
		scaleVector(-1, x, n);
		addVector(x_diff, x_tilde, x, n);
		//printVector(x_diff, n);
		scaleVector(-1, x, n);
		for(j=0; j<n; j++) delta += fabs(x_diff[j]);
		delta /= n;
		if(i==0) delta_max = delta;
		else{
			if(delta > delta_max) delta_max = delta;
		}
		//printf("delta_max = %lf\n", delta_max);
	}

	printf("Max ||x_tilde - x|| = %.42lf\n", delta_max);

	//TEST Max ||Id - AA-¹||

	double *Id = allocateMatrix(n,n);
	double *A_1 = allocateMatrix(n,n); //A^(-1)
	double *AA = allocateMatrix(n,n); // AA-¹
	double *AA_diff = allocateMatrix(n,n); // Id-AA-¹
	identityMatrix(Id, n, n);
	delta_max = 0.;
	
	for(i=0; i<100; i++){

		randomMatrix(A, n, n);
		invertMatrix(A_1, A, n);
		matrixMatrixProduct(AA, A, A_1, n, n, n);
	
		delta = 0;
		unsigned int j;

		scaleMatrix(-1, AA, n, n);
		addMatrix(AA_diff, Id, AA, n, n);
		scaleMatrix(-1, AA, n, n);

		for(j=0; j<n*n; j++) delta += fabs(AA_diff[j]);
		delta /= n*n;
		if(i==0) delta_max = delta;
		else{
			if(delta > delta_max) delta_max = delta;
		}
		//printf("delta_max = %lf\n", delta_max);
	}

	printf("Max ||Id - AA-¹|| = %.42lf\n", delta_max);

	return 0;
}

/* Max ||x_tilde - x|| = 0.000000046582972664660492254537054726201695
Max ||Id - AA-¹|| = 0.000000000061101734291687103490660333229194
*/


/* ___________________________________________________________________________________________*/

/* Fonctions manipulation vecteurs/matrices */

double *allocateVector(unsigned int m){
	double* TabV = (double*)calloc(m, sizeof(double));
	if(TabV == NULL) exit(1);
	return TabV;

}

double* allocateMatrix(unsigned int m, unsigned int n){
	double* TabM = (double*)calloc(m*n, sizeof(double));
	if(TabM == NULL) exit(1);
	return TabM;
}

void freeVector(double* v){
	free(v);
}

void freeMatrix(double* A){
	free(A);
}


void printVector(double *v, unsigned int m){
	unsigned int i;
	for(i = 0; i < m; i++) printf("%+1.18e\n",v[i]);
}


void printMatrix(double *A, unsigned int m, unsigned int n){
	int cpt = 0;
	unsigned int i;
	for(i = 0; i < m*n; i++){
		printf("%+1.18e\t",A[i]);
		cpt ++;
		if(cpt%n==0){
			cpt = 0;
			printf("\n");
		}
	}
}

void elementVector(double* v, unsigned int k, unsigned int m){
	unsigned int i;
	for(i = 0; i < m; i++){
		if(i == k) v[i] = 1;
		else v[i] = 0;
	}
}


double randomValue(){
	return 2*((double)(rand()))/((double)RANDMAX)-1;
}

void randomMatrix(double *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i = 0;i < m*n; i++) A[i] = randomValue()/4;
	for(j = 0; j < MIN(m,n); j++){
			A[LIN(j,j,m,n)] = (randomValue()+1)/4 + 0.75;
	}
}

void randomVector(double *v, unsigned int m){
	unsigned int i;
	for(i = 0; i < m; i++) v[i] = randomValue();
}

void setMatrixColumn(double* A, double* v, unsigned int k, unsigned int m, unsigned int n){
	unsigned int p;
	for(p = 0; p < m; p++){
			A[LIN(p,k,m,n)] = v[p];
	}
}

void scaleVector(double alpha, double* v, unsigned int m){
	unsigned int i;
	for(i=0; i<m; i++) v[i] *= alpha;
}


void scaleMatrix(double alpha, double* A, unsigned int m, unsigned int n){
	unsigned int i;
	for(i=0; i<m*n; i++) A[i] *= alpha;
}

void identityMatrix(double *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i = 0;i < m*n; i++) A[i] = 0;
	for(j = 0; j < MIN(m,n); j++){
			A[LIN(j,j,m,n)] = 1;
	}
}

void addMatrix(double* C, double* A, double* B, unsigned int m, unsigned int n){
	unsigned int i;
	for(i=0; i<m*n; i++) C[i] = A[i] + B[i];
}

void addVector(double* w, double* u, double* v, unsigned int m){
	unsigned int i;
	for(i=0; i<m; i++) w[i] = u[i] + v[i];
}

void copyVector(double* v,double* u, unsigned int m){
	unsigned int i;
	for(i = 0;i < m; i++) v[i] = u[i];
}

void copyMatrix(double* B,double* A, unsigned int m, unsigned int n){
	unsigned int i;
	for(i = 0;i < m*n; i++) B[i] = A[i];
}

void matrixMatrixProduct(double* C, double* A, double* B, unsigned int m, unsigned int n, unsigned int p){
	unsigned int i, j ,k;
	double somme;
	for(i=0; i<m; i++){
		for(k=0;k<p;k++){
			somme = 0;
			for(j=0; j<n;j++){
				somme += A[LIN(i,j,m,n)]*B[LIN(j,k,n,p)];
			}
			C[LIN(i,k,m,p)]=somme;
		}
	}
}

void matrixVectorProduct(double* v, double* A, double* u, unsigned int m, unsigned int n){
	unsigned int i, j;
	double somme;
	for(i = 0; i < m; i++){
		somme = 0;
		for(j = 0; j < n; j++){
			somme += A[LIN(i,j,m,n)]*u[j];	
		}
		v[i] = somme;
	}
}

void gaussianElimination(double *A, double *b, unsigned int n){
	
	unsigned int i,j,k;
	double c;
	for(i=0;i<n;i++){
		for(j=i+1;j<n;j++){
			c = A[LIN(j,i,n,n)]/A[LIN(i,i,n,n)];
			for(k=i;k<n;k++){
				A[LIN(j,k,n,n)]-= c*A[LIN(i,k,n,n)];
			}
			b[j] -= c*b[i];
		}
	}
}

void solveTriangularUpper(double* x, double* A, double* b, unsigned int n){
	// On suppose n déjà triangularisé
	x[n-1] = b[n-1]/A[LIN(n-1,n-1,n,n)];
	unsigned int i,j;
	for(i=n-2; i>=0; i--){
		for(j=i+1;j<n;j++){
			x[i]+= x[j]*A[LIN(i,j,n,n)];
		}
		x[i] = (b[i]-x[i])/A[LIN(i,i,n,n)];

		if(i==0) break;
	}
}
