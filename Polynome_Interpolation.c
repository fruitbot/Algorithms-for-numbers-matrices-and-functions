# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <time.h>
# include <math.h>
# include <inttypes.h>
# include <string.h>


#define RANDMAX 10000000
#define LIN(i,j,m,n) ((i)*(n)+(j))
#define MIN(a,b) ((a)>(b)? (b):(a))

/* ----------------------------SIGNATURE FONCTIONS MANIPULATION MATRICES -------------*/

double *allocateVector(unsigned int m);
double* allocateMatrix(unsigned int m, unsigned int n);
void freeVector(double* v);
void freeMatrix(double* A);
void copyVector(double* v,double* u, unsigned int m);
void copyMatrix(double* B,double* A, unsigned int m, unsigned int n);
void elementVector(double* v, unsigned int k, unsigned int m);
void identityMatrix(double *A, unsigned int m, unsigned int n);
double randomValue();
void randomVector(double *v, unsigned int m);
void randomMatrix(double *A, unsigned int m, unsigned int n);
void readVector(double *v, unsigned int m);
void readMatrix(double *A, unsigned int m, unsigned int n);
void printVector(double *v, unsigned int m);
void printMatrix(double *A, unsigned int m, unsigned int n);
void luDecomposition(double *A, unsigned int n);
void setMatrixColumn(double* A, double* v, unsigned int k, unsigned int m, unsigned int n);
void solveTriangularLower(double* x, double* A, double* b, unsigned int n);
void solveTriangularUpper(double* x, double* A, double* b, unsigned int n);
void finalSolveSystemLU(double* x, double* A, double* b, double *scratch, unsigned int n);
void solveSystemLU(double* x, double* A, double* b, unsigned int n);
double evaluatePolynomial(double *c, double x, unsigned int n);
void printPolynomialApproxError(FILE *fd, double *c, unsigned int n, double (*f)(double), double a, double b, unsigned int k);

/* ------------------------------------TME ------------------------------------------*/


void evaluateFonction(double *y, double (*f)(double), double *x, unsigned int n){
	unsigned int i;
	for(i=0; i<n; i++){
		y[i] = f(x[i]);
	}
}

void equidistantPoints(double *x, double a, double b, unsigned int n){

	unsigned int j;
	for(j=0; j<n; j++){
		x[j] = a+j*(((double)(b-a))/(n-1));	
	}
}

void interpolatePoints(double *c, double *x, double *y, unsigned int n){
	
	double *V = allocateMatrix(n,n);
	double *Vi = allocateVector(n);
	

	//Matrice de Vandermonde
	unsigned int i, j;
	for(i=0; i< n; i++){
		Vi[i] = 1;
	}
	setMatrixColumn(V, Vi, 0, n, n);
	for(i=1; i<n; i++){
		for(j=0; j<n; j++){
			Vi[j] = Vi[j]*x[j];
		}
		setMatrixColumn(V, Vi, i, n, n);
	}


	solveSystemLU(c, V, y, n);

	freeMatrix(V);
	freeVector(Vi);
}

void interpolateFunction(double *c, double *x, double (*f)(double), unsigned int n){
	double *y = allocateVector(n);
	evaluateFonction(y, f, x, n);
	interpolatePoints(c, x, y, n);
	freeVector(y);
}

void interpolateFunctionEquidistantPoints(double *c,double (*f)(double), double a, double b, unsigned int n){
	double *x = allocateVector(n);

	unsigned int i;
	for(i=0; i<n; i++) x[i] = 0;
	equidistantPoints(x, a, b, n);
	interpolateFunction(c, x, f, n);
	freeVector(x);
}



/* -------------------- FONCTIONS TESTS ---------------------------*/

double my_cos(double d){
	return cos(d);
}

double my_sin(double d){
	return sin(d);
}

double my_tan(double d){
	return tan(d);
}

double my_acos(double d){
	return acos(d);
}

double my_asin(double d){
	return asin(d);
}

double my_atan(double d){
	return atan(d);
}

double my_exp(double d){
	return exp(d);
}


/* ------------------------- MAIN -----------------------------------*/

int main(){

	/*Test interpolatePoints*/
	double *c = allocateVector(4);
	double x[] = {-1,0,1,2};
	double y[] = {6,-24,6,-6};
	interpolatePoints(c, x, y, 4);
	printVector(c, 4);

	/*QUESTION 2: Test interpolateFunctionEquidistantPoints*/

	FILE *fd = fopen("poly.dat", "w");

	
	unsigned int n[] = {3,7,17,23,42};
	unsigned int i;
	
	for(i=0; i<5; i++){	
		double *c1 = allocateVector(n[i]);
		interpolateFunctionEquidistantPoints(c1, my_exp, -1, 1, n[i]);
		printPolynomialApproxError(fd, c1, n[i], my_exp, -1, 1, 5);
		freeVector(c1);
	}
	

	/*
	unsigned int n=42;

	
	double *c1 = allocateVector(n);
	interpolateFunctionEquidistantPoints(c1, my_exp, -1, 1, n);
	printPolynomialApproxError(fd, c1, n, my_exp, -1, 1, 1000);
	freeVector(c1);
	*/

	//Désallocation mémoire

	freeVector(c);
	fclose(fd);

	return 0;
}


/* ---------------------FONCTIONS MANIPULATION DE MATRICES-------------------------*/


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

void copyVector(double* v,double* u, unsigned int m){
	unsigned int i;
	for(i = 0;i < m; i++) v[i] = u[i];
}

void copyMatrix(double* B,double* A, unsigned int m, unsigned int n){
	unsigned int i;
	for(i = 0;i < m*n; i++) B[i] = A[i];
}

void elementVector(double* v, unsigned int k, unsigned int m){
	unsigned int i;
	for(i = 0; i < m; i++){
		if(i == k) v[i] = 1;
		else v[i] = 0;
	}
}

void identityMatrix(double *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i = 0;i < m*n; i++) A[i] = 0;
	for(j = 0; j < MIN(m,n); j++){
			A[LIN(j,j,m,n)] = 1;
	}
}


double randomValue(){
	return 2*((double)(rand()))/((double)RANDMAX)-1;
}

void randomVector(double *v, unsigned int m){
	unsigned int i;
	for(i = 0; i < m; i++) v[i] = randomValue();
}

void randomMatrix(double *A, unsigned int m, unsigned int n){
	unsigned int i,j;
	for(i = 0;i < m*n; i++) A[i] = randomValue()/4;
	for(j = 0; j < MIN(m,n); j++){
			A[LIN(j,j,m,n)] = (randomValue()+1)/4 + 0.75;
	}
}

void readVector(double *v, unsigned int m){
	double d;
	unsigned int i;	
	printf("Saisie du vecteur (double):");
	for(i = 0; i < m; i++){	
		if(scanf("%lf",&d)!=1) exit(1);
		else v[i] = d;
	}
}

void readMatrix(double *A, unsigned int m, unsigned int n){
	double d;
	unsigned int i;
	printf("Saisie de la matrice (double):");
	for(i = 0; i < m*n; i++){	
		if(scanf("%lf",&d)!=1) exit(1);
		else A[i] = d;
	}
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

double maximumAbsVector(double* v, unsigned int m){
	double max = v[0];
	unsigned int i;
	for(i=0;i<m;i++){
		if(fabs(max)<fabs(v[i])) max = v[i];
	}
	return max;
}

double maximumAbsMatrix(double* A, unsigned int m, unsigned int n){
	double max = A[0];
	unsigned int i;
	for(i=0;i<m*n;i++){
		if(max<fabs(A[i])) max = A[i];
	}
	return max;
}

void setMatrixColumn(double* A, double* v, unsigned int k, unsigned int m, unsigned int n){
	unsigned int p;
	for(p = 0; p < m; p++){
			A[LIN(p,k,m,n)] = v[p];
	}
}

void setMatrixRow(double* A, double* v, unsigned int k, unsigned int m, unsigned int n){
	unsigned int p;
	for(p = 0; p < n; p++){
			A[LIN(k,p,m,n)] = v[p];
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

void addVector(double* w, double* u, double* v, unsigned int m){
	unsigned int i;
	for(i=0; i<m; i++) w[i] = u[i] + v[i];
}

void addMatrix(double* C, double* A, double* B, unsigned int m, unsigned int n){
	unsigned int i;
	for(i=0; i<m*n; i++) C[i] = A[i] + B[i];
}

double scalarProduct(double* u, double* v, unsigned int m){
	double ps = 0;
	unsigned int i;
	for(i=0; i<m; i++) ps += u[i] * v[i];
	return ps;
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

double evaluatePolynomial(double *c, double x, unsigned int n){
	double y;
	unsigned int i;
	y = 0.0;
	for(i=n; i>0; i--){
		y = c[i]+x*y;
	}
	return y;
}

void printPolynomialApproxError(FILE *fd, double *c, unsigned int n, double (*f)(double), double a, double b, unsigned int k){
	unsigned int j;
	double x_j, delta_j;
	for(j=0; j<k; j++){
		x_j = a+j*(((double)(b-a))/(k-1));
		delta_j = evaluatePolynomial(c, x_j, n) - f(x_j);
		fprintf(fd, "%+1.18e\t%+1.18e\n", x_j, delta_j);
	}
}
