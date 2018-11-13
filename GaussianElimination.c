# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <math.h>
# include <inttypes.h>

#define RANDMAX 10000000
#define LIN(i,j,m,n) ((i)*(n)+(j))
#define MIN(a,b) ((a)>(b)? (b):(a))

/* Signatures fonctions de manipulation de matrices et vecteurs*/

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


/* Exercice 3 */

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

void solveSystemGauss(double* x, double* A, double* b, unsigned int n){
	double *b1 = allocateVector(n);
	double *A1 = allocateMatrix(n, n);
	
	copyVector(b1, b, n);
	copyMatrix(A1, A, n, n);
	
	gaussianElimination(A1, b1, n);
	solveTriangularUpper(x, A1, b1, n);
}



int main(){

	//Initialisation à la main
	//readVector(TabV1, 4);
	//readMatrix(TabM1, 4, 4);
	
	//Initialisation avec Exo3 1. en TD
	double TabV1[] = {7,29,17,-23};
	double TabM1[] = {2,1,1,-3,6,2,5,-8,4,3,3,-9,-2,-2,-5,10};

	printf("TabV1: \n"); printVector(TabV1, 4);
	printf("TabM1: \n"); printMatrix(TabM1, 4, 4);

	//Test gaussianElimination
	gaussianElimination(TabM1, TabV1, 4);	
	printf("TabV1: \n"); printVector(TabV1, 4);
	printf("TabM1: \n"); printMatrix(TabM1, 4, 4);

	//Test solveTriangularUpper
	double *X = allocateVector(4);
	solveTriangularUpper(X, TabM1, TabV1, 4);
	printf("X: \n"); printVector(X, 4);

	//Test solveSystemGauss
	double *X1 = allocateVector(4);
	solveSystemGauss(X1, TabM1, TabV1, 4);
	printf("X1: \n"); printVector(X1, 4);


	freeVector(X);
	freeVector(X1);

	//A dé-commenter si on a initialisé TabV1 et TabM1 à la main
	//freeVector(TabV1);
	//freeMatrix(TabM1);
	
	return 0;

}


/* Fonctions manipulation matrices et vecteurs*/

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
