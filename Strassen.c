# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <time.h>
# include <math.h>
# include <inttypes.h>

#define RANDMAX 10000000
#define LIN(i,j,m,n) ((i)*(n)+(j))
#define MIN(a,b) ((a)>(b)? (b):(a))

double* allocateMatrix(unsigned int m, unsigned int n);
void freeMatrix(double* A);
void printMatrix(double *A, unsigned int m, unsigned int n);
void scaleMatrix(double alpha, double* A, unsigned int m, unsigned int n);
void addMatrix(double* C, double* A, double* B, unsigned int m, unsigned int n);
void matrixMatrixProduct(double* C, double* A, double* B, unsigned int m, unsigned int n, unsigned int p);
double randomValue();
void randomMatrix(double *A, unsigned int m, unsigned int n);


void splitMatrix(double *B, double *C, double *D, double *E, double *A, unsigned int m, unsigned int n, unsigned int p, unsigned int q){
	unsigned int i,j;

	for(i=0;i<p;i++){
		for(j=0;j<q;j++){
			B[LIN(i,j,p,q)] = A[LIN(i,j,m,n)];
		}
	}

	for(i=0;i<p;i++){
		for(j=0;j<n-q;j++){
			C[LIN(i,j,p,n-q)] = A[LIN(i,q+j,m,n)];
		}
	}

	for(i=0;i<m-p;i++){
		for(j=0;j<q;j++){
			D[LIN(i,j,m-p,q)] = A[LIN(i+p,j,m,n)];
		}
	}

	for(i=0;i<m-p;i++){
		for(j=0;j<n-q;j++){
			E[LIN(i,j,m-p,n-q)] = A[LIN(i+p,j+q,m,n)];
		}
	}	
}



void composeMatrix(double *A, double *B, double *C, double *D, double *E, unsigned int m, unsigned int n, unsigned int p, unsigned int q){
	unsigned int i,j;

	for(i=0;i<p;i++){
		for(j=0;j<q;j++){
			A[LIN(i,j,m,n)] = B[LIN(i,j,p,q)];
		}
	}

	for(i=0;i<p;i++){
		for(j=0;j<n-q;j++){
			A[LIN(i,q+j,m,n)] = C[LIN(i,j,p,n-q)];
		}
	}

	for(i=0;i<m-p;i++){
		for(j=0;j<q;j++){
			A[LIN(i+p,j,m,n)] = D[LIN(i,j,m-p,q)];
		}
	}

	for(i=0;i<m-p;i++){
		for(j=0;j<n-q;j++){
			A[LIN(i+p,j+q,m,n)] = E[LIN(i,j,m-p,n-q)];
		}
	}

}

void matrixMatrixProductStrassen(double *C,double *B,double *A,unsigned int k){

	double pow0 = pow(2.,(double)k);
	if(k<=2){
		matrixMatrixProduct(C, B, A, pow0, pow0, pow0);
	}
	// k'
	else{

	double pow1 = pow(2.,(double)(k-1));
	double *A11 = allocateMatrix(pow1,pow1);
	double *A12 = allocateMatrix(pow1,pow1);
	double *A21 = allocateMatrix(pow1,pow1);
	double *A22 = allocateMatrix(pow1,pow1);

	splitMatrix(A11, A12, A21, A22, A, pow0, pow0, pow1, pow1);

	double *B11 = allocateMatrix(pow1,pow1);
	double *B12 = allocateMatrix(pow1,pow1);
	double *B21 = allocateMatrix(pow1,pow1);
	double *B22 = allocateMatrix(pow1,pow1);

	splitMatrix(B11, B12, B21, B22, B, pow0, pow0, pow1, pow1);

	double *P11 = allocateMatrix(pow1,pow1);
	double *P12 = allocateMatrix(pow1,pow1);
	double *P21 = allocateMatrix(pow1,pow1);
	double *P22 = allocateMatrix(pow1,pow1);

	double *T1 = allocateMatrix(pow1,pow1);
	double *T2 = allocateMatrix(pow1,pow1);
	double *T3 = allocateMatrix(pow1,pow1);
	double *T4 = allocateMatrix(pow1,pow1);
	double *T5 = allocateMatrix(pow1,pow1);
	double *T6 = allocateMatrix(pow1,pow1);
	double *T7 = allocateMatrix(pow1,pow1);

	//Calcul de T1
	double *T11 = allocateMatrix(pow1,pow1);
	double *T12 = allocateMatrix(pow1,pow1);

	addMatrix(T11, B11, B22, pow1, pow1);
	addMatrix(T12, A11, A22, pow1, pow1);
	matrixMatrixProductStrassen(T1, T11, T12, k-1);

	//Calcul de T2
	double *T21 = allocateMatrix(pow1,pow1);

	addMatrix(T21, B21, B22, pow1, pow1);
	matrixMatrixProductStrassen(T2, T21, A11, k-1);

	//Calcul de T3
	double *T31 = allocateMatrix(pow1,pow1);
	scaleMatrix(-1, A22, pow1, pow1);
	addMatrix(T31, A12, A22, pow1, pow1);
	matrixMatrixProductStrassen(T3, B11, T31, k-1);
	scaleMatrix(-1, A22, pow1, pow1);
	
	//Calcul de T4
	double *T41 = allocateMatrix(pow1,pow1);
	scaleMatrix(-1, A11, pow1, pow1);
	addMatrix(T41, A21, A11, pow1, pow1);
	matrixMatrixProductStrassen(T4, B22, T41, k-1);
	scaleMatrix(-1, A11, pow1, pow1);

	//Calcul de T5
	double *T51 = allocateMatrix(pow1,pow1);

	addMatrix(T51, B11, B12, pow1, pow1);
	matrixMatrixProductStrassen(T5, T51, A22, k-1);

	//Calcul de T6
	double *T61 = allocateMatrix(pow1,pow1);
	double *T62 = allocateMatrix(pow1,pow1);


	scaleMatrix(-1, B11, pow1, pow1);

	addMatrix(T61, B21, B11, pow1, pow1);
	addMatrix(T62, A11, A12, pow1, pow1);
	matrixMatrixProductStrassen(T6, T61, T62, k-1);
	scaleMatrix(-1, B11, pow1, pow1);

	//Calcul de T7
	double *T71 = allocateMatrix(pow1,pow1);
	double *T72 = allocateMatrix(pow1,pow1);

	scaleMatrix(-1, B22, pow1, pow1);

	addMatrix(T71, B12, B22, pow1, pow1);
	addMatrix(T72, A21, A22, pow1, pow1);
	matrixMatrixProductStrassen(T7, T71, T72, k-1);
	scaleMatrix(-1, B22, pow1, pow1);

	//Calcul des P11
	
	scaleMatrix(-1, T5, pow1, pow1);
	addMatrix(P11, T1, T4, pow1, pow1);
	addMatrix(P11, P11, T5, pow1, pow1);
	addMatrix(P11, P11, T7, pow1, pow1);
	scaleMatrix(-1, T5, pow1, pow1);

	//Calcul des P22
	
	scaleMatrix(-1, T2, pow1, pow1);
	addMatrix(P22, T1, T2, pow1, pow1);
	addMatrix(P22, P22, T3, pow1, pow1);
	addMatrix(P22, P22, T6, pow1, pow1);
	scaleMatrix(-1, T2, pow1, pow1);

	//Calcul des P12
	
	addMatrix(P12, T3, T5, pow1, pow1);

	//Calcul des P21
	
	addMatrix(P21, T2, T4, pow1, pow1);

/*	printMatrix(P11, pow1, pow1);
	printMatrix(P12, pow1, pow1);
	printMatrix(P21, pow1, pow1);
	printMatrix(P22, pow1, pow1);
*/

	//Composition matrice finale
	composeMatrix(C, P11, P12, P21, P22, pow0,pow0,pow1,pow1);

	freeMatrix(A11);
	freeMatrix(A12);
	freeMatrix(A21);
	freeMatrix(A22);
	freeMatrix(B11);
	freeMatrix(B12);
	freeMatrix(B21);
	freeMatrix(B22);
	freeMatrix(P11);
	freeMatrix(P12);
	freeMatrix(P21);
	freeMatrix(P22);
	freeMatrix(T1);
	freeMatrix(T2);
	freeMatrix(T3);
	freeMatrix(T4);
	freeMatrix(T5);
	freeMatrix(T6);
	freeMatrix(T7);
	freeMatrix(T11);
	freeMatrix(T12);
	freeMatrix(T21);
	freeMatrix(T31);
	freeMatrix(T41);
	freeMatrix(T51);
	freeMatrix(T61);
	freeMatrix(T62);
	freeMatrix(T71);
	freeMatrix(T72);
	}
}	


int main(){

	srand(time(NULL));

	unsigned int k = 1 + rand()%10;
	double powk = pow(2., (double)k);

	double *A = allocateMatrix(powk, powk);
	double *B = allocateMatrix(powk, powk);	
	double *C = allocateMatrix(powk, powk); //C prendra la valeur de la A*B par la méthode naïve
	double *D = allocateMatrix(powk, powk);	//D prendra la valeur de A*B par la méthode de Strassen

	double diff = 0.; //moyenne des erreurs entre les 2 méthodes pour un certain k
	double moyDiff = 0.; //moyenne des erreurs sur 100 tests

	unsigned int i, j;
	for(i=0; i<100; i++){
		randomMatrix(A, pow(2,k), pow(2,k)); //Matrice aléatoire A
		randomMatrix(B, pow(2,k), pow(2,k)); //Matrice aléatoire B


		//Calcul de A*B par la méthode naïve
		matrixMatrixProduct(C, A, B, powk, powk, powk);	
	
		//Calcul de A*B par la méthode de Strassen
		matrixMatrixProductStrassen(D, A, B, k);

		//Calcul de la moyenne des différences entre les valeurs trouvées pour chaque méthode
		for(j=0; j<powk*powk; j++) moyDiff += fabs(C[j]-D[j]);
		diff /= powk*powk;
		moyDiff += diff;
	}
	
	moyDiff /= 100;
	printf("La moyenne des erreurs entre la multplication naïve et la multiplication de Strassen est de: %.63f\n", moyDiff);
		
	return 0;
}







double* allocateMatrix(unsigned int m, unsigned int n){
	double* TabM = (double*)calloc(m*n, sizeof(double));
	if(TabM == NULL) exit(1);
	return TabM;
}

void freeMatrix(double* A){
	free(A);
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

void scaleMatrix(double alpha, double* A, unsigned int m, unsigned int n){
	unsigned int i;
	for(i=0; i<m*n; i++) A[i] *= alpha;
}

void addMatrix(double* C, double* A, double* B, unsigned int m, unsigned int n){
	unsigned int i;
	for(i=0; i<m*n; i++) C[i] = A[i] + B[i];
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
