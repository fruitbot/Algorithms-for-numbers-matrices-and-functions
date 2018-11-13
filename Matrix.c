# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <math.h>
# include <inttypes.h>

#define RANDMAX 10000000
#define LIN(i,j,m,n) ((i)*(n)+(j))
#define MIN(a,b) ((a)>(b)? (b):(a))


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
			A[LIN(p,k-1,m,n)] = v[p];
	}
}

void setMatrixRow(double* A, double* v, unsigned int k, unsigned int m, unsigned int n){
	unsigned int p;
	for(p = 0; p < n; p++){
			A[LIN(k-1,p,m,n)] = v[p];
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


int main(){
	double* TabV1 = allocateVector(5); 
	double* TabV2 = allocateVector(5);
	
	double* TabM1 = allocateMatrix(3, 5);
	double* TabM2 = allocateMatrix(3, 5);
	
	//Initialisation avec fonctions Random
	randomVector(TabV1, 5);
	randomMatrix(TabM1, 3, 5);
	
	printf("TabV1: \n"); printVector(TabV1, 5);
	printf("TabM1: \n"); printMatrix(TabM1, 3, 5);

	//Initialisation à la main
	readVector(TabV2, 5);
	readMatrix(TabM2, 3, 5);
	
	printf("TabV2: \n"); printVector(TabV2, 5);
	printf("TabM2: \n"); printMatrix(TabM2, 3, 5);
		
	
	elementVector(TabV1, 4, 5);
	identityMatrix(TabM1, 3, 5);

	copyVector(TabV2, TabV1, 5);
	copyMatrix(TabM2, TabM1, 3, 5);

	printf("TabV1: \n"); printVector(TabV1, 5);
	printf("TabM1: \n"); printMatrix(TabM1, 3, 5);

	printf("TabV2: \n"); printVector(TabV2, 5);
	printf("TabM2: \n"); printMatrix(TabM2, 3, 5);

	//Test des setMatrix
	printf("Test des setMatrix\n");
	printf("TabM1: \n"); printMatrix(TabM1, 3, 5);
	setMatrixColumn(TabM1, TabV1, 2, 3, 5);
	setMatrixRow(TabM1, TabV1, 2, 3, 5);

	printf("TabM1: \n"); printMatrix(TabM1, 3, 5);

	//Test multiplication par un scalaire
	printf("Test multiplication par un scalaire\n");
	scaleVector(2, TabV1,5);
	scaleMatrix(2, TabM1, 3,5);

	printf("TabV1: \n"); printVector(TabV1, 5);
	printf("TabM1: \n"); printMatrix(TabM1, 3, 5);

	//Test add
	printf("Test add\n");
	double* TabV3 = allocateVector(5); double* TabM3 = allocateMatrix(3,5);
	
	addVector(TabV3, TabV1, TabV2, 5);
	addMatrix(TabM3, TabM1, TabM2, 3,5);
	
	printf("TabV3: \n"); printVector(TabV3, 5);
	printf("TabM3: \n"); printMatrix(TabM3, 3, 5);

	//Test produit scalaire de vecteurs
	printf("Test produit scalaire de vecteurs\n");
	printf("TabV1: \n"); printVector(TabV1, 5);
	printf("TabV2: \n"); printVector(TabV2, 5);

	printf("TabV1.TabV2= \n"); printf("%lf\n",scalarProduct(TabV1, TabV2, 5));

	//Test produit matrice vecteur
	printf("Test produit matrice vecteur\n");
	double* TabV4 = allocateVector(3);
	printf("TabV1: \n"); printVector(TabV1, 5);
	printf("TabM1: \n"); printMatrix(TabM1, 3, 5);
	matrixVectorProduct(TabV4, TabM1, TabV1, 3, 5);

	printf("TabV4: \n"); printVector(TabV4, 3);

	//Test produt matrice matrice 
	printf("Test produit matrice matrice\n");
	double* TabM4 = allocateMatrix(3, 5);
	double* TabM5 = allocateMatrix(5, 4);
	double* TabM6 = allocateMatrix(3, 4);
	//randomMatrix(TabM4, 3, 5);
	//randomMatrix(TabM5, 5, 4);
	readMatrix(TabM4, 3, 5);
	readMatrix(TabM5, 5, 4);
	printf("TabM4: \n"); printMatrix(TabM4, 3, 5);
	printf("TabM5: \n"); printMatrix(TabM5, 5, 4);

	matrixMatrixProduct(TabM6, TabM4, TabM5, 3, 5, 4);
	printf("TabM6: \n"); printMatrix(TabM6, 3, 4);
	
	freeVector(TabV1);
	freeVector(TabV2);
	freeVector(TabV3);
	freeVector(TabV4);
	freeMatrix(TabM1);
	freeMatrix(TabM2);
	freeMatrix(TabM3);
	freeMatrix(TabM4);
	freeMatrix(TabM5);
	freeMatrix(TabM6);	

	return 0;
}

