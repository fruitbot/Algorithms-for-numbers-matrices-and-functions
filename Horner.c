# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <time.h>
# include <math.h>
# include <inttypes.h>

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


double f(double d){
	return 0;
}

int main(){
	FILE* fd = fopen("poly.dat", "w+");
	unsigned int n;
	printf("Donnez un degré n\n");
	scanf("%u", &n);
	double *c = (double*)malloc(sizeof(double)*(n+1));
	unsigned int i, k;
	printf("Entrez les coefficients\n");
	for(i=0; i<n+1; i++){
		scanf(" %lf", &c[i]);
	}
	
	printf("Entrez [a, b] et k\n");
	double a, b;
	scanf("%lf", &a);
	scanf("%lf", &b);
	scanf("%u", &k);
	/*double *f = (double*)malloc(sizeof(double)*(n+1));
	for(i=0; i<n; i++){ 
		f[i]=0;
	}*/
	
	printPolynomialApproxError(fd, c, n, f, a, b, k);
	free(c);
	fclose(fd);
	return 0;
}
