# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <time.h>
# include <math.h>
# include <inttypes.h>


/* ------------------------------------ FONCTIONS EXERCICE 1 ------------------------------------------------------*/


int findZeroBisection(double *z, double(*f)(double), double a, double b, double eps){

	if(f(a)*f(b)>0) return 0;
	
	double a0 = a, b0 = b;
	*z = (a0 + b0)/2.;

	while(fabs(f(*z))>eps){
		if(f(a0)*f(*z)>=0) a0 = *z;
		else b0 = *z;
		*z = (a0 + b0)/2.;
	}

	return 1;
}

int findZeroNewton(double *z, double(*f)(double), double(*fprime)(double), double a, double b, double eps){

	double xk = (a+b)/2.;

	while(fabs(f(xk))>eps){
		
		if((fprime(xk)==0)||(xk<a)||(xk>b)){
			printf("Je suis ici\n");
			return 0;
		}	
		else xk -= ((double)f(xk))/fprime(xk);

	}
	*z = xk;

	return 1;

}


int findZero(double *z, double(*f)(double), double(*fprime)(double), double a, double b, double eps){

	double a0 = a, b0 = b;
	double xk = (a+b)/2.;

	while(fabs(f(xk))>eps){
		if((fprime(xk)==0)||(xk<a0)||(xk>b0)){
			if(f(a)*f(b)>0) return 0;
			if(f(a0)*f(xk)>=0) a0 = xk;
			else b0 = xk;
			xk = (a0 + b0)/2.;			
		}
		else xk -= ((double)f(xk))/fprime(xk);
	}

	*z = xk;
	return 1;
}



/* ------------------------------------------------- FONCTIONS TESTS ----------------------------------------------------*/

/* NOTE: On met tous les polynômes sous la forme de polynômes de Horner */

double f(double x){
	return exp(x)-4;
}

double fprime(double x){
	return exp(x);
}

double g(double x){
	return sin(x-3);
}

double gprime(double x){
	return cos(x-3);
}

double h(double x){
	return 3 + x*(-18 + x*(162*x));
}

double hprime(double x){
	return -18 + x*(486*x);
}

double p(double x){
	return x*(x*(2-x*(1+x*(2-x))));
}


double pprime(double x){
	return x*(4-x*(3+x*(8-5*x)));
}

double q(double x){
	return atan(x+1);
}

double qprime(double x){
	return 1./(2+x*(2+x));
}

/* ------------------------------------------------ TEST MAIN ------------------------------------------------------------ */

int main(){

	double eps = 0.000030517578125;
	double z;

	//Test findZeroBisection
	printf("TEST findzeroBisection\n\n");
	if(findZeroBisection(&z, f, 1, 1.5, eps)) printf("Le zéro de f est %lf\n", z);
	else printf("Erreur pour f\n");
	if(findZeroBisection(&z, g, -0.5, 0.5, eps)) printf("Le zéro de g est %lf\n", z);
	else printf("Erreur pour g\n");
	if(findZeroBisection(&z, h, -1./3, 1, eps)) printf("Le zéro de h est %lf\n", z);
	else printf("Erreur pour h\n");
	if(findZeroBisection(&z, p, -1.5, 0.5, eps)) printf("Le zéro de p est %lf\n", z);
	else printf("Erreur pour p\n");
	if(findZeroBisection(&z, q, -2, 20, eps)) printf("Le zéro de q est %lf\n", z);
	else printf("Erreur pour q\n");

	//Test findZeroNewton
	printf("\n\nTEST findZeroNewton\n\n");
	if(findZeroNewton(&z, f, fprime, 1, 1.5, eps)) printf("Le zéro de f est %lf\n", z);
	else printf("Erreur pour f\n");
	if(findZeroNewton(&z, g, gprime, -0.5, 0.5, eps)) printf("Le zéro de g est %lf\n", z);
	else printf("Erreur pour g\n");
printf("ca donne: %d\n", findZeroNewton(&z, h, hprime, -1./3, 1, eps));
	if(findZeroNewton(&z, h, hprime, -1./3, 1, eps)) printf("Le zéro de h est %lf\n", z);
	else printf("Erreur pour h\n");
	if(findZeroNewton(&z, p, pprime, -1.5, 0.5, eps)) printf("Le zéro de p est %lf\n", z);
	else printf("Erreur pour p\n");
	if(findZeroNewton(&z, q, qprime, -2, 20, eps)) printf("Le zéro de q est %lf\n", z);
	else printf("Erreur pour q\n");

	//Test findZero
	printf("\n\nTEST findZero\n\n");
	if(findZero(&z, f, fprime, 1, 1.5, eps)) printf("Le zéro de f est %lf\n", z);
	else printf("Erreur pour f\n");
	if(findZero(&z, g, gprime, -0.5, 0.5, eps)) printf("Le zéro de g est %lf\n", z);
	else printf("Erreur pour g\n");
	if(findZero(&z, h, hprime, -1./3, 1, eps)) printf("Le zéro de h est %lf\n", z);
	else printf("Erreur pour h\n");
	if(findZero(&z, p, pprime, -1.5, 0.5, eps)) printf("Le zéro de p est %lf\n", z);
	else printf("Erreur pour p\n");
	if(findZero(&z, q, qprime, -2, 20, eps)) printf("Le zéro de q est %lf\n", z);
	else printf("Erreur pour q\n");
	
	return 0;
	

}
