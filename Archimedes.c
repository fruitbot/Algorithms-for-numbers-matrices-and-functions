# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
#define _USE_MATH_DEFINES
# include <math.h>

double archimedes1(unsigned int k){
	if(k<2) return k;
	
	double L2k = sqrt(2);
	if(k==2) return L2k;
	int i;

	for(i = 3; i <= k; i++){
		L2k = sqrt(2*(1-sqrt(1-pow(L2k,2)/4.)));
		//printf("k=L2k1 = %f\n", L2k);
	}
		
	return L2k;

}

double archimedes2(unsigned int k){
	if(k<2) return k;
	
	double L2k = sqrt(2);
	if(k==2) return L2k;
	int i;

	for(i = 3; i <= k; i++){
		L2k = sqrt(4*((pow(L2k,2)/4.)/(2*(1+sqrt(1-pow(L2k,2)/4.)))));
	}

	return L2k;

}


int main(){
	int k;

	for(k = 10; k<=30; k++){
		printf("k = %d\n", k);	
		printf("Formule 1 : pi = %.25f\n", archimedes1(k)*pow(2,k-1));
		printf("Formule 2 : pi = %.25f\n", archimedes2(k)*pow(2,k-1));
		printf("M_PI : pi = %.25f\n", M_PI);
	}
	

	/* Commentaire pour Ex7, Q3:
	On teste les approximations citées dans notre réponse à la q3 de l'exercice 7.
	En prenant un Ln très petit:
	TEST FORMULE 1:
	printf("%.50f\n", sqrt(2*(1-sqrt(1-pow(0.0000000149011611938476562500000000000000,2)/4.))));
	>>> 1.00000000000000000000000000000000000000000000000000
	TEST FORMULE 2:
	printf("%.50f\n", 0.0000000149011611938476562500000000000000*pow(2,28));
	>>> 0.00000000000000000093132257461547857362925210760893
	*/

	return 0;
}
