# include <stdio.h>
# include <stdint.h>
# include <inttypes.h>

/* Pour la fonction ci-dessous on va implémenter un tableau comme vu en cours et en td, dont les colonnes correspondront respectivement à u1, u2, u3,, v1, v2, v3, q, t1, t2, t3 .
Au final, on aura u = t1; v = t2; pgcd = dernier reste non-nul = v3*/

void ext_eucl_div(int64_t *u, int64_t *v, int64_t *g, int64_t a, int64_t b){
	int64_t u1, u2, u3 , v1, v2, v3, q, t1, t2, t3;
	int tour = 0;
	do{
		if(tour == 0){
			u1 = 1; u2 = 0; u3 = a; v1 = 0; v2 = 1; v3 = b;
		}
		else{
			u1 = v1; u2 = v2; u3 = v3; v1 = t1; v2 = t2; v3 = t3;
		}
		q = u3/v3;
		t1 = u1 - q*v1;
		t2 = u2 - q*v2;
		t3 = u3%v3;
		tour++;
	} while(t3>=1);

	*u = v1;
	*v = v2;
	*g = v3;

}

int64_t gcd(int64_t a, int64_t b){
	int64_t u, v, g;
	ext_eucl_div(&u, &v, &g, a, b);
	return g;
}

int modular_inverse(int64_t *i, int64_t a, int64_t m){

	if(gcd(a,m)!=1) return 0;

	int64_t u, v, g;
	ext_eucl_div(&u, &v, &g, a, m);	
	*i = u;
	return 1;
}



int main(){
	int64_t a, b, u = 0, v = 0, pgcd = 1, m, i;
	printf("Calcul du PGCD: veuillez choisir deux entiers a et b:\na= ");
	scanf("%" PRId64,&a);
	printf("b= ");
	scanf("%" PRId64,&b);
	ext_eucl_div(&u, &v, &pgcd, a, b);
	printf("D'après la relation de Bezout, %" PRId64, a);
	printf(" * %" PRId64, u);
	printf(" + %" PRId64, b);
	printf(" * %" PRId64, v);
	printf(" = %" PRId64 "\n", gcd(a,b));
	
	printf("Calcul de l'inverse modulaire: veuillez choisir deux entiers a et m:\na= ");
	scanf("%" PRId64,&a);
	printf("m= ");
	scanf("%" PRId64,&m);
	if(modular_inverse(&i,a,m)){
		printf("L'inverse modulaire de %" PRId64, a);
		printf(" modulo %" PRId64, m);
		printf(" existe et est égal à: %" PRId64 "\n", i);
	} else {
		printf("Pas d'inverse modulaire.");
	}

	return 0;
}
