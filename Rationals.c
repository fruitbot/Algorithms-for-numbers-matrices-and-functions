# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <inttypes.h>


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
	ext_eucl_div(&u, &v, &g, abs(a), abs(b));
	return g;
}

typedef struct __ratio_struct_t{
	int64_t num;
	int64_t den;
} ratio_t;


ratio_t createRatio(int64_t a, int64_t b){
	if(b == 0){
		printf("Erreur : a divise par 0\n");
		exit(1);
	}
	ratio_t r;
	int64_t pgcd = gcd(a, b);
	r.num = a/pgcd;
	r.den = b/pgcd;
	return r;
}

ratio_t addRatio(ratio_t a, ratio_t b){
	int64_t num, den;
	num = a.num*b.den + b.num*a.den;
	den = a.den * b.den;
	return createRatio(num, den);
}


ratio_t subRatio(ratio_t a, ratio_t b){
	int64_t num, den;
	num = a.num*b.den - b.num*a.den;
	den = a.den * b.den;
	return createRatio(num, den);
}

		
ratio_t mulRatio(ratio_t a, ratio_t b){
	return createRatio(a.num*b.num, a.den*b.den);
}

ratio_t divRatio(ratio_t a, ratio_t b){
	ratio_t r = createRatio(b.den,b.num);
	return mulRatio(a,r);
}

double approxRatio(ratio_t a){
	return ((int)((((double)a.num)/(double)a.den)*100))/100.0;
}

void printRatio(ratio_t a){
	printf("%" PRId64, a.num);
	printf("/%" PRId64 "\n", a.den);
}

ratio_t computeS(unsigned int n){

	ratio_t res;
	res.num = 1;
	res.den = 1;
	ratio_t res1;
	if(n == 0) return (res);

	int j = 1;
	ratio_t somme;

	while(j <= n){
		res.num = 1;
		res.den = 1;
		int k = 1;
		while(k <= j){
			res1.num = 1;
			res1.den = k;
			res = mulRatio(res,res1);
			k++;
		}
		if(j == 1 ) somme = res; 
		else{
			somme = addRatio(somme, res);
		}

		j++;
	}	

	return somme;
}

ratio_t computeA(unsigned int n){

	if(n == 0) return createRatio(11,2);
	if(n == 1) return createRatio(61,11);
	int i;
	ratio_t a0 = createRatio(11,2), a1 = createRatio(61,11);
	ratio_t an;
	for(i = 2; i <= n; i++){
		an = subRatio(createRatio(111,1),divRatio(subRatio(createRatio(1130,1),divRatio(createRatio(3000,1),a0)),a1));
		a0 = a1;
		a1 = an;
	}

	return an;

}



int main(){

	// Test opérations sur les fractions

	ratio_t a = {2,4}, b = {5,7};
	ratio_t ex = addRatio(a,b);
	printf("Addition: ");
	printRatio(ex);

	ex = subRatio(a,b);
	printf("Soustraction: ");
	printRatio(ex);

	ex = mulRatio(a,b);
	printf("Multiplication: ");
	printRatio(ex);

	ex = divRatio(a,b);
	printf("Division: ");
	printRatio(ex);

	// Test double précision
	printf("Précision machine: %lf\n", ((double)b.num)/b.den);
	printf("Approx ratio: %lf\n", approxRatio(b));

	// Test série s
	unsigned int n = 3;
	printRatio(computeS(n));

	// Test suite a
	unsigned int n1 = 3;
	printRatio(computeA(n1));
	
	int i=0;
	while(i<= 20){
		printf("Pour n = %d: \n", i);
		printRatio(computeA(i));
		printf("%f \n",approxRatio(computeA(i)));
		i += 1;
	}
	
	/* On voit que la suite an tend vers 5,78 (pour un n maximum de n=7, puis les nombres du dénominateur et
	 du numérateur dépassent la plage possible (val max = 2⁶⁴-1 environ 1*10¹⁹) en 64 bit*/

	return 0;
}
