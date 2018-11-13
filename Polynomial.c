# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <time.h>
# include <math.h>
# include <inttypes.h>
# include <string.h>

#define MAX(a,b) ((a)>(b)? (a):(b))
#define MIN(a,b) ((a)>(b)? (b):(a))


typedef struct __ratio_struct_t{
	int64_t num;
	int64_t den;
} ratio_t;


typedef struct __poly_struct_t *poly_t;
struct __poly_struct_t{
	unsigned int deg;
	ratio_t *coeffs;
};


/* -------------------------------------- SIGNATURES FONCTIONS TME3 ----------------------------------------------------- */


void ext_eucl_div(int64_t *u, int64_t *v, int64_t *g, int64_t a, int64_t b);
int64_t gcd(int64_t a, int64_t b);
ratio_t createRatio(int64_t a, int64_t b);
ratio_t addRatio(ratio_t a, ratio_t b);
ratio_t subRatio(ratio_t a, ratio_t b);
ratio_t mulRatio(ratio_t a, ratio_t b);

ratio_t divRatio(ratio_t a, ratio_t b);
void printRatio(ratio_t a);
void polyPrint(poly_t p);



/*--------------------------------------------------- TME 10 ----------------------------------------------------------- */


poly_t polyFromRatioArray(ratio_t *c, unsigned int d){

	poly_t p = malloc(sizeof(struct __poly_struct_t));
	p->deg = d;
	ratio_t *coeff = calloc(d+1, sizeof(ratio_t));
	unsigned int i;
	for(i=0; i<=d; i++) coeff[i] = c[i];
	p->coeffs = coeff;
	return p;
}

poly_t polyFromRatio(ratio_t c){
	return polyFromRatioArray(&c, 0);
}

poly_t polyFromInt(int64_t c){
	ratio_t c2 = createRatio(c,(int64_t)1);
	return polyFromRatioArray( &c2, 0);
}


poly_t polyFromMonomialRatio(ratio_t c, unsigned int k){

	ratio_t *coeff = calloc(k+1, sizeof(ratio_t));
	coeff[k] = c;	
	unsigned int i;
	ratio_t zero = createRatio(0, 1);
	for(i = 0; i < k; i++) coeff[i] = zero;
	poly_t m = polyFromRatioArray(coeff, k);
	free(coeff);
	return m;
}

poly_t polyFromIdentity(){
	poly_t pid = (struct __poly_struct_t*)malloc(sizeof(struct __poly_struct_t));
	pid = polyFromMonomialRatio(createRatio(1,1), 0);
	return pid;
}

poly_t polyFromCopy(poly_t p){
	return polyFromRatioArray(p->coeffs, p->deg);
}


void polyFree(poly_t p){
	free(p->coeffs);
	free(p);
}

unsigned int polyGetDegree(poly_t p){
	return p->deg;
}

int polyIsConstantZero(poly_t p){
	ratio_t sum = createRatio(0, 1);
	unsigned int i;
	for (i = 0; i <= p->deg; i++) {
  		sum = addRatio(sum, p->coeffs[i]);
	}
	if (sum.num != 0) {
		return 0;
	}
	return 1;
}

poly_t polyAdd(poly_t p, poly_t q){

	poly_t a = (struct __poly_struct_t*)malloc(sizeof(struct __poly_struct_t));
	if((p->deg - q->deg) != 0){
		unsigned int min_d = MIN(p->deg, q->deg);
		unsigned int max_d = MAX(p->deg, q->deg);
		a->deg = max_d;
		ratio_t *coeff = malloc(sizeof(ratio_t)*((a->deg)+1));
		unsigned int i;
		for(i = 0; i <= min_d; i++){
			coeff[i]=addRatio(p->coeffs[i], q->coeffs[i]);
		}
		for(i = min_d+1; i <= a->deg; i++){
			if(max_d == p->deg) coeff[i]=p->coeffs[i];
			else coeff[i]=q->coeffs[i];
		}
		a->coeffs = coeff;
	}
	else{
		unsigned int i;
		a->deg = p->deg;
		ratio_t *coeff = malloc(sizeof(ratio_t)*((a->deg)+1));
		for(i=0; i <= p->deg; i++) coeff[i] = addRatio(p->coeffs[i], q->coeffs[i]);
		a->coeffs = coeff;
	}

	return a; 
}

poly_t polySub(poly_t p, poly_t q){

	unsigned int min_d = MIN(p->deg, q->deg);
	unsigned int max_d = MAX(p->deg, q->deg);
	poly_t s;
	if(max_d == p->deg){
		s = polyFromCopy(p);
	}
	else{
		 s = polyFromCopy(q);
	}
	if((p->deg - q->deg) != 0){
		s->deg = max_d;
		ratio_t *coeff = malloc(sizeof(ratio_t)*((s->deg)+1));
		unsigned int i;
		for(i = 0; i <= min_d; i++){
			coeff[i]=subRatio(p->coeffs[i], q->coeffs[i]);
		}
		for(i = min_d+1; i <= s->deg; i++){
			if(max_d == p->deg) coeff[i]=p->coeffs[i];
			else coeff[i]=q->coeffs[i];
		}
		s->coeffs = coeff;
	}
	else{
		unsigned int i;
		s->deg = p->deg;
		ratio_t *coeff = malloc(sizeof(ratio_t)*((s->deg)+1));
		for(i=0; i <= p->deg; i++) coeff[i] = subRatio(p->coeffs[i], q->coeffs[i]);
		s->coeffs = coeff;	
	}
	return s;
}


poly_t polyMul(poly_t p, poly_t q){

	poly_t m = (struct __poly_struct_t*)malloc(sizeof(struct __poly_struct_t));
	unsigned int i, j;
	ratio_t *coeff = malloc(sizeof(ratio_t)*((m->deg)+1));

	for(i = 0; i <= p->deg; i++){
		for(j = 0; j <= q->deg; j++){
		coeff[i+j] = createRatio(0,1);
		}
	}
	m->deg = p->deg + q->deg;
	for(i = 0; i <= p->deg; i++){
		for(j = 0; j <= q->deg; j++){
			coeff[i+j]= addRatio(coeff[i+j], mulRatio(p->coeffs[i], q->coeffs[j]));
		}
	}
	m->coeffs = coeff;
	return m;
}


void polyPrint(poly_t p){
	unsigned int i;
	printf("(");
	printRatio(p->coeffs[0]);
	printf(")");
	for(i=1; i<= p->deg; i++){
		if(p->coeffs[i].num != 0){
			printf("+ ");
			printf("(");
			printRatio(p->coeffs[i]);
			printf(")");
			printf("X^%u ", i);
		}else continue;
	}
	printf("\n");
}

poly_t polyDifferentiate(poly_t p){
	
	poly_t temp = polyFromMonomialRatio((p->coeffs)[p->deg],p->deg-1);
	temp->coeffs[p->deg-1] = mulRatio(temp->coeffs[p->deg-1], createRatio(p->deg,1));
	polyPrint(temp);
	unsigned int i;
	for(i=1; i < (p->deg); i++){
		temp->coeffs[i-1] = mulRatio(p->coeffs[i], createRatio(i,1));
		}
	return temp;
}

void polyGetLeadingCoeff(ratio_t *c, unsigned int *d, poly_t *r, poly_t p){
	*r = (struct __poly_struct_t*)malloc(sizeof(struct __poly_struct_t));
	*d =  polyGetDegree(p);
	*c = p->coeffs[*d];
	poly_t m = polyFromMonomialRatio(*c, *d);
	*r =  polySub(p, m);
}



void polyDiv(poly_t *q, poly_t *r, poly_t a, poly_t b){  // b est diviseur
	poly_t copieA = malloc(sizeof(struct __poly_struct_t*));
	memcpy(&copieA, &a, sizeof(copieA));

	*q = (struct __poly_struct_t*)malloc(sizeof(struct __poly_struct_t));
	(*q)->deg = a->deg - b->deg;
	(*q)->coeffs = malloc(sizeof(ratio_t)*(((*q)->deg)+1));

	memcpy(r, &a, sizeof(*r));

	unsigned int cpt = 0;
	while(!(polyIsConstantZero(*r)) && cpt <= (*q)->deg){
		ratio_t leadA, leadB;
		poly_t rA, rB;
		polyGetLeadingCoeff(&leadA, &((copieA)->deg), &rA, copieA);
		polyGetLeadingCoeff(&leadB, &(b->deg), &rB, b);

		ratio_t coeffTmp = divRatio(leadA, leadB);
		(*q)->coeffs[(*q)->deg - cpt] = coeffTmp;

		poly_t pm = polyFromMonomialRatio(coeffTmp, (*q)->deg - cpt);
		
		poly_t prod = polyMul(pm, b);

		*r = polySub(copieA, prod);
		copieA =polyFromCopy(*r);

		cpt++;
	}
}



int main(){

	ratio_t r[3] = {createRatio(3,1), createRatio(2,1), createRatio(1,2)};

	printf("\nTest polyFromRatioArray(ratio_t *c, unsigned int d): \n");
	poly_t p1 = polyFromRatioArray(r, 2); polyPrint(p1);
	
	printf("\n\nTest polyFromRatio(ratio_t c): \n");
	poly_t p2 = polyFromRatio(createRatio((int64_t)2,(int64_t)5)); polyPrint(p2);
	
	printf("\n\nTest polyFromInt(int64_t c): \n");
	poly_t p3 = polyFromInt(4); polyPrint(p3);
	poly_t p0 = polyFromInt(0); polyPrint(p0);
	
	printf("\n\nTest polyFromMonomialRatio(ratio_t c, unsigned int k): \n");
	poly_t p4 = polyFromMonomialRatio(r[2], 3); polyPrint(p4);
	
	printf("\n\nTest polyFromIdentity(): \n");
	poly_t p5 = polyFromIdentity(); polyPrint(p5);
	
	printf("\n\nTest polyGetDegree(poly_t p): \n");
	unsigned int deg = polyGetDegree(p1); printf("deg est %u\n", deg);
	
	printf("\n\nTest polyIsConstantZero(p3) et p0: \n");
	if(polyIsConstantZero(p3)) printf("Erreur p3 est zero\n");
	printf("p0 is constant 0 ? :");
	if(polyIsConstantZero(p0)) printf("P0 c'est 0\n");
	

	polyFree(p0); polyFree(p1); polyFree(p2); polyFree(p3); polyFree(p4); polyFree(p5);


	//Tests opérations sur polynomes

	ratio_t r1[5] = {createRatio(1,1), createRatio(2,1), createRatio(1,2), createRatio(1,1), createRatio(1,2)};
	ratio_t r2[3] = {createRatio(3,1), createRatio(2,1), createRatio(1,2)};

	printf("\n\nTests opérations sur polynomes\n");
	poly_t poly1 = polyFromRatioArray(r1, 4);
	poly_t poly2 = polyFromRatioArray(r2, 2); 
	polyPrint(poly1); polyPrint(poly2);

	printf("\n\nTest polyAdd(poly_t p, poly_t q)\n");
	polyPrint(polyAdd(poly1, poly2));
	
	printf("\n\nTest polySub(poly_t p, poly_t q)\n");
	polyPrint(polySub(poly1, poly2));
	
	printf("\n\nTest polyMul(poly_t p, poly_t q)\n");
	polyPrint(polyMul(poly1, poly2));
	
	printf("\n\nTest polyDifferentiate(poly_t p)\n");
	polyPrint(polyDifferentiate(poly1));
	
	printf("\n\nTest polyGetLeadingCoeff(ratio_t *c, unsigned int *d, reste, p1)\n");
	printf("p1: "); polyPrint(p1);
	ratio_t c;
	unsigned int d;
	poly_t reste;
	polyGetLeadingCoeff(&c, &d, &reste, p1);
	printf("c: "); printRatio(c);
	printf("   d: %u\n", d);
	printf("reste: "); polyPrint(reste); printf("\n\n");
	
	/* //Test pour la fonction polyDiv: ne marche pas
	printf("\n\nTest polyDiv(poly_t *q, poly_t *r, poly_t a, poly_t b)\n");
	poly_t pq, pr;
	printf("poly1/poly2 = \n\n");
	polyDiv(&pq, &pr, poly1, poly2);
	printf("pq: ");polyPrint(pq);
	printf("pr: ");polyPrint(pr);
	*/
	
	/* On a écrit la fonction polyDiv mais elle nous renvoie toujours l'erreur:

poly: malloc.c:2372: sysmalloc: Assertion `(old_top == (((mbinptr) (((char *) &((av)->bins[((1) - 1) * 2])) - __builtin_offsetof (struct malloc_chunk, fd)))) && old_size == 0) || ((unsigned long) (old_size) >= (unsigned long)((((__builtin_offsetof (struct malloc_chunk, fd_nextsize))+((2 *(sizeof(size_t))) - 1)) & ~((2 *(sizeof(size_t))) - 1))) && ((old_top)->size & 0x1) && ((unsigned long) old_end & pagemask) == 0)' failed.
Abandon

	*/


	polyFree(poly1); polyFree(poly2);// polyFree(reste); polyFree(pq); polyFree(pr); 
	
	return 0;
}

/* ---------------------------------------------------------- FONCTIONS DU TME3 ---------------------------------------------------------------------- */

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
	ext_eucl_div(&u, &v, &g, llabs(a), llabs(b));
	return g;
}

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

void printRatio(ratio_t a){
	printf("%" PRId64, a.num);
	printf("/%" PRId64, a.den);
}
