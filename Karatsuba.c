# include <stdio.h>
# include <stdlib.h>
# include <stdint.h>
# include <inttypes.h>
# include <math.h>
# include <time.h>

/* Par hypothèse, on prend 0 <= a,b <= 1 */

uint64_t multiplication1Bit(uint64_t a, uint64_t b){

	return a&b;
/*	if(a==1&&b==1) return (uint64_t)1;
	return (uint64_t)0; */
}

/* Il s'agit de coder Karatsuba avec les 4 appels récursifs de taille n/2, complexité en n² */

uint64_t multiplicationNaiveAux(uint64_t a, uint64_t b, unsigned int n){

	if(n==1) return multiplication1Bit(a,b);
	uint64_t x = a, y = b;
	uint64_t xh, xl, yh, yl;
	int nhalf = n>>1; //n/2
	xh = x>>nhalf;
	xl = a - (xh<<nhalf);
	yh = y>>nhalf;
	yl = b - (yh<<nhalf);
	
	return ((multiplicationNaiveAux(xh,yh,nhalf))<<n) 
		+ ((multiplicationNaiveAux(xh,yl,nhalf)+multiplicationNaiveAux(xl,yh,nhalf))<<(nhalf)) 
		+ multiplicationNaiveAux(xl,yl,nhalf);
}


/* Maintenant codé sur 32 bits */

uint64_t multiplicationNaive(uint32_t a, uint32_t b){

	return multiplicationNaiveAux((uint64_t) a, (uint64_t) b, 32u);
}


/* Karatsuba avec 3 appels récursifs de taille n/2, complexité en n*log2(3) */

uint64_t multiplicationKaratsubaAux(uint64_t a, uint64_t b, unsigned int n){

	if(n==1) return multiplication1Bit(a,b);
	uint64_t x = a, y = b;
	uint64_t xh, xl, yh, yl, xhtmp, yhtmp;
	unsigned int n0 = n, n2 = n;
	unsigned int nhalf = n0>>1; //n/2

	xh = x>>nhalf;
	xhtmp = xh; //On garde la valeur de xh et yh avant décalage
	xl = a - (xh<<nhalf);
	yh = y>>nhalf;
	yhtmp = yh;
	yl = b - (yh<<nhalf);

	uint64_t t1, t2, t1h, t1l, t2h, t2l;
	uint64_t t1htmp, t2htmp;
	t1 = xhtmp+xl; //sur (n/2)+1 bits
	t1h = t1>>(nhalf); //retenue 1 bit
	t1htmp=t1h;
	t1l = t1- (t1h<<(nhalf)); //reste sur n/2 bits
	
	t2 = yhtmp+yl; //sur (n/2)+1 bits
	t2h = t2>>(nhalf); //retenue 1 bit
	t2htmp=t2h;
	t2l = t2- (t2h<<(nhalf)); //reste sur n/2 bits
	
	// (x*y) = (xh.2^(n/2) + xl)(yh.2^(n/2) + yl) = xh.yh.2^n + (xhyl+xlyh)* 2^(n/2) + xl.yl
	//Donc (x*y) = xh.yh*(2^n -2^(n/2)) + xl.yl(1-2^(n/2)) + (xh+xl)(yh+yl)*2^(n/2) 

	// (xh+xl)(yh+yl) = t1h.t2h.2^n + t1h.t2l.2^(n/2) + t1l.t2h.2^(n/2) + t1l.t2l
	uint64_t res = multiplicationKaratsubaAux(t1l,t2l,nhalf);
	
	if(t1htmp) res += (t2l<<nhalf);
	if(t2htmp) res += (t1l<<nhalf);
	res += (multiplication1Bit(t1htmp,t2htmp))<<n2 ; //res = (xh+xl)(yh+yl)
	
	uint64_t res1, res2, res3; //x.y = 2^n.res1+ 2^(n/2).res2 + res3
	res1 = multiplicationKaratsubaAux(xhtmp,yhtmp,nhalf);
	res3 = multiplicationKaratsubaAux(xl,yl,nhalf);
	res2 = res - res1 - res3;
	
	return (res1<<n2) + (res2<<(nhalf)) + res3;

}

uint64_t multiplicationKaratsuba(uint32_t a, uint32_t b){

	return multiplicationKaratsubaAux((uint64_t) a, (uint64_t) b, 32u);
}

uint64_t multiplicationMachine(uint32_t a, uint32_t b){

	return ((uint64_t) a)*((uint64_t) b);
}


uint32_t random32bits(){
	
	uint32_t r1, r2;
	r1 = lrand48();
	r1 = r1&(0xffff);
	r2 = lrand48();
	r2 = r2&(0xffff);
	return (r1<<16)|r2;	
}




int main(){
	srand48(time(NULL));
	int i;
	for(i = 0; i < 300; i++){
		
		uint32_t a = random32bits(), b = random32bits();

		if((multiplicationNaive(a,b)==multiplicationKaratsuba(a,b))&&(multiplicationKaratsuba(a,b)==multiplicationMachine(a,b))){
			printf("true\n");
		} else {
			printf("false\n");
			break;
		}
	}
	return 0;

}
