#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "smalljac.h"
#include "cstd.h"

/*
    Copyright (c)  2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// simple command-line program to compute the L-polynomial of a J(C/F_p)

int main (int argc, char *argv[])
{
	clock_t start,end;
	unsigned long p;
	long a[3];
	int sts, flags;

	if ( argc < 3 ) { puts ("lpoly curve q [a1-only]");   puts ("   \"lpoly [1,2,3,4,5] 65537\"");    puts("    \"lpoly x^5+3x^3-2x+1 10000019\""); return 0; }

	p = atoll(argv[2]);
	flags = ( argc > 3 ? atoi(argv[3]) : 0 );
	if ( flags > 0 ) flags = SMALLJAC_A1_ONLY;

	start = clock();
	sts = smalljac_Lpoly (a, argv[1], p, flags);
	end = clock();
	if ( sts <= 0 ) { 
		switch (sts) {
		case 0: printf ("Specified curve has bad reduction at %lu\n", p);  break;
		case SMALLJAC_PARSE_ERROR: printf ("Unable to parse curve string: %s\n", argv[1]); break;
		case SMALLJAC_UNSUPPORTED_CURVE: puts ("Specified curve not supported\n");  break;
		case SMALLJAC_SINGULAR_CURVE: puts ("Specified curve is singular\n");  break;
		case SMALLJAC_INVALID_PP: puts ("Invalid modulus (not prime or out of supported range)");  break;
		default: printf ("smalljac_Lpoly returned error %d\n", sts);
		}
		return 0;
	}
	if ( flags&SMALLJAC_A1_ONLY ) {
		printf ("a1 coefficient of L_p(T) is %ld\n", a[0]);
	} else {
		switch ( sts ) {
		case 1: printf ("L_q(T) = %ld*T^2 + %ld*T + 1\n", p, a[0]);  break;
		case 2: printf ("L_q(T) = %ld*T^4 + %ld*T^3 + %ld*T^2 + %ld*T + 1\n", p*p, p*a[0], a[1], a[0]);  break;
		case 3: printf ("L_q(T) = %ld*T^6 + %ld*T^5 + %ld*T^4 + %ld*T^3 + %ld*T^2 + %ld*T + 1\n", p*p*p, p*p*a[0], p*a[1], a[2], a[1], a[0]);  break;
		default: printf ("Lpoly returned unexpected value %d\n", sts);
		}
	}
	if ( end > start ) printf ("%.3f secs\n", (double)(end-start)/CLOCKS_PER_SEC);
	return 0;
}
