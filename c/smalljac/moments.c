#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "mpzutil.h"
#include "smalljac.h"
#include "smalljac_internal.h"
#include "cstd.h"
#include "polyparse.h"
#include "ff_poly.h"
#include "mpzpolyutil.h"

/*
    Copyright (c)  2010-2013 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// Program to demonstrate the use of the smalljac interface

#define MOMENTS		11
#define MAXM			4096
#define SPLIT_POLY_MAX_DEGREE	100

struct filter_param {
	struct {
		unsigned long mod;
		char filter[MAXM];
	} mod;
	struct {
		unsigned long a;
		unsigned long b;
	} ab;
	struct {
		int d;
		mpz_t poly[SPLIT_POLY_MAX_DEGREE+1];
		int partial;
	} splitpoly;
	struct {
		int factors;
		int pattern[SMALLJAC_MAX_DEGREE+1];
	} factors;
};

int filter_mod (smalljac_curve_t curve, unsigned long p, void *arg)
{
	struct filter_param *params = (struct filter_param *) arg;

	return params->mod.filter[p % params->mod.mod] == 1;
}

int filter_ab (smalljac_curve_t curve, unsigned long p, void *arg)
{
	ff_t t, inv;
	long k;
	struct filter_param *params = (struct filter_param *) arg;

	if ( params->mod.mod && ! params->mod.filter[p % params->mod.mod] ) return 0;
	ff_setup_ui(p);

	if ( (p-1) % 4 ) return 0;
	_ff_set_i(t, params->ab.a);
	k = (p-1) / 4;
	ff_exp_ui(&inv, &t, k);

	if ( params->ab.b ) {
		_ff_set_i(t, params->ab.b);
		if ( ! _ff_equal(inv, t) ) return 0;
	} else {
		if ( _ff_one(inv) ) return 0;
		ff_negate(inv);
		if ( _ff_one(inv) ) return 0;
	}

	return 1;
}

int filter_splitpoly (smalljac_curve_t curve, unsigned long p, void *arg)
{
	struct filter_param *params = (struct filter_param *) arg;

	ff_t f[SPLIT_POLY_MAX_DEGREE+1];
	int i;

	if ( params->mod.mod && ! params->mod.filter[p % params->mod.mod] ) return 0;
	ff_setup_ui(p);
	for ( i = 0 ; i <= params->splitpoly.d ; i++ )
		_ff_set_mpz (f[i], params->splitpoly.poly[i]);

	i = ff_poly_count_roots(f, params->splitpoly.d);
	return ( params->splitpoly.partial ? i : (i==params->splitpoly.d?1:0) );
/*	if ( ! i ) return 0;
	if ( ! params->splitpoly.partial && i != params->splitpoly.d ) return 0;
	
	for ( i = 0 ; i <= params->splitpoly.d ; i++ )
		_ff_set_mpz (f[i], params->splitpoly.poly[i]);

	if ( params->splitpoly.partial ) {
		return ff_poly_count_roots(f, params->splitpoly.d);
	} else {
		return 1;
	}
*/
}

int filter_factors (smalljac_curve_t curve, unsigned long p, void *arg)
{
	struct filter_param *params = (struct filter_param *) arg;

	int d;
	ff_t f[SPLIT_POLY_MAX_DEGREE+1], inv;
	int pattern[SMALLJAC_MAX_DEGREE+1];
	int i, j;

	if ( params->mod.mod && ! params->mod.filter[p % params->mod.mod] ) return 0;
	ff_setup_ui(p);
	if ( ! smalljac_Qcurve_reduce (f, curve) ) return 0;
	d = smalljac_curve_degree(curve);
	if ( ! _ff_one(f[d]) ) {
		_ff_invert(inv,f[d]);
		for ( j = 0 ; j <= d ; j++ ) ff_mult(f[j],f[j],inv);
	}

	if ( params->factors.factors > 0 ) {
		// skip prime if f doesn't split into specific number of factors
		if ( ff_poly_count_factors(f,d) != params->factors.factors ) return 0;
	} else {
		i = ff_poly_factorization_pattern (pattern, f, d);
		if ( i != -params->factors.factors ) return 0;
		for ( i = 1 ; i <= d ; i++ ) if ( pattern[i] != params->factors.pattern[i] ) break;
		if ( i <= d ) return 0;
	}
	return 1;
}

int main (int argc, char *argv[])
{
	//clock_t start_cpu;
	time_t start_time, end_time;
	unsigned long start, end;
	smalljac_curve_t curve;
	char *s, *t;
	int i, j, m, n, err, genus;
	char filter1[MAXM], filter2[MAXM], buf[256], *STgroup;
	double *moments;
	struct filter_param filter_param;
	int (*filter_func)(smalljac_curve_t, unsigned long, void *) = NULL;
	
	if ( argc < 4 ) {
		puts ("moments start end curve [n filter split]");
		puts ("    \"moments 1 1000000 x^5-3x^4+19x^3+4x^2+56x-12 0 11[1,3,4,5,9]\"\n    \"moments 1000 2000 [1,0,-1,37,42]\""); 
		puts ("The parmater n indicates that moments should be computed for L-poly coeffs a_1, a_2, ..., a_n.  By default n is set to the genus.");
		puts ("The filter \"11[1,3,4,5,9]\" restricts the scan to primes congruent to 1,3,4,5, or 9 mod 11.");
		puts ("The split parameter specifies a polynomial and restricts scan to primes p for which the polynomial factors completely mod p");
		puts ("You may specify a filename as the curve parameter to use precomputed lpoly data (in which case you get everything, not just a1).");
		return 0;
	}
	start = atol_exp (argv[1]);
	end = atol_exp (argv[2]);

	memset (&filter_param,0, sizeof(filter_param));
	if ( argc > 5 ) {
		for ( s = argv[5] ; *s && *s != '^' ; s++ );
		if ( *s == '^' ) {
			*s++='\0';
			filter_param.ab.a = atol(argv[6]);
			if ( filter_param.ab.a ) {
				filter_param.ab.b = atoi(s);
				if ( filter_param.ab.b ) {
					printf ("filtering for primes where %ld^((p-1)/4)=%ld\n", filter_param.ab.a, filter_param.ab.b);
				} else {
					printf ("filtering for primes where %ld^((p-1)/4)!=+/-1\n", filter_param.ab.a);
				}
			}
			filter_func = filter_ab;
		} else {
			for ( s = argv[5]+1 ; *s && *s != '[' ; s++ );
			if ( strstr(argv[5],"[") ) {
				n = atoi(argv[5]);
				if ( n <=1 || n > MAXM ) { printf ("Invalid filter modulus %d, MAXM = %d.\n", n, MAXM);  return 0; }
				filter_param.mod.mod = n;
				memset(&filter_param.mod.filter, 0, sizeof(filter_param.mod.filter));
				if ( ! *s ) { puts ("Invalid filter.");  return 0; }	
				do {
					n = atoi(++s);
					if ( !n || n >= filter_param.mod.mod ) { puts("Invalid filter.");  return 0; }
					filter_param.mod.filter[n] = 1;
					while ( *s && *s != ',' ) s++;
				} while (*s );
				printf ("filtering p mod %lu\n", filter_param.mod.mod);
				filter_func = filter_mod;
			} else if ( strstr(argv[5],"q") ) {
				s = argv[5];
				for ( s = argv[5] ; *s ; s++ ) {
					while ( *s && ! isdigit(*s) && *s != '-' ) s++;
					if ( ! *s ) break;
					for ( t = s ; *t == '-' || isdigit(*t) ; t++ );
					if ( *t != 'q' ) { printf ("Invalid quadratic residue filter %s\n", argv[5]); return 0; }
					n = atoi(s);  s = t;
					m = qr_mod (filter1, MAXM, n);
					if ( ! m )  { printf ("Exceeded MAXM = %d.\n", MAXM);  return 0; }
					if ( filter_param.mod.mod ) {
						memcpy (filter2, filter_param.mod.filter, filter_param.mod.mod);
						filter_param.mod.mod = qr_mod_merge (filter_param.mod.filter, MAXM, filter1, m, filter2, filter_param.mod.mod);
						if ( ! filter_param.mod.mod ) { printf ("Exceeded MAXM = %d\n", MAXM); return 0; }
					} else {
						memcpy (filter_param.mod.filter, filter1, m);  filter_param.mod.mod = m;
					}
					printf ("filtering for kron(%d,p) = 1\n", n);
				}
				filter_func = filter_mod;
			} else if ( strcmp(argv[5],"0") != 0 ) {
				printf ("Don't  know how to parse filter expression %s\n", argv[5]); return 127;
			}
		}
	}

	// Filtering with splitting polynomials and counting by degree
	if ( argc > 6 ) {
		if ( filter_func == filter_mod ) printf ("overriding kron/mod filtering by using splitpoly filtering\n");
		if ( strstr(argv[6],"x") || strstr(argv[6],"X") ) {
			filter_param.splitpoly.d = poly_parse_max_degree(argv[6]);
			if ( filter_param.splitpoly.d > 0 ) {
				if ( filter_param.splitpoly.d > SPLIT_POLY_MAX_DEGREE ) {
					printf ("Exceeded SPLIT_POLY_MAX_DEGREE = %d\n", SPLIT_POLY_MAX_DEGREE);
					return 0;
				}

				for ( i = 0 ; i <= filter_param.splitpoly.d ; i++ )
					mpz_init (filter_param.splitpoly.poly[i]);

				mpz_poly_parse (filter_param.splitpoly.poly, filter_param.splitpoly.d, argv[6]);
				if ( argv[6][strlen(argv[6])-1] == '+' ) {
					filter_param.splitpoly.partial = 1;
					printf ("filtering for primes that split into at least one degree-one prime in the number field defined by %s\n", argv[6]); 
				} else {
					printf ("filtering for primes that split completely in the number field defined by %s\n", argv[6]);
				}
			} else {
				filter_param.splitpoly.d = 0;
			}

			filter_func = filter_splitpoly;
		} else if ( argv[6] [0]== '[' ) {
				filter_param.factors.factors = 0;
				for ( i = 0 ; i <= SMALLJAC_MAX_DEGREE ; i++ ) filter_param.factors.pattern[i] = 0;
				for ( s = argv[6]+1 ; *s ; ) {
					if ( ! isdigit(*s) ) { printf ("Invalid factorization pattern %s specified\n", argv[6]); return 0; }
					i = atoi (s);
					if ( i > SMALLJAC_MAX_DEGREE )  { printf ("Invalid factorization pattern %s specified\n", argv[6]); return 0; }
					filter_param.factors.pattern[i]++;
					filter_param.factors.factors--;
					while ( isdigit(*s) ) s++;
					while ( *s && ! isdigit(*s) ) s++;
				}
				if ( ! filter_param.factors.factors )  { printf ("Invalid factorization pattern %s specified\n", argv[6]); return 0; }
				printf ("filtering for primes with factorization pattern %s\n", argv[6]);
				filter_func = filter_factors;
		} else {
			filter_param.factors.factors = atoi(argv[6]);
			if ( filter_param.factors.factors ) {
				printf ("filtering for primes where f(x) splits into %d factors\n", filter_param.factors.factors);
			}
			filter_func = filter_factors;
		}
	}

	/*
	// if curve parameter contains the string "data" assume it is a file of Lpoly data, not a curve specification
	if ( strstr(argv[3],"data") ) {
		printf ("Reading Lpoly data from file %s\n", argv[3]);
		start_time = time(0);
		result = smalljac_Lpolys_from_file (argv[3], start, end, flags, moments, (void*)&context);
		if ( result < 0 ) {  printf ("smalljac_Lpolys_from_file returned error %ld\n", result);  return 0; }
		context.cpu_musecs = clock()-start_cpu;  end_time = time(0);
		goto done;
	}*/
	
	// create the curve - give an explanation of the error if there is a problem
	curve = smalljac_curve_init (argv[3], &err);
	if ( ! curve ) { 
		switch (err) {
		case SMALLJAC_PARSE_ERROR: printf ("Unable to parse curve string: %s\n", argv[3]); break;
		case SMALLJAC_UNSUPPORTED_CURVE: puts ("Specified curve not supported, check equation\n");  break;
		case SMALLJAC_SINGULAR_CURVE: puts ("Specified curve is singular\n");  break;
		default: printf ("smalljac_curve_init returned error %d\n", err);
		}
		return 0;
	}

	n = genus = smalljac_curve_genus(curve);
	if ( argc > 4 ) {
		i = atoi(argv[4]);
		if ( i > n ) { printf ("n cannot exceed the genus, which is %d\n", n); return 0; }
		if ( i ) n = i;
	}
	m = MOMENTS;

	moments = (double *)malloc(sizeof(double) * n * m);

	if ( n == genus && genus < 3 ) { STgroup = buf; 	STgroup[0] = '\0'; } else { STgroup = 0; }

	start_time = time(0);
	smalljac_moments(curve, start, end, moments, n, m, STgroup, filter_func, &filter_param);
	end_time = time(0);

	for ( i = 0 ; i < n ; i++ ) {
		printf("a%d Moments: 1,  ", i + 1);
		for ( j = 1 ; j < m ; j++ ) {
			printf("%.3f,  ", moments[i * m + j]);
		}
		printf("\n");
	}

	if ( STgroup ) {
		if ( strlen(STgroup) ) {
			printf ("ST group provisionally identified as %s\n", STgroup);
		} else {
			printf ("Unable to identify ST group, consider increasing prime bound\n");
		}
	}

	printf("Total elapsed time: %ld seconds\n", end_time - start_time);
}
