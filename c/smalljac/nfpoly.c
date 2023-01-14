#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>
#include "ff_poly.h"
#include "nfpoly.h"
#include "mpzpolyutil.h"
#include "cstd.h"

/*
    Copyright (c)  2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

int nf_poly_init (nf_poly *nf, char *expr)
{
	char buf[4096];
	char *s, *t, x, z, *c[NFPOLY_MAX_DEGREE+1];
	int i, j, d, sign;

	if ( sizeof(ff2k_t) != sizeof(ff2k_t) ) { err_printf ("Fatal error in nf_poly_init: types ff2k_t and ff_t differ in size\n"); abort(); }

	if ( strlen(expr)+1 > sizeof(buf) ) { err_printf ("number field polynomial expression too large to parse in nf_poly_init\n"); abort(); }
	for ( t = buf, s = expr ; *s ; s++ ) if ( ! isspace(*s) ) *t++ = *s;			// strip whitespace to simplify parsing
	
	// make sure expr looks like [F(x))]/(G(z)) (where the coefficients of F are polys in Z[z])
	if ( buf[0] != '[' ) return -2;
	for ( s = buf+1 ; *s && *s != ']' ; s++ );
	if ( ! *s ) return -2;
	if ( *(++s) != '/' || *(++s) != '(' ) return -2;
	// find the variable name z in G (should be first character since G is monic) and the variable name x in F
	z = *(++s); if ( ! isalpha(z) ) return -2;

	// assume the first term in G is the highest degree term and get its degree
	if ( *(s+1) == '^' ) d = atoi(s+2); else d = 1;
	memset (nf, 0, sizeof(*nf));
	mpz_init (nf->Disc);  mpz_init (nf->D);
	nf->G = mem_alloc ((d+1)*sizeof(*nf->G));
	nf->g = mem_alloc ((d+1)*sizeof(*nf->g));
	for ( i = 0 ; i <= d ; i++ ) mpz_init (nf->G[i]);

	// parse the poly G(z)
	for ( t = s ; *t && *t != ')' ; t++ );
	if ( *t != ')' ) { nf_poly_clear (nf); return -2; }
	*t = '\0';
	nf->n = mpz_poly_parse (nf->G, d, s);
	if ( nf->n != d ) { nf_poly_clear (nf); return -2; }
	
	// for quadratic G, compute the discriminant of N and make it maximal at 2
	if ( nf->n == 2 ) {
		mpz_mul (nf->Disc, nf->G[1], nf->G[1]);
		mpz_mul_2exp (nf->D, nf->G[0], 2);
		mpz_sub (nf->Disc, nf->Disc, nf->D);
		mpz_set (nf->D, nf->Disc);
		while ( mpz_divisible_ui_p (nf->D, 4) ) mpz_div_2exp (nf->D, nf->D, 2);
		if ( ! mpz_congruent_ui_p (nf->D, 1, 4) ) mpz_mul_2exp (nf->D, nf->D, 2);
	} else { // otherwise just compute the discriminant the hard way
        mpz_t *w;
        w = mem_alloc ((d*d+2)*sizeof(mpz_t));
        for ( i = 0 ; i < d*d+2 ; i++ ) mpz_init (w[i]);
        mpz_poly_discriminant (nf->Disc, nf->G, nf->n, w);
        for ( i = 0 ; i < d*d+2 ; i++ ) mpz_clear (w[i]);
        mem_free (w);
    }

	// determine whether we think F is in Weierstrass form (if it contains a comma, we assume it is), if not find the variable name
	if ( strchr (buf, ',') ) {
		nf->WSflag = 1;
		x = ',';
	} else {
		nf->WSflag = 0;
		for ( t = buf+1; *t && (*t == z || ! isalpha(*t)) ; t++ );
		if ( ! *t ) { nf_poly_clear (nf); return -2; }
		x = *t;
	}
    
    // as a sanity check, make sure that no alpha characters other than those we expect appear
    if ( nf->WSflag ) {
        for ( s = buf ; *s ; s++ ) if ( isalpha(*s) && *s != z ) { nf_poly_clear(nf); return -2; }
    } else {
        for ( s = buf ; *s ; s++ ) if ( isalpha(*s) && *s != z && *s != x ) { nf_poly_clear(nf); return -2; }
    }
    
	// now break up the expr for F(x) into coeffiicient polys, computing the degrees and offsets of each
	s = buf+1;
	if ( *s == ']' ) { nf_poly_clear (nf); return -2; }
	for ( i = 0 ; i <= NFPOLY_MAX_DEGREE ; i++ ) { nf->cd[i] = -1; c[i] = 0; }

	i = -1;
	while ( *s != ']' && *s != '/' ) {
		for ( t = s ; *t && *t != x && *t != ']' ; t++ );
		if ( ! *t ) { nf_poly_clear (nf); return -2; }
		// determine the degree i of the current term f_i(z)*x^i, for WS form, degree is simply an index 0,1,2,3,4 into the coefficient array
		if ( nf->WSflag ) {
			++i;
		} else {
			if ( *t == ']' ) i = 0;  else if ( *(t+1) == '^' ) i = atoi(t+2);  else i = 1;
		}
		*t = '\0';
		if ( *(t-1) == '*' ) *(t-1) = '\0';
		if ( c[i] ) { nf_poly_clear (nf); return -2; }			// don't deal with repeated coefficients
		c[i] = s;
		// determine the degree of f_i(z)
		if ( *c[i] == '0' && ! *(c[i]+1) ) {					// try to catch zeros in WS form (otherwise, there is no reason they should be there)
			nf->cd[i] = -1;
		} else {
			for ( d = 0 ; *s ; s++ ) {
				if  ( *s == z ) if ( ! d ) d = 1;
				if ( *s == '^' ) { j = atoi(s+1); if ( j > d ) d = j; }
			}
			nf->cd[i] = d;
		}
		s = t+1;
		if ( *s == '^' ) for ( s++ ; isdigit(*s) ; s++ );
	}
	
	// determine the degree of f and validate WS form
	if ( ! nf->WSflag ) for ( i = NFPOLY_MAX_DEGREE ; nf->cd[i] == -1 ; i-- );
	else if ( nf->cd[0]==-1 && nf->cd[1]==-1 && nf->cd[2]==-1 ) { c[0] = c[3]; nf->cd[0] = nf->cd[3]; c[1] = c[4]; nf->cd[1] = nf->cd[4]; i = 1; }
	nf->d = i;
	if ( nf->WSflag && nf->d != 1 && nf->d !=4 ) { err_printf ("Coefficient list in Weierstrass form must have length 2 or 5\n"); nf_poly_clear (nf); return -2; }
	
	// determine the total number of integer coefficients, and offsets of f_i coefficients
	for ( i = 0, nf->cc = 0 ; i <= nf->d ; i++ ) nf->cc += nf->cd[i]+1;
	nf->co[0] = 0;
	for ( i = 1 ; i <= nf->d ; i++ ) nf->co[i] = nf->co[i-1]+nf->cd[i-1]+1;
	if ( ! nf->cc ) return -1;
	nf->F = mem_alloc (nf->cc*sizeof(*nf->F));
	for ( i = 0 ; i < nf->cc ; i++ ) mpz_init (nf->F[i]);
	nf->f = mem_alloc (nf->cc*sizeof(nf->f));
	for ( i = 0 ; i <= nf->d ; i++ ) {
		if ( nf->cd[i] >= 0 ) {
			s = c[i];
			sign = 1;
			if ( *s == '-'  &&  (*(s+1) == '(' || !*(s+1)) ) { sign = -1; s++; }
			else if ( *s == '+' && (*(s+1) == '(' || !*(s+1))  ) s++;
			if ( ! *s ) mpz_set_ui (nf->F[nf->co[i]], 1);			// empty string should be treated as 1 here
			else mpz_poly_parse (nf->F + nf->co[i], nf->cd[i], s);
			if ( sign < 0 ) for ( j = 0 ; j <= nf->cd[i] ; j++ ) mpz_neg (nf->F[nf->co[i]+j], nf->F[nf->co[i]+j]);
//printf ("coeff %d, expr %s deg %d :", i, s, nf->cd[i]); for  ( k = 0 ; k <= nf->cd[i] ; k++ ) gmp_printf ("%Zdz^%d ", nf->F[nf->co[i]+k], k); puts ("");
		}
	}
	if ( nf->d == 3 )
		if ( nf->cd[3] != 0 || nf->cd[2] != -1 || mpz_cmp_ui(nf->F[nf->co[3]], 1) != 0 ) { err_printf ("Cubics should be in the form x^3+A*x+B, or use Weierstrass [a1,a2,a3,a4,a6]\n"); nf_poly_clear (nf); return -2; }
	
	for ( i = 0 ; i <= nf->d  && nf->cd[i] <= 0 ; i++ );
	nf->Qflag = ( i > nf->d ? 1 : 0 );
	nf->r = mem_alloc (nf->n*sizeof(*nf->r));
	nf->k = 0;  nf->p = 0;
//if ( nf->WSflag ) printf ("nf_poly initialized genus 1 Weierstrass equation over number field of degree %d, Qflag = %d\n", nf->n, nf->Qflag);
//else printf ("nf_poly initialized degree %d poly over number field of degree %d, Qflag = %d\n", nf->d, nf->n, nf->Qflag);
	return ( nf->WSflag ? 3 : nf->d );
}


void nf_poly_clear (nf_poly *nf)
{
	int i;

	mpz_clear (nf->Disc); mpz_clear (nf->D);
	if ( nf->F ) { for ( i = 0 ; i < nf->cc ; i++ ) mpz_clear (nf->F[i]); mem_free (nf->F); }
	if ( nf->f ) mem_free (nf->f);
	if ( nf->G ) { for ( i = 0 ; i <= nf->n ; i++ ) mpz_clear (nf->G[i]); mem_free (nf->G); }
	if ( nf->g ) mem_free (nf->g);
	if ( nf->r ) mem_free (nf->r);
}

int nf_poly_reduce_setup_char2 (nf_poly *nf, int e)
{
	int i, k;
	
	ff2k_setup (e);
	// optimized code for quadratic fields
	if ( nf->n == 2 ) {
		i = mpz_kronecker_ui (nf->D, 2);
		if ( e == 1 ) {
			nf->k = i+1;
			if ( nf->k > 0 ) {
				ff2k_poly_set_mpz (nf->g, nf->G, nf->n);
				k = ff2k_poly_distinct_roots (nf->r, nf->g, nf->n);
				if ( k != nf->k ) { err_printf ("unexpected number of distinct roots of defining poly for quadratic field in char 2\n"); abort(); }
				if ( nf->k == 2 && _ff2k_zero(nf->r[1]) ) { _ff2k_set (nf->r[1], nf->r[0]); _ff2k_set_zero (nf->r[0]); }	// sort images of degree-1 prime ideals in [0,1] (i.e. flip 1,0 to 0,1 if necessary).
			}
		} else {
			nf->k = ( i < 0 ? 1 : 0 );
			if ( nf->k ) {
				ff2k_poly_set_mpz (nf->g, nf->G, nf->n);
				k = ff2k_poly_distinct_roots (nf->r, nf->g, nf->n);
				if ( ! k ) { err_printf ("Unable to obtain a root of g in F_4\n"); abort(); }
			}
		}
	} else {
		if ( e == 1 ) {
			ff2k_poly_set_mpz (nf->g, nf->G, nf->n);
			nf->k = ff2k_poly_distinct_roots (nf->r, nf->g, nf->n);
			if ( nf->k == 2 && _ff2k_zero(nf->r[1]) ) { _ff2k_set (nf->r[1], nf->r[0]); _ff2k_set_zero (nf->r[0]); }		// sort images of degree-1 prime ideals in [0,1] (i.e. flip 1,0 to 0,1 if necessary).
		} else {
			err_printf ("Support for degree-%d primes in degree-%d number fields not implemented in nfpoly.c\n", e, nf->n);  abort();
		}
	}
	if ( nf->k ) ff2k_poly_set_mpz (nf->f, nf->F, nf->cc-1);
	nf->p = 2;  nf->e = e;
	return nf->k;	
}


// setup reduction for degree e primes above _ff_p (sets the current finite field if it is not already set)
// currently e must be 1 unless e=nf->n=2
int nf_poly_reduce_setup (nf_poly *nf, long p, int e)
{
	ff_t w[2];
	
    if ( mpz_divisible_ui_p (nf->Disc, p) ) return 0;
	if ( p == 2 ) return nf_poly_reduce_setup_char2 (nf, e);
	if ( p != _ff_p ) ff_setup_ui (p);
	if ( e == 1 ) {
		ff_poly_set_mpz (nf->g, nf->G, nf->n);
		nf->k = ff_poly_distinct_roots (nf->r, nf->g, nf->n);
		ff_sort (nf->r, nf->k); 		// sort the roots according to their image in Fp -> Z/pZ -> [0..(p-1)]
	} else if ( e == 2 && nf->n == 2 ) {
		_ff_set_mpz (w[0], nf->Disc);  _ff_set_zero(w[1]);
		if ( ! ff2_sqrt (w, w) ) { err_printf ("Unable to compute sqrt of an element of Fp in Fp^2 !!!\n"); abort(); }
		if ( _ff_zero(w[1]) ) {
			nf->k = 0;
		} else {
			_ff_set_mpz (nf->r[0], nf->G[1]); ff_negate(nf->r[0]); _ff_set_zero(nf->r[1]);
			ff2_add (nf->r, nf->r, w); ff2_scalar_mult (nf->r, _ff_half, nf->r);
			nf->k = 1;
		}
	} else {
		 err_printf ("Support for degree-%d primes in degree-%d number fields not implemented in nfpoly.c\n", e, nf->n);  abort();
	}
	if ( nf->k ) ff_poly_set_mpz (nf->f, nf->F, nf->cc-1);
	nf->p = p;  nf->e = e;
	return nf->k;
}

int nf_poly_reduce_char2 (ff_t w[], nf_poly *nf, int e, int j)
{
	int i;
	
	if ( ! nf->WSflag || nf->d != 4 ) return 0;
	for ( i = 0 ; i <= nf->d ; i++ ) ff2k_poly_eval (w+i, nf->f+nf->co[i], nf->cd[i], nf->r+j);
	return nf->d;
}

// returns the reduction of F at the j-th degree d prime,  f should have space for d*(nf->d+1) elements
int nf_poly_reduce (ff_t f[], nf_poly *nf, long p, int e, int j)
{
	int i;
	
//printf ("nf_poly reduce p=%ld, e=%d\n", p, e);
	if ( p != nf->p || nf->e != e ) { err_printf ("You must call nf_poly_reduce_setup for degree %d primes above p=%ld before calling nf_poly_reduce\n", e, p); exit (0); }
	if ( p > 2 && p != _ff_p ) { err_printf ("Finite field change since call to nf_poly_reduce_setup, you can't switch fields between calls (p=%ld, _ff_p=%ld\n", p, _ff_p); exit (0); }
	if ( p == 2 ) return nf_poly_reduce_char2 (f, nf, e, j);
	if ( j < 0 || j >= nf->k ) { err_printf ("Invalid prime index %d\n", j);  abort(); }
	if ( e == 1 ) {
//printf ("r=%ld\n", _ff_get_ui(nf->r[0])); 
		if ( nf->WSflag ) {
			if ( nf->d == 1 ) {
				// note that [a4,a6] --> [f1,f0]
				_ff_set_one (f[3]); _ff_set_zero (f[2]);  ff_poly_eval (f+1, nf->f+nf->co[0], nf->cd[0], nf->r+j);  ff_poly_eval (f, nf->f+nf->co[1], nf->cd[1], nf->r+j);
			} else {
				ff_t w[5];
				if ( nf->d != 4 ) { err_printf ("Expected 5 Weierstrass coefficients in nf_poly_reduce, only found %d\n", nf->d+1); abort(); }
				for ( i = 0 ; i <= nf->d ; i++ ) ff_poly_eval (w+i, nf->f+nf->co[i], nf->cd[i], nf->r+j);
//printf ("[%ld,%ld,%ld,%ld,%ld]\n", _ff_get_ui(w[0]), _ff_get_ui(w[1]), _ff_get_ui(w[2]), _ff_get_ui(w[3]), _ff_get_ui(w[4]));

				if ( _ff_p == 3 ) ff_poly_med_weierstrass (f, w);
				else ff_poly_short_weierstrass (f, w);
			}
			i = 3;
		} else {
			for ( i = 0 ; i <= nf->d ; i++ ) 	ff_poly_eval (f+i, nf->f+nf->co[i], nf->cd[i], nf->r+j);
			for ( i = nf->d ; i >= 0 && _ff_zero (f[i]) ; i-- );
		}
//printf ("%ld(%d): ", _ff_p, j); ff_poly_print (f,i);
		return i;
	} else if ( e== 2 && nf->n == 2 ) {
//printf ("r=%ld*z+%ld\n", _ff_get_ui(nf->r[1]), _ff_get_ui(nf->r[0])); 
		if ( nf->WSflag ) {
			if ( nf->d == 1 ) {
				ff2_set_one (f+6); ff2_set_zero (f+4);  ff2_poly_eval_ff (f+2, nf->f+nf->co[0], nf->cd[0], nf->r+j);  ff2_poly_eval_ff (f, nf->f+nf->co[1], nf->cd[1], nf->r+j);
			} else {
				ff_t w[10];
				for ( i = 0 ; i <= nf->d ; i++ ) ff2_poly_eval_ff (w+e*i, nf->f+nf->co[i], nf->cd[i], nf->r+e*j);
				if ( _ff_p == 3 ) ff2_poly_med_weierstrass (f, w);
				else ff2_poly_short_weierstrass (f, w);
			}
			i = 3;
		} else {
			for ( i = 0 ; i <= nf->d ; i++ )  ff2_poly_eval_ff (f+e*i, nf->f+nf->co[i], nf->cd[i], nf->r+e*j);
			for ( i = nf->d ; i >= 0 && ff2_zero (f+2*i) ; i-- );
		}
		return i;
	}
	 err_printf ("Support for degree-%d primes in nfpoly.c not supported for degree-%d number fields\n", e, nf->n);  abort();
}

#ifdef TEST_NFPOLY

int main (int argc, char *argv[])
{
	nf_poly nf[1];
	ff_t f[2*NFPOLY_MAX_DEGREE+2];
	long p;
	int i, k;
	
	if ( argc < 2 ) { puts("Please specify poly over a number field in the form [f(x)]/(g(z)), where f is in (Z[z])[x] and g is in Z[z]"); return 0; }
	if ( nf_poly_init (nf, argv[1]) == -2 ) { printf ("nf_poly_init failed, likely due to a syntax error\n"); return 0; }
	printf ("Degree of number field is %d, degree of f(x) is %d\n", nf->n, nf->d);
	if ( argc > 2 ) {
		p = atol (argv[2]);
		ff_setup_ui (p);
		k = nf_poly_reduce_setup (nf, p, 1);
		printf ("%d degree 1 primes above p=%ld\n", k, p);
		for ( i = 0 ; i < k ; i++ ) {
			nf_poly_reduce (f, nf, p, 1, i);
			printf ("Reduction at %d of %d degree 1 primes above p=%ld is: ", i+1, k, p);  ff_poly_print(f,nf->d); puts ("");
		}
		if ( k == 0 && nf->n == 2 ) {
			k = nf_poly_reduce_setup (nf, p, 2);
			printf ("%d degree 2 primes above p=%ld\n", k, p);
			for ( i = 0 ; i < k ; i++ ) {
				nf_poly_reduce (f, nf, p, 2, i);
				printf ("Reduction at %d of %d degree 1 primes above p=%ld is: ", i+1, k, p);  ff2_poly_print(f,nf->d); puts ("");
			}
		}
	}
	nf_poly_clear (nf);
}

#endif
