#ifndef _NFPOLY_INCLUDE_
#define _NFPOLY_INCLUDE_

/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff_poly.h"

#define NFPOLY_MAX_DEGREE   15

#ifdef __cplusplus
extern "C" {
#endif

// for convenience we use the same variables in odd prime fields and prime fields of order 2^k (for k<32), it is assumed that types ff2k_t and ff_t are the same
// a size check is performed in nf_poly_init

// data structure to support polynomials over number fields and their reduction to finite fields
typedef struct nf_poly_struct {
	mpz_t *F;												// array of coefficients of F in (Q[z])[x].  coefficients are stored in order of increasing degree in x and then z
	mpz_t *G;												// array of n+1 coefficients of the poly defining the number field Q[z]/(G(z)), must be monic and irreducible (not verified)
	int d;													// degree of F (or 5 if WSflag is set)
	int n;													// degree of G (at most POLY_MAX_DEGREE)
	int cc;													// total number of integer coefficients stored in F
	int cd[NFPOLY_MAX_DEGREE+1];							// cd[i] is the degree of x^i in F(x), zero coefficients have degree -1
	int co[NFPOLY_MAX_DEGREE+1];							// co[i] is the offset of the coefficient poly for x^i in F(x)
	ff_t *g;												// reduction of G mod p
	ff_t *f;												// reduction of F in Fp^k (where 1 <= k <= n)
	ff_t *r;												// roots of G in Fp, each corresponds to a degree-1 prime ideal, ordered via the map Fp->Z/pZ->[0,p-1]
	int k;													// number of roots of G in Fp
	int e;													// degree of primes for which reduction has been setup
	long p;													// prime at which reduction has been setup
	mpz_t Disc;												// Disc is the discriminant of G
	mpz_t D;												// for n=2, D = Disc made maximal at 2 (i.e., powers of 4 removed as long as D=0,1 mod 4)
	int Qflag;												// true if the coefficients of F all lie in Q
	int WSflag;												// F is actually in Weierstrass form [a1,a3,a2,a4,a6].
} nf_poly;

int nf_poly_init (nf_poly *nfp, char *expr);
void nf_poly_clear (nf_poly *nfp);
int nf_poly_reduce_setup (nf_poly *nfp, long p, int e);     // returns 0 if p divides Disc, otherwise returns the number of degree e primes above p
int nf_poly_reduce (ff_t f[], nf_poly *nfp, long p, int e, int j);

#ifdef __cplusplus
}
#endif

#endif
