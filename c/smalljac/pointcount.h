#ifndef _POINTCOUNT_INCLUDE_
#define _POINTCOUNT_INCLUDE_

/*
    Copyright (c) 2007-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/


#include <gmp.h>
#include "ff_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
	This module implements pointcounting for hyperelliptic curves over finite fields of odd
	characteristic based on finite differences, as described in [KedlayaSutherland2007].
	
	Genus 1 curves y^2=f(x) where deg(f) is 3 or 4 are also supported (strictly speaking
	these are elliptic curves, not hyperelliptic curves).
	
	None of these functions checks for good reduction, however,
	GOOD REDUCTION IS ASSUMED...
*/

#define POINTCOUNT_MULTI_X		32
#define POINTCOUNT_MAX_TINYP	(1<<8)

// Allocates memory for map of residues - must call first
void pointcount_init (unsigned maxp);

// Note that precompute returns integer values, independent of p
// These need to be reduced mod p prior to calling any of the functions below
void pointcount_precompute (mpz_t D[], mpz_t f[], int degree);
void pointcount_precompute_long (long D[], long f[], int degree);

// handy inline for doing reduction (note sign handling, d[i] must have nonnegative values)
static inline void pointcount_reduce (unsigned long d[], long D[], int degree, long p)
	{ register long x; register int i; for ( i = 0 ; i <= degree ; i++ ) { x=D[i]%p;  d[i] = (x < 0 ? x+p : x); } }

// Standard hyperelliptic point counting functions for curves of the form y^2=f(x)  - default is degree 2g+1
unsigned pointcount_g1 (unsigned long D[4], unsigned p);
unsigned pointcount_g1d4 (unsigned long D[5], unsigned p, unsigned long f4);
unsigned pointcount_g2 (unsigned long D[6], unsigned p);
unsigned pointcount_g2d6 (unsigned long D[7], unsigned p, unsigned long f6);
unsigned pointcount_g3 (unsigned long D[8], unsigned p);
unsigned pointcount_g3d8 (unsigned long D[9], unsigned p, unsigned long f8);
unsigned pointcount_g4 (unsigned long D[10], unsigned p);
unsigned pointcount_g4d10 (unsigned long D[10], unsigned p, unsigned long f10);

// Point counting over Picard curves y^3 = f(x)
unsigned pointcount_pd4 (unsigned long D[5], unsigned p);

// Point counting in F_p^2 for genus 2 curves
unsigned pointcount_g1_d2 (unsigned long f[4], unsigned p);
unsigned pointcount_g1d4_d2 (unsigned long f[5], unsigned p);
unsigned pointcount_g2_d2 (unsigned long f[6], unsigned p);
unsigned pointcount_g2d6_d2 (unsigned long f[7], unsigned p);

// This use half as much memory but are slower - use when p/8 > L2 cache
unsigned pointcount_big_g1 (unsigned long D[4], unsigned p);
unsigned pointcount_big_g2 (unsigned long D[6], unsigned p);
unsigned pointcount_big_g2d6 (unsigned long D[7], unsigned p, unsigned long f6);
unsigned pointcount_big_g3 (unsigned long D[8], unsigned p);
unsigned pointcount_big_g3d8 (unsigned long D[9], unsigned p, unsigned long f8);
unsigned pointcount_big_g4 (unsigned long D[10], unsigned p);

// These routines return point counts on 32 curves f(x), f(x)+1, ..., f(x)+32
int pointcount_multi_g2 (unsigned pts[], unsigned long D[6], unsigned p);
int pointcount_multi_g2d6 (unsigned pts[], unsigned long D[6], unsigned p, unsigned long f6);
int pointcount_multi_g3 (unsigned pts[], unsigned long D[8], unsigned p);
int pointcount_multi_g3d8 (unsigned pts[], unsigned long D[6], unsigned p, unsigned long f8);

// naive polynomial evaluation, provided for testing and to handle small cases
unsigned pointcount_slow (ff_t f[], int d, unsigned p);
unsigned pointcount_tiny (unsigned long f[], int d, unsigned p);
unsigned pointcount_tiny_pd4 (unsigned long f[], int d, unsigned p);

// naive polynomial evaluation for pointcounting over F_p^2 - only used for p <= POINTCOUNT_MAX_TINYP
unsigned pointcount_slow_d2 (ff_t f[], int d, unsigned p);
unsigned pointcount_tiny_d2 (unsigned long f[], int d, unsigned p);

// naive polynomial evaluation for pointcounting over F_p^3 - only used for p <= POINTCOUNT_MAX_TINYP
unsigned pointcount_tiny_d3 (unsigned long f[], int d, unsigned p);

#ifdef __cplusplus
}
#endif

#endif
