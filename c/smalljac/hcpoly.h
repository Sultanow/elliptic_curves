#ifndef _HCPOLY_INCLUDE_
#define _HCPOLY_INCLUDE_

/*
     Copyright (c) 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <assert.h>
#include "ff_poly.h"
#include "ecurve_ff2.h"
#include "hecurve.h"

typedef struct hc_poly_struct {
	ff_t f[2*(HECURVE_MAX_DEGREE+1)];	// poly of degree d=2g+1 or d=2g+2 over F_p^n, where n is 1 or 2.
	int d;							// degree of f
	int g;							// (d-1)/2
	int n;							// degree of number field, currently can be at most 2
	int two_rank;						// cached 2-rank of the Jacobian (negative value indicates it has not been computed yet)  only supported for prime fields
} hc_poly;

static inline void hc_poly_init (hc_poly *hc, ff_t f[], int d, int n)
	{ assert (d <= HECURVE_MAX_DEGREE && n <= 2); for ( int i = 0 ; i <= d*n ; hc->f[i] = f[i] );  hc->d = d;  hc->g = (d-1)/2;  hc->n = n; hc->two_rank = -1;}
		
// the functions below are currently only supported over prime fields
int hc_poly_compute_two_rank (hc_poly *hc);
static inline int hc_poly_two_rank (hc_poly *hc)
	{ if ( hc->two_rank >= 0 ) hc_poly_compute_two_rank (hc);  return hc->two_rank; }
		
static inline void hc_poly_copy (hc_poly *hc1, hc_poly *hc2) { *hc1 = *hc2; }
static inline void hc_poly_print (hc_poly *hc)
	{ if ( hc->n == 1 ) ff_poly_print(hc->f, hc->d); else if ( hc->n == 2 ) ff2_poly_print (hc->f, hc->d); else printf ("Can't currently print poly over number fields of degree %d\n", hc->n); }

// put hyperelliptic curve into a suitable form for efficient computation (see hc_poly.c for details)
void hc_poly_standardize (hc_poly *c);
		
// initialize an hc_poly with integer (mpz_t) coefficients (necessarily defined over Q)
static inline void hc_poly_set_mpz (hc_poly *hc, mpz_t F[], int d)
{
	hc->d = d;  hc->g = (d-1)/2; hc->n = 1;
	if ( (d&1) ) _ff_set_zero (hc->f[d+1]);								// this is necessary because the hecurve code doesn't look at hc->d, it will look at hc->f[2g+2] (which is hc->f[d+1] when d is odd)
	ff_poly_set_mpz (hc->f, F, d);										// reduce integer coefficients mod p
	if ( (d&1) && _ff_one(hc->f[d]) && _ff_zero(hc->f[d-1]) ) return;			// if f is monic, depressed, and of odd degree, we are good to go
	hc_poly_standardize (hc);											// otherwise standarize the poly
}

static inline void hc_poly_twist (hc_poly *t, hc_poly *c)
{
	assert (c->n == 1);
	hc_poly_copy (t, c);
	ff_poly_twist (t->f, c->f, c->d);									// note that ff_poly_twist does not change the factorization pattern but may change roots if the degree is odd
	if ( (t->d&1) && _ff_one(t->f[t->d]) && _ff_zero(t->f[t->d-1]) ) return;		// if f is monic, depressed, and of odd degree, we are good to go
	hc_poly_standardize (t);											// otherwise standarize the poly (note that this may be necessary even if c was already standardized)
}

static inline int hc_poly_discriminant_nonzero (hc_poly *hc)
{
	if ( hc->n == 1 ) {
		if ( hc->d == 3 && _ff_one(hc->f[3]) && _ff_zero(hc->f[2]) )
			{ ff_t D;  return ff_poly_x3axb_disc (&D, hc->f); }
		else return ff_poly_discriminant_nonzero (hc->f, hc->d);
	} else if ( hc->n == 2 ) {
		if ( hc->d == 3 && ff2_one(hc->f+6) ) {
			ff_t D[2];
			if ( _ff_p == 3 ) return ff2_poly_disc_char3 (D, hc->f);
			if ( ff2_zero(hc->f+4) ) return ff2_poly_x3axb_disc (D, hc->f);
		}
		return ff2_poly_discriminant_nonzero (hc->f, hc->d);
	}
	err_printf ("Unhandled case number field degree %d in hc_poly_discriminant_nonzero\n", hc->n); exit (0);
}

void hc_integral_model (mpz_t den, mpz_t f[], mpq_t F[], int d);

#ifdef __cplusplus
}
#endif

#endif
