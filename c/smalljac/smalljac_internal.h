#ifndef _SMALLJAC_INTERNAL_
#define _SMALLJAC_INTERNAL_

/*
    Copyright (c) 2007-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>
#include "ff_poly.h"
#include "smalljac.h"
#include "jac.h"
#include "nfpoly.h"
#include "igusa.h"

#ifdef __cplusplus
extern "C" {
#endif

// Internal smalljac declarations that are not part of the public interface

#define SMALLJAC_MAX_CONSTRAINTS		64
#define SMALLJAC_MAX_COEFFICIENTS		15					// must be at least max(15,SMALLJAC_MAX_DEGREE+1) (the 15 is for plane quartics)
#define SMALLJAC_MAX_NF_DEGREE		12					// maximum degree of number fields supported by smalljac (increase as needed, it just takes up more memory)
#define SMALLJAC_MAX_BAD_PRIMES		32
#define SMALLJAC_MAX_FACTOR_BITS		64					// don't waste time trying to factor anything hard
#define SMALLJAC_SMALL_INTERVAL		1000				// don't bother factoring for small intervals

#define SMALLJAC_SPECIAL_X6PA			1					// y^2=x^6+a
#define SMALLJAC_SPECIAL_X5PAX			2					// y^2=x^5+ax
#define SMALLJAC_SPECIAL_X5PA			3					// y^2=x^5+a (not yet fully supported, falls back to generic code for p=1 mod 5)
#define SMALLJAC_SPECIAL_PICARD			4					// y^3 = quartic(x) (Picard curve, genus 3)
#define SMALLJAC_SPECIAL_CM			5					// set for elliptic curves defined over Q with CM
#define SMALLJAC_SPECIAL_X7PAX			6					// y^2=x^7+ax
#define SMALLJAC_SPECIAL_X8PA			7					// y^2=x^8+a
#define SMALLJAC_SPECIAL_X6P1_TWIST	8					// twist of y^2=x^6+1
#define SMALLJAC_SPECIAL_X5PX_TWIST	9					// twist of y^2=x^5+x
#define SMALLJAC_SPECIAL_CM2SQUARE	10					// y^2=x^6-5x^4-5x^2+1 (Jacobian Q-isogenous to the square of an elliptic curve with CM discriminant -2)
#define SMALLJAC_SPECIAL_FK_TWIST		11					// special plane quartic (twist of Fermat or klein curve)

#define SMALLJAC_CURVE_FLAG_DELTA		0x1					// indicates Delta values have been precomputed
#define SMALLJAC_CURVE_FLAG_WS		0x2					// indicates genus 1 curve specified in Weierstrass form

#define SMALLJAC_CURVE_ELLIPTIC			1					// elliptic curve
#define SMALLJAC_CURVE_HYPERELLIPTIC	2					// hyperelliptic curve y^2+h(x)y=f(x) (curves of genus 1 also allowed)
#define SMALLJAC_CURVE_PICARD			3					// y^3=quartic, distinguished from other plane quartics
#define SMALLJAC_CURVE_PLANE_QUARTIC	4					// plane quartic

#define SMALLJAC_GROUP_FLAG			0x1000

typedef struct smalljac_curve_struct {
	mpz_t f[SMALLJAC_MAX_COEFFICIENTS];						// Integer poly f(x) for which the curve is represented as y^2=f(x) for all p > 2g+1 (y^3=f(x) for Picard curves, coeffs of quartic poly in x and y for plane quartics)
	mpz_t Deltas[SMALLJAC_MAX_DEGREE+1];					// Deltas[k] = (\Delta^k f)(0), the k-th difference poly of f evaluated at 0
	mpz_t D;												// discriminant of the curve -- primes dividing D are of bad reduction,  for p > 2g+1 these will be the same as the primes dividing disc=disc(f)
	mpz_t disc;												// discriminant of the poly f(x) (this is not necessarily the same as the discriminant of the curve)
	mpz_t F[2*SMALLJAC_MAX_GENUS+3];						// for elliptic and hyperelliptic curves over Q, the curve is y^2 + H(x)*y = F(x)
	mpz_t H[SMALLJAC_MAX_GENUS+2];
	int dF, dH;
	int f_inits, Delta_inits, F_inits, H_inits;						// # of elements initiailized
	int degree;												// degree of f(x)
	int genus;												// genus of the curve
	unsigned flags;											// mask of boolean flags indicating the state of initialization of other parameters
	unsigned ws2;											// for Weierstrass specified curves, holds the coefficients mod 2 in bottom 5 bits
	unsigned ws3;											// for Weierstrass specified curves, holds the coefficients mod 3 in bottom 10 bits (2 bits per)
	long a[2*SMALLJAC_MAX_GENUS];							// L_p(T) coefficients (or group structure coefficients) most recently computed
	int n;													// number of coefficients
	int type;												// set to one of the SMALLJAC_CURVE_xxxxx values defined above (elliptic/hyperelliptic/picard/plane-quartic)
	int special;											        // flag for special curves, e.g. x^6+a
        int special_curve_id;                                                                            // identifier used to distinguish various special curves of the same type
	char str[SMALLJAC_CURVE_STRING_LEN];					// pointer to original curve specification - null terminated
	char *nfstr;												// points to number field specificiation in str, if present, null ow
	nf_poly *nfp;												// pointer to number field poly structure, used for curves defined over number fields
	int nfd;													// degree of number field (1 for Q -- in which case nfp is null)
	int Qflag;												// set whenever the curve is defined over Q (even if it is being considered over a number field)
	hc_poly hc[1];											// local reduction of curve at a prime over p (over Q this is just the curve reduced mod p)
	long q;													// current prime (or prime power) being processed
	long pts;												// pointcount over F_q, if performed, zero o.w.
} smalljac_curve;

int padic_charpoly(long a[], long f[], int n, unsigned long p);	// c -> c++ interface function for David Harvey's frobenius() code - used in genus 3 only

int smalljac_x3pa_Lpoly (long a[], ff_t f0);
int smalljac_x3pax_Lpoly (long a[], ff_t f1);
int smalljac_x5pax_Lpoly (long a[], ff_t f1);
int smalljac_x6pa_Lpoly (long a[], ff_t f0);
int smalljac_cm2square_Lpoly (long a[]);
int smalljac_x7pax_a1 (long a[], ff_t f1);
int smalljac_x8pa_a1 (long a[], ff_t f1);
int smalljac_FKtwist_lookup (mpz_t C[15]);
int smalljac_FKtwist_Lpoly (long a[], int id);

int smalljac_internal_Lpoly_Q (long a[], smalljac_curve *sc, long p, unsigned long flags);
int smalljac_internal_Lpoly_nf (long a[], hc_poly *c, long p, unsigned long flags);
int smalljac_tiny_Lpoly (long a[], smalljac_curve *sc, int p, unsigned long flags);
int smalljac_generic_Lpoly (long a[], hc_poly *hc, long pts, unsigned long flags);
int smalljac_padic_Lpoly (long a[], smalljac_curve *sc, long p, unsigned long flags);
unsigned long smalljac_pointcount_modp (smalljac_curve *sc,  long p);

int smalljac_Lpoly_extend (long a[], int n, long p, int h);						        // extend coefficients for prime field p to extension field of size q = p^h			

static inline unsigned long smalljac_curve_max_p (smalljac_curve *sc)
{
	if ( sc->special && (sc->flags&SMALLJAC_A1_ONLY) ) return (1L<<44);		// MAX_ENUM_PRIME
	if ( sc->special && !(sc->flags&SMALLJAC_GROUP) ) return (1L<<44);			// MAX_ENUM_PRIME
	return smalljac_max_p (sc->genus);
}

// the following functions update an existing curve structure to reflect new curve parameters without reallocating.
int smalljac_curve_set_str (smalljac_curve *sc, char *str, int *err);
int smalljac_curve_set_mpz (smalljac_curve *sc, mpz_t f[], int degree,  char *str);		// note str is not validated
int smalljac_curve_set_i (smalljac_curve *sc, long f[], int degree, char *str);			// ditto

smalljac_curve *smalljac_curve_alloc ();										        // simply allocates an unitialized curve structure to be used above

// This function does not check for bad reduction (TODO: fix this!)
static inline int smalljac_Qcurve_reduce (ff_t f[], smalljac_curve *sc)
	{ if ( ! (sc->Qflag) ) return 0; ff_poly_set_mpz (f, sc->f, sc->degree); return 1; }

static inline int smalljac_curve_degree (smalljac_curve *sc)
	{ return sc->degree; }
	
static inline int smalljac_curve_lc_divisible_p (smalljac_curve *sc,  unsigned long p)	        // returns 1 if leading coefficient is divisible by p
	{ return mpz_divisible_ui_p (sc->f[sc->degree], p); }

static inline int smalljac_curve_cc_divisible_p (smalljac_curve *sc,  unsigned long p)	// returns 1 if  constant coefficient is divisible by p
	{ return mpz_divisible_ui_p (sc->f[0], p); }

static inline void smalljac_curve_igusa_invariants (smalljac_curve *sc, mpq_t I[3])
	{ mpq_poly_igusa_inv (I, sc->f, sc->degree); }
	
static inline const char *smalljac_curve_string (smalljac_curve *sc)
	{ return sc->str; }
	
void smalljac_init (void);		// will automatically be called when needed
	
int smalljac_analyze_poly_string (char *str, int *ws, int *nf);

// computes L_p(1) - does not check for overflow
static inline unsigned long smalljac_Lp1_ui (long a[], int genus, unsigned long p)
{
	switch (genus) {
	case 1: return (unsigned long)((long)(p+1)+a[0]);
	case 2: return (unsigned long)((long)(p*p+1)+(long)(p+1)*a[0]+a[1]);
	case 3: return (unsigned long)((long)(p*p*p+1)+(long)(p*p+1)*a[0]+(long)(p+1)*a[1]+a[2]);
	default: printf ("unhandled genus %d\n", genus);  exit (0);
	}
}

#ifdef __cplusplus
}
#endif

#endif
