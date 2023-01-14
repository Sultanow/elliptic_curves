#ifndef _MPZPOLYUTIL_INCLUDE_
#define _MPZPOLYUTIL_INCLUDE_

/*
    Copyright 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
	This module contains code for converting polynomials to/from strings, and some very basic code for computing
	discriminants, and resultants, (optimized for a special few cases but slow in general, only use for low degree polys).
*/

#define MPZ_POLY_DISC_MAX_DEGREE		10					// max degree supported by discriminant and resultant functions

struct mpz_rational_function_struct {
	int d_num, d_den;
	mpz_t *num, *den;
	mpq_t y;
};
typedef struct mpz_rational_function_struct mpz_rational_function_t[1];

/*
	Syntax for polynomial and rational function expressions:

		ratexpr = parpoly [ "/" parpoly ]		(NB: numerator and denominator *must* be parenthesized, even if they consist of a single term)
		polyexpr = parpoly | poly
		parpoly = "[" poly "]" | "(" poly ")"
		poly = term | term sign poly
		term = scoeff | [scoeff] var exp
		scoeff = [sign] coeff
		coeff = "" | zcoeff | qcoeff
		zcoeff = digit+
		qcoeff = zcoeff "/" zcoeff
		exp = "" | "^" digit+
		sign = "+" | "-"
		digit = 0-9
		var = a-zA-Z 					(NB: every occurence of var must be the same character)

	Note that parentheses ("[" or "(") are only permitted at the topmost level, one cannot specify polynomial arithmetic in a poly expression.
	Also note that terms may either be seperated by a single sign and have unsigned coefficients, or be separated by a sign and have signed coefficients (so + + and + - are okay).
	Parsing stops at the first character that breaks one of the rules above.
*/

int mpz_rational_function_init (mpz_rational_function_t f, char *expr);
void mpz_rational_function_clear (mpz_rational_function_t f);

// functions for working with Weierstrass coefficients  [a1,a2,a3,a4,a6] specifying an elliptic curve
static inline int mpz_ws_from_str(mpz_t W[5], char *cstr) { return ( gmp_sscanf (cstr, "[%Zd,%Zd,%Zd,%Zd,%Zd]", W[0], W[1], W[2], W[3], W[4]) == 5 ? 1 : 0 ); }
void mpq_ws_make_integral (mpq_t W[5]);
void mpz_short_ws (mpz_t A, mpz_t B, mpz_t a1, mpz_t a2, mpz_t a3, mpz_t a4, mpz_t a6);
static inline void mpz_short_ws_jinv_num_den (mpz_t num, mpz_t den, mpz_t A, mpz_t B)
{
	mpz_mul(num,A,A); mpz_mul(num,num,A); mpz_mul_2exp(num,num,2);			// num = 4A^3
	mpz_mul(den,B,B); mpz_mul_ui(den,den,27); mpz_add(den,den,num);			// den = 4A^3+27B^2
	mpz_mul_ui(num,num,1728);												// num = 1728 * 4A^3
}
static inline int mpz_short_ws_jinv (mpq_t J, mpz_t A, mpz_t B)
{
	mpz_t num, den;
	
	mpz_init(num); mpz_init(den);
	mpz_short_ws_jinv_num_den (num, den, A, B);
	if ( ! mpz_sgn(den) ) { mpz_clear (num); mpz_clear(den); return 0; }
	mpq_set_num(J,num);  mpq_set_den(J,den);  mpq_canonicalize (J);
	mpz_clear (num); mpz_clear(den);
	return 1;
}
static inline void mpz_short_ws_from_jinv (mpz_t A, mpz_t B, mpq_t J)
{
	mpz_t X;
	
	if ( !mpq_sgn(J) ) { mpz_set_ui(A,0); mpz_set_ui(B,1);  return; } 
	if ( mpz_cmp_ui(mpq_denref(J),1)==0 && mpz_cmp_ui(mpq_numref(J),1728)==0 ) { mpz_set_ui(A,1); mpz_set_ui(B,0); return; }
	mpz_init (X);
	mpz_mul_ui (X, mpq_denref(J), 1728);  mpz_sub (X, X, mpq_numref(J));			// X = den(J)(1728-J)
	mpz_mul_ui (A, mpq_numref(J), 3); mpz_mul (A, A, X);							// A = 3J(1728-J) * den(J)^2
	mpz_add (B, mpq_numref(J), mpq_numref(J)); mpz_mul (B,B,X);  mpz_mul (B,B,X);	// B = 2J(1728-J)^2 * den(J)^3
	mpz_clear (X);
	return;
}

int mpz_jinv_has_cm (mpz_t J);
static inline int mpq_jinv_has_cm (mpq_t J)
{
	if ( mpz_cmp_ui (mpq_denref(J), 1) != 0 ) return 0;
	return mpz_jinv_has_cm (mpq_numref(J));
}	

static inline int mpz_short_ws_has_cm (mpz_t A, mpz_t B)
{
	mpz_t num, den;
	int cm;
	
    if ( !mpz_sgn(A) || !mpz_sgn(B) ) return 1;
	mpz_init(num); mpz_init(den);
	mpz_short_ws_jinv_num_den (num, den, A, B);
	if ( ! mpz_divisible_p (num, den) ) return 0;
	mpz_divexact (num, num, den);
	cm = mpz_jinv_has_cm (num);
	mpz_clear (den);  mpz_clear (num);
	return cm;
}

static inline int mpz_ws_has_cm (mpz_t a1, mpz_t a2, mpz_t a3, mpz_t a4, mpz_t a6)
{
	mpz_t A, B;
	int cm;
	
	mpz_init(A); mpz_init(B);
	mpz_short_ws (A, B, a1, a2, a3, a4, a6);
	cm = mpz_short_ws_has_cm (A, B);
	mpz_clear (A); mpz_clear (B);
	return cm;
}

// poly parse functions return degree (-1 for 0 poly), -2 if deg exceeds maxd, and -3 in the case of an error
int mpq_poly_parse (mpq_t f[], int maxd, char *expr);
int mpz_poly_parse (mpz_t f[], int maxd, char *expr);

// for plane quartics f[] must have space for 15 coefficients
int mpq_poly_parse_plane_quartic (mpq_t f[], char *expr);
int mpz_poly_parse_plane_quartic (mpz_t f[], char *expr);

static inline void mpz_poly_eval (mpz_t y, mpz_t f[], int d, mpz_t x)				// y cannot coincide with x
{
	register int i;
	
	if ( d < 0 ) { mpz_set_ui(y,0); return; }
	mpz_set (y, f[d]);
	for ( i = d-1 ; i >= 0 ; i-- ) { mpz_mul (y, y, x);  mpz_add (y, y, f[i]); }
}

static inline void mpz_poly_eval_q (mpq_t y, mpz_t f[], int d, mpq_t x)				// y cannot coincide with x
{
	mpq_t z;
	register int i;
	
	if ( d < 0 ) { mpq_set_ui(y,0,1); return; }
	mpq_init (z);
	mpq_set_z (y, f[d]);
	for ( i = d-1 ; i >= 0 ; i-- ) { mpq_mul (y, y, x); mpq_set_z (z, f[i]);  mpq_add (y, y, z); }
}

static inline int mpz_rational_function_eval (mpq_t y, mpz_rational_function_t f, mpq_t x)
{
	mpz_poly_eval_q (y, f->num, f->d_num, x);
	if ( ! f->den ) return 1;
	mpz_poly_eval_q (f->y, f->den, f->d_den, x);
	if ( ! mpq_sgn (f->y) ) return 0;
	mpq_div (y, y, f->y);
	return 1;
}

void mpz_poly_print (mpz_t f[], int df);
int mpz_poly_sprint (char *buf, mpz_t f[], int df);

static inline void mpz_rational_function_print (mpz_rational_function_t f)
{
	mpz_poly_print (f->num, f->d_num);
	if ( f->den ) { printf (" / "); mpz_poly_print (f->den, f->d_den); }
}

unsigned long i_poly_set_mpq (long f[], mpq_t F[], int d);									// sets f[i] to numberator of F[i] with common denominator and returns denominator or zero if overflow
static inline void ui_poly_set_mpz_mod_p (unsigned long f[], mpz_t F[], int d, unsigned long p)
	{ for ( int i = 0 ; i <= d ; i++ ) f[i] = mpz_fdiv_ui (F[i], p); }
static inline void i_poly_set_mpz (long f[], mpz_t F[], int d)
	{ for ( int i = 0 ; i <= d ; i++ ) f[i] = mpz_get_si (F[i]); }

int mpz_poly_depress_monic_inplace (mpz_t f[], int d);
void mpz_poly_discriminant (mpz_t D, mpz_t f[], int d, mpz_t *w);								// thread safe generic discriminant code -- w must point to d*d+2 initialized mpz_t's
void mpz_fast_poly_discriminant (mpz_t disc, mpz_t f[], int d);								// optimized NON-THREAD-SAFE code, checks for lots of special cases, avoids reallocing, etc...
void mpz_poly_resultant (mpz_t R, mpz_t f[], int d_f, mpz_t g[], int d_g);

int mpz_hyperelliptic_curve_parse (mpz_t f[], int *pdf, mpz_t h[], int *pdh, int maxd, char *str);		// maxd bounds the degree of f, (maxd+1)/2 bounds the degree of h, returns genus (0 for failure)
void mpz_hyperelliptic_curve_discriminant (mpz_t D, mpz_t f[], int df, mpz_t h[], int dh);			// uses thread safe discriminant code
int mpz_hyperelliptic_curve_normalize (mpz_t g[], mpz_t f[], int df, mpz_t h[], int dh);				// sets g to 4*f+h^2, or just to f if h is zero (g is allowed to alias f but not h)
static inline int mpz_hyperelliptic_curve_sprint (char *s, mpz_t f[], int df, mpz_t h[], int dh)
{
	char *t = s + mpz_poly_sprint (s, f, df) ;
	char *u = t + mpz_poly_sprint (t, h, dh);
	*(t-1) = ','; *t = ' ';
	return u-s;
}

// g=f is ok
static inline void mpz_poly_derivative (mpz_t g[], mpz_t f[], int d)
	{ register int i;  if ( d<= 0 ) return;  mpz_set (g[0],f[1]);  for ( i = 2 ; i <= d ; i++ ) mpz_mul_ui(g[i-1],f[i],i); }

// g=f is ok
static inline void mpz_poly_reduce (mpz_t g[], mpz_t f[], int d, mpz_t m)
	{ register int i;  if ( d < 0 ) return;  for ( i = 0 ; i <= d ; i++ ) mpz_mod (g[i], f[i], m); }

// y = f(x) mod p
static inline void mpz_poly_eval_mod (mpz_t y, mpz_t f[], int d, mpz_t x, mpz_t m)				// y cannot coincide with x
{
	register int i;
	
	if ( d < 0 ) { mpz_set_ui(y,0); return; }
	mpz_mod (y, f[d], m);
	for ( i = d-1 ; i >= 0 ; i-- ) { mpz_mul (y, y, x);  mpz_add (y, y, f[i]); mpz_mod (y,y,m); }
}

#ifdef __cplusplus
}
#endif

#endif
