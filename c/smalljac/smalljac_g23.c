#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gmp.h>
#include "mpzutil.h"
#include "ff_poly.h"
#include "hecurve.h"
#include "jac.h"
#include "jacorder.h"
#include "smalljac.h"
#include "smalljac_internal.h"
#include "smalljac_g23.h"
#include "pointcount.h"
#include "cstd.h"

/*
    Copyright (c) 2007-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

static inline void mpz_add_si (mpz_t rop, mpz_t op1, long op2) { if ( op2 < 0 ) mpz_sub_ui (rop, op1, -op2); else mpz_add_ui (rop, op1, op2); }
static inline void mpz_sub_si (mpz_t rop, mpz_t op1, long op2) { if ( op2 < 0 ) mpz_add_ui (rop, op1, -op2); else mpz_sub_ui (rop, op1, op2); }

int smalljac_distinguish_group_orders (mpz_t N[], int k, hc_poly c[1]);

/*
    This module contains genus 2 and genus 3 specific code for deriving the
    numerator of the zeta function, P(T) = L_p(T), when given certain
    partial information such as the coefficients mod p or the values P(1) and P(-1).
*/

/*
	Derive P(T) from P(1) and P(-1) in genus 2.  No group operations required.  Assumes p < 2^31
	Sanity check coefficients - caller may not be sure about P(1) and P(-1) and may use this to narrow down the possibilities
*/
int smalljac_genus2_charpoly (long a[3], long p, double sqrtp, long P1, long PN1)
{
	if ( p > (1L<<31) ) { err_printf ("Overflow in smalljac_genus2_charpoly, p = %ld\n", p);  return 0; }
	a[0] = (P1-PN1)/(2*(p+1));
	a[1] = ((P1+PN1) - 2*(p*p+1)) / 2;
	if ( (p*p+1) + (p+1)*a[0] + a[1] !=  P1 ) return 0;
	if ( (p*p+1) - (p+1)*a[0] + a[1] != PN1 ) return 0;
	if ( abs(a[0]) > 4*sqrtp ) return 0;										// standard Weil bound
	if ( a[1] > 2*p + a[0]*a[0]/4 || a[1] < -2*p+a[0]*a[0]/2 ) return 0;			// these bounds follow from Prop. 4 of [KedlayaSutherland2007]
	if ( abs(a[2]) > 20*p*sqrtp ) return 0;
	return 1;
}


/*
	Derive P(T) from P(1), P(-1) and #C/F_p in genus 3.  No group operations required.  Assumes p<2^20
	Sanity check coefficients - caller may not be sure about P(1) and P(-1) and may use this to narrow down the possibilities
*/
int smalljac_genus3_charpoly (long a[3], long p, double sqrtp, long P1, long PN1, unsigned long pts)
{
	if ( p > (1L<<20) ) { err_printf ("Overflow in smalljac_genus3_charpoly, p = %ld\n", p);  return 0; }
	a[0] = (long)pts - (long)(p+1);
	a[2] = ((P1-PN1) - 2*(p*p+1)*a[0]) / 2;
	a[1] = (P1 - (p*p*p+1) - (p*p+1)*a[0] - a[2])/(p+1);
	if ( (p*p*p+1) + (p*p+1)*a[0] + (p+1)*a[1] + a[2] !=  P1 ) return 0;
	if ( (p*p*p+1) - (p*p+1)*a[0] + (p+1)*a[1] - a[2] != PN1 ) return 0;
	if ( abs(a[0]) > 6*sqrtp ) return 0;												// standard Weil bound
	if ( a[1] > 3*p + a[0]*a[0]/3 || a[1] < -3*p +a[0]*a[0]/2 ) return 0;					// these bounds on a2 follow from Prop. 4 of [KedlayaSutherland2007]
	if ( abs(a[2]) > 20*p*sqrtp ) return 0;
	return 1;
}


/*
	Derive P(T) from P(1) in genus 2, using O(\lg(p)) gops in the twist.
*/

int smalljac_genus2_charpoly_from_P1 (long a[2], long P1, long Min, long Max, hc_poly c[1])
{
	hc_poly twist;
	jac_t h, g;
	unsigned long n;
	long e, o, to, to0, a0, p;
	int i,k;
	
	n = jac_gops;
	hc_poly_twist (&twist, c);

	// We want to compute a and b such that P(T) = 1 + a1T + a2T^2 + pa1T^3 + p^2T^4
	// We are given P(1) and we know that P(-1) is the order of the twist
	// We also know that a1 >= (P(1) - (p^2+1) - 6p)/(p+1) and a <= (P(1) - (p^2+1) + 6p)/(p+1)
	// which gives a range of at most 12 values for a1, and thus at most 12 possibilities for P(-1) = P(1) - 2(p+1)a1
	// all of which differ by a multiple of 2(p+1).  We just need to distinguish P(-1) among these 11 possibilities.
	p = _ff_p;
	a0 = (P1 - p*p + 6*p - 1)/(p+1);						// max value for a1 = (P(1) - (p^2+6p+1)) / (p+1) gives min value for twist order to
	to = P1 - 2*(p+1)*a0;
	to0 = to;
	i = 0;
	while ( to < Min ) { i++; to += 2*(p+1); }				// Skip twist orders not in Weil interval
	for ( k = 0 ; k < SMALLJAC_RETRIES ; k++ ) {
		_jac_random(g, twist);
		jac_exp_ui (&h, &g, 2*(p+1), &twist);
		if ( ! _jac_is_identity(h) ) break;
	}
	if ( k == SMALLJAC_RETRIES ) {
		// 2*(p+1) is an exponent for the group - hopefully only one candidate twist order is compatible, since they are all exponents
		o = 0;
		for ( ; i <= 12 ; i++, to += 2*(p+1) ) {
			if ( to > Max ) break;
			if ( ui_compatible (to, 2*(p+1)) ) {
				if ( o ) { /*err_printf ("%7u: More than one twist order implied by P(1) = %ld is a 2(p+1) compatible exponent of the twist (%ld,%ld).\n", _ff_p, P1, o, to);*/  return 0; }
				o = to;
			}
		}
	} else {
		// Ok we have an element whose order contains a prime not in 2*(p+1).  This should almost always uniquely determine the twist
		// since two candidate twist orders differ by i*2*(p+1) where i <= 10
		o = 0;
		jac_exp_ui (&g, &g, to, &twist);
		for ( ; i <= 12 ; i++, to += 2*(p+1) ) {
			if ( to > Max ) break;
			if ( _jac_is_identity (g) ) {
				if ( o ) {	// we know that e=(to-o) is an exponent of g - try to find an element for which e is not an exponent
					e = to-o;
					for ( k = 0 ; k < SMALLJAC_RETRIES ; k++ ) {
						_jac_random(g, twist);
						jac_exp_ui (&h, &g, e, &twist);
						if ( ! _jac_is_identity(h) ) break;
					}
					if ( k == SMALLJAC_RETRIES ) {/*err_printf ("%7u: More than one twist order implied by P(1) = %ld is an exponent of the twist (%ld,%ld).\n", _ff_p, P1, o, to);*/  return 0; }
					jac_exp_ui (&h, &g, o, &twist);
					if ( _jac_is_identity(h) ) {
						// to can't be an exponent, and o is probably the twist order, reset g and h to continue the search
						jac_exp_ui (&h, &g, 2*(p+1), &twist);
						jac_exp_ui (&g, &g, to, &twist);
					} else {
						o = 0;
						jac_exp_ui (&h, &g, 2*(p+1), &twist);
						jac_exp_ui (&g, &g, to, &twist);
						if ( _jac_is_identity (g) ) o = to;
					}
				} else {
					o = to;
				}
			}
			_jac_mult (g, g, h, twist);
		}
	}
	if ( ! o ) { err_printf ("%7lu: None of the implied twist orders for P(1) = %ld are exponents with to0 = %ld in Weil interval [%ld,%ld]\n", _ff_p, P1, to0, Min, Max);  return 0; }
	a[0] = (P1 - o) / (2*(p+1));		// a1 = P(1)-P(-1)/(2(p+1))
	a[1] = (P1+o-2*(p*p+1))/2;		// a2 = (P(1)+P(-1) - 2(p^2+1))/2
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}


/*
	Derive P(T) from P(1) in genus 3, using O(p^{1/4}) gops in the twist.
*/

int smalljac_genus3_charpoly_from_P1 (long o[3], long P1, long pts, long Min, long Max, hc_poly c[1])
{
	static int init;
	static mpz_t Tmin, Tmax, e[2];
	hc_poly twist;
	jac_t g;
	long p, tmin, tmax, to, spc, bmin, bmax, k;
	unsigned long n, e0;
	double x;
	int constraints[2];
	int cnt;
	
	n = jac_gops;
	if ( ! init ) { mpz_init (Tmin);  mpz_init (Tmax);  mpz_init (e[0]);  mpz_init (e[1]);  init = 1; }

	hc_poly_twist (&twist, c);

	p = (long) _ff_p;
	o[0] = pts - 1 - p;
	x = sqrt((double)p);
	spc = (long)ceil(x);
	bmin = (P1 - (p*p*p+1) - (p*p+1)*o[0] - 20*p*spc) / (p+1);
	tmin = 2*(p*p*p+1) - P1+ 2*(p+1)*bmin;
	bmax = (P1 - (p*p*p+1) - (p*p+1)*o[0] + 20*p*spc) / (p+1) + 1;	
	tmax = 2*(p*p*p+1) - P1+ 2*(p+1)*bmax;
	if ( tmin < Min ) tmin += _ui_ceil_ratio((Min-tmin),2*(p+1)) * 2*(p+1);
	tmax = 2*(p*p*p+1) - P1+ 2*(p+1)*15*p;
	if ( tmax > Max ) tmax -= (tmax-Max)/(2*(p+1)) * 2*(p+1);
	_jac_random (g,twist);
	mpz_set_ui (Tmin, tmin);  mpz_set_ui (Tmax, tmax);
	cnt = jac_search (e, &g, 2*(p+1), Tmin, Tmax, 0, &twist);
	if ( ! cnt ) { err_printf ("%7lu: Search failed in genus3_charpoly with P(1) = %ld in [%ld,%ld]\n", _ff_p, P1, Min, Max);  return 0; }
	if ( cnt > 1 ) {
		k =  2*(p*p*p+1) - P1;		// note k may be big
		constraints[0] = ( k < 0 ? (int) (2*(p+1) - ((-k)%(2*(p+1)))) : (int) (k%(2*(p+1))) );
		constraints[1] = -1;
		if ( jac_order ((unsigned long *)&to, Min, Max, -o[0], 0, 0, constraints, &twist) != 1 ) {
			gmp_err_printf ("%7lu: Ambiguous result in genus3_charpoly with P(1) = %ld in [%ld,%ld] (%Zd,%Zd), order computation for twist failed\n", _ff_p, P1, Min, Max, e[0], e[1]);
			return 0;
		}
	} else {
		e0 = mpz_get_ui (e[0]);
		to = _ui_ceil_ratio(tmin,e0)*e0;
	}
	
	if ( ! smalljac_genus3_charpoly (o, p, x, P1, to, pts) ) { err_printf ("%7lu: Unable to derive consistent Lpolynomial from P(1) = %lu and P(-1) = %lu given a1 = %ld\n", p, P1, to, o[0]);  return 0; }
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}

/*
	Derive P(T) given P(T) mod p in genus 2, using O(lg(p)) group operations.
*/
int smalljac_genus2_charpoly_from_Pmodp (long a[2], hc_poly c[1])
{
	static mpz_t e0, e1, Min, Max;		// use mpz's so  we can handle Jacobians with order > 2^64
	static int init;
	jac_t g, h;
	unsigned long n;
	long m, p, w, a1, a2, a2min, a2max;
	double x;
	
	if ( ! init ) { mpz_init (Min);  mpz_init (Max);  mpz_init (e0);  mpz_init (e1); init = 1; }
	assert ( p < (1L<<40) );		// for convenience, we asssume 8*p^(3/2) < 2^63

	p = (long) _ff_p;
	a1 = a[0];  a2 = a[1];
	while ( a2 < 0 ) a2 += p;		// make a2 positive
	x = sqrt(p);
	while ( a1 > -4.0*x ) a1 -= p;		// set a1 to min value in [-4sqrt(p),4sqrt(p)] congruent to a[0] mod p
	a1 += p;
	if ( a1 > 4.0*x ) { err_printf ("%9lu: Invalid trace %ld modulo p, violates Weil bounds\n", p, a1); return 0; }
	if ( a1 + p < 4.0*x ) { err_printf ("%9lu: Trace %ld is ambiguous modulo p\n", p, a1); return 0; }
	a2min = a1 < 0 ? floor(-2.0*a1*x-2.0*p) : floor(2.0*a1*x-2.0*p);
	a2max = ceil ((a1*a1)/4.0 +2.0*p);
	
	n = jac_gops;
	mpz_set_ui (Min, p);  mpz_mul_ui (Min, Min, p);  mpz_add_ui (Min, Min, 1); mpz_set (Max, Min);
	mpz_set_si (e0, (p+1)*a1+a2min); mpz_add (Min, Min, e0);				// Set Min = p^2+1+(p+1)a1+a2min
	mpz_set_si (e0, (p+1)*a1+a2max); mpz_add (Max, Max, e0);				// Set Max = p^2+1+(p+1)a1+a2max
	m = mpz_fdiv_ui (Min, p);
	w = (1+p+a1+a2) % p;
	mpz_add_ui (Min, Min, m <= w ? w-m : p-(m-w));	// Make Min congruent to P(1) mod p
	mpz_set (e0, Min);
	_jac_random(g, c[0]);
	jac_exp_mpz (&h, &g, e0, c);
	jac_exp_ui (&g, &g, p, c);
	mpz_set_ui (e1, 0);
	while ( mpz_cmp (e0, Max) <= 0 ) {
		if ( _jac_is_identity (h) ) {
			if ( mpz_sgn(e1) ) { gmp_err_printf ("%9lu: Ambiguous result in genus2_charpoly_from_Pmodp %Zd and %Zd are both possible exponents\n", p, e1, e0);  return 0; }
			mpz_set (e1, e0);
		}
		_jac_mult (h, h, g, c[0]);
		mpz_add_ui (e0, e0, p);
	}
	if ( ! mpz_sgn(e1) ) { gmp_err_printf ("%9lu: Search failed in genus2_charpoly_from_Pmodp, a[1] = %ld, Min = %Zd, Max = %Zd\n", _ff_p, a[1], Min, Max);  return 0; }
	
	// set a[1] = a2 = e1 - (p*p+1+(p+1)*a1)
	mpz_set_ui (e0, p); mpz_mul_ui (e0,e0,p);  mpz_sub (e1,e1,e0); mpz_set_si (e0,1+(p+1)*a1);  mpz_sub(e1,e1,e0);
	a[1] = mpz_get_si (e1);
	smalljac_charpoly_gops += jac_gops-n;
	return 1;
}

// Bounds on a3 given a1,a2 are take from Theorem 1.1 of Haloui, "The characteristic polynomials of abelian varieties of dimension 3 over finite fields", JNT 130 (2010) 2745-2752
static inline long a3min (long a1, long a2, long p)
{
	register double x, y;
	
	if ( a2 == -p ) return -2*a1*p;	// if a2 == -p then we the L-poly has 2 real roots and we must have a3 = -2*a1*p
	x = sqrt(a1*a1-3*a2+9*p);
	x = -2.0*(a1*a1*a1)/27.0 + (a1*a2)/3.0 + (double)p*a1 - 2.0/27.0*x*x*x;
	y = sqrt(p);
	y = -2.0*p*a1 - 2*y*a2 - 2*p*y;
	return (long) floor ((x > y ? x : y));	// err on the low side in case of rounding errors
}

static inline long a3max (long a1, long a2, long p)
{
	register double x, y;
	
	if ( a2 == -p ) return -2*a1*p;	// if a2 == -p then we the L-poly has 2 real roots and we must have a3 = -2*a1*p
	x = sqrt(a1*a1-3*a2+9*p);
	x = -2.0*(a1*a1*a1)/27.0 + (a1*a2)/3.0 + (double)p*a1 + 2.0/27.0*x*x*x;
	y = sqrt(p);
	y = -2.0*p*a1 + 2*y*a2 + 2*p*y;
	return (long) ceil ((x < y ? x : y));		// err on the high side in case of rounding errors
}

/*
	Derive P(T) given P(T) mod p in genus 3, using O(p^{1/4}) group operations.

	We assume p is big enough to avoid any unpleasant special cases but less than 2^31
*/
int smalljac_genus3_charpoly_from_Pmodp (long a[3], hc_poly c[1])
{
	mpz_t e[2], N[2], w1, w2, w3, w4, w5, w6, Min, Max, *M;
	hc_poly ct[1];
	jac_t g, gt, h;
	long p, a1, a2, a2Min, a2Max, a3Min, a3Max, good_a2, spc, t;
	unsigned long n;
	double x;
	int i, j, k, r, cnt, sts;

	assert ( SMALLJAC_GENUS == 3 );
	mpz_init (e[0]); mpz_init (e[1]); mpz_init (w1); mpz_init(w2); mpz_init(w3); mpz_init(w4); mpz_init (w5); mpz_init (w6); mpz_init (Min); mpz_init (Max); mpz_init (N[0]); mpz_init (N[1]);

	sts = 0;
	
	p = (long) _ff_p;
	x = sqrt((double)p);
	spc = (long)ceil(x);

	cnt = 0;
	n = jac_gops;
	
	// compute a1 exactly, do a pointcount if necessary
	while ( a[0] < 0 ) a[0] += p;
	while ( a[0] >= p ) a[0] -= p;
	if ( a[0] > 6*x ) a[0] -= p;
	assert ( a[0] > -6*x && a[0] < 6*x );
	if ( a[0] - p > -6*x || a[0] + p < 6*x ) {
		unsigned long f[9];
		assert ( c->d <= 8 );
		for ( i = 0 ; i <= c->d ; i++ ) f[i] = _ff_get_ui (c->f[i]);
		a1 = pointcount_tiny (f, c->d, p) - p - 1;
		assert ( (a1-a[0])%p == 0 );
		a[0] = a1;
	}
	r = hc_poly_compute_two_rank (c);
	
	// put a2 and a3 in [0,p-1]
	while ( a[1] < 0 ) a[1] += p;
	while ( a[1] >= p ) a[1] -= p;
	while ( a[2] < 0 ) a[2] += p;
	while ( a[2] >= p ) a[2] -= p;
	// Set a2 bounds based on Proposition 4 of KedlayaSutherland 2007, we have -p <= a2 <= 3p + a1^2/3, and for |a1| > 2 we also have 4|a1| - 9 < a2
	t = 4 * ( a[0] < 0 ? - a[0] : a[0] ) - 9*spc;
	if ( t > 0 ) {
		for ( a2Min = a[1] ; a2Min < t ; a2Min += p );
	} else {
		a2Min = -p+a[1];
	}
	a2Max = 3*p + _ui_ceil_ratio(a[0]*a[0], 3);
	mpz_set_ui (w1, p);  mpz_mul_ui (w3, w1, p);  mpz_set (w2, w3); mpz_mul_ui (w3, w3, p);  mpz_add_ui (w3, w3, 1);  mpz_add_ui (w2, w2, 1);		// w3 = p^3+1, w2 = p^2+1
	mpz_mul_si (w2, w2, a[0]);  mpz_add (w1, w2, w3);  mpz_add (w3, w3, w3);				// w1 = (p^3+1) + (p^2+1)a1, w3 = 2(p^3+1)
	mpz_set_ui (w2, p+1);  mpz_mul_si (w2, w2, a2Min);  mpz_add (w1, w1, w2);			// w1 = (p^3+1) + (p^2+1)a1 + (p+1)*a2Min
	mpz_set_ui (e[1], 0);  good_a2 = 0;
	_jac_random(g, c[0]);
	hc_poly_twist (ct, c);
	_jac_random(gt, ct[0]);
	mpz_set_ui (w2, p+1);  mpz_mul_ui (w2, w2, p);									// w2 = p^2+p
//gmp_printf ("p=%ld, a2min=%ld, a2max=%ld, a2 mod p = %ld, a1=%ld, w1=%Zd\n", p, a2Min, a2Max, a[1], a[0], w1);
	for ( a2 = a2Min ; a2 <= a2Max ; a2 += p ) {
//printf ("p=%ld, a2=%ld\n", p, a2);
		a3Min = a3min (a[0], a2, p); t = a3Min % p;  if ( t < 0 ) t += p;  t = a3Min - t + a[2]; if ( t < a3Min ) t += p;  a3Min = t;
		a3Max = a3max (a[0], a2, p); t = a3Max % p;  if ( t < 0 ) t += p; t = a3Max - t + a[2];  if ( t > a3Max ) t -= p;  a3Max = t;
//printf ("p=%ld, a3min=%ld, a3max=%ld, a3 mod p = %ld\n", p, a3Min, a3Max, a[2]);
		if ( a3Min > a3Max ) { mpz_add (w1, w1, w2); continue; }
		mpz_add_si (Min, w1, a3Min);  mpz_add_si (Max, w1, a3Max);
//gmp_printf ("initial p=%ld, r=%d, w1=%Zd, w2=%Zd, Min = %Zd, Max = %Zd\n", p, r, w1, w2, Min, Max);
		while ( ! mpz_divisible_2exp_p (Min, r) ) mpz_add_ui (Min, Min, p);
		while ( ! mpz_divisible_2exp_p (Max, r) ) mpz_sub_ui (Max, Max, p);
		if ( mpz_cmp (Min, Max) > 0 ) { mpz_add (w1, w1, w2); continue; }
//gmp_printf ("p=%ld, r=%d, Min = %Zd, Max = %Zd\n", p, r, Min, Max);
		i = jac_search (e, &g, p*(1L<<r), Min, Max, 0, c);
/*
if ( ! i ) {
	mpz_set_ui (w5, p*(1L<<r));
	for ( mpz_set (w6, Min) ; mpz_cmp (w6, Max) <= 0 ; mpz_add (w6, w6, w5) ) {
		jac_exp_mpz (&h, &g, w6, c);
		if ( _jac_is_identity(h) ) gmp_printf ("*** p=%ld, jac_search missed exponent %Zd ***\n", p, w6);
	}
}
*/
		if ( i ) {
			mpz_set_ui (w4, p+1);  mpz_mul_si (w4, w4, a2);  mpz_add (w4, w4, w4);  mpz_add (w4, w4, w3);	// w4 = 2(p^3+1+(p+1)a2)
			if ( i == 1) {
				mpz_sub (w5, w4, e[0]);
				jac_exp_mpz (&h, &gt, w5, ct);
				if ( _jac_is_identity (h) ) {
					if ( cnt ) {
						mpz_set (N[1], e[0]);
						cnt = smalljac_distinguish_group_orders (N, 2, c);
						assert ( cnt <= 1 );
						if ( cnt && mpz_cmp(N[0],e[0]) == 0 ) good_a2 = a2;
					} else {
						mpz_set (N[0], e[0]);
						good_a2 = a2;
						cnt = 1;
					}
				}
			} else {
				// this case should be very rare, don't worry about handling it particularly efficiently
				mpz_sub (w5, e[1], e[0]);  mpz_sub (w6, Max, Min);
				mpz_fdiv_q (w6, w6, w5);  k = mpz_get_ui (w6)+2;
				M = malloc (k*sizeof(*M));
				for ( i = 0 ; i < k ; i++ ) mpz_init (M[i]);
				mpz_init_set (M[0], e[0]);
				mpz_init_set (M[1], e[1]);
				mpz_set (w5, e[1]);
				for ( i = 2 ; i < k ; i++ ) {
					mpz_add (w5, w5, w4);
					if ( mpz_cmp (w5, Max) > 0 ) break;
					mpz_init_set (M[i], w5);
				}
				k = i;
//gmp_printf ("%9lu: Found %d candidate orders (cnt=%d), checking twist compatibility\n", p, k, cnt);
				for ( i = 0 ; i < k ; i++ ) {
					mpz_sub (w5, w4, M[i]);
					jac_exp_mpz (&h, &gt, w5, ct);
					if ( ! _jac_is_identity (h) ) mpz_set_ui (M[i], 0);
				}
				for ( i = j = 0 ; j < k ; j++ ) if ( mpz_sgn(M[j]) ) mpz_set (M[i++], M[j]);
				if ( i ) {
//					gmp_printf ("%9lu: Found %d twist-compatible orders\n", p, i);
					if ( cnt ) mpz_set (M[i++], N[0]);
					if ( i > 1 ) {
						printf ("Distinguishing %d candidate orders\n", i);
						i = smalljac_distinguish_group_orders (M, i, c);
						printf ("Eliminated all but %d candidates\n", i);
						assert ( i <= 1 );
						if ( i == 1 && cnt && mpz_cmp(N[0],M[0]) != 0 ) { mpz_set (N[0], M[0]); good_a2 = a2; }
						cnt = i;
					} else {
						mpz_set (N[0], M[0]);  good_a2 = a2;
					}
				}
				for ( i = 0 ; i < k ; i++ ) mpz_clear (M[i]);
				free (M);
			}
		}
		mpz_add (w1, w1, w2);
	}
	if ( !cnt ) {
		gmp_err_printf ("%9lu: Search failed in genus3_charpoly_from_Pmodp a2Min = %ld, a2Max = %ld, Min = %Zd, Max = %Zd\n", _ff_p, a2Min, a2Max, Min, Max);
		gmp_err_printf ("%9lu: padic_lpoly coefficients: a1=%ld, a2=%ld, a3=%ld\n", _ff_p, a[0], a[1], a[2]);
		goto done;
	}
	a[1] = good_a2;
	mpz_set_ui (w1, p);  mpz_mul_ui (w1, w1, p);  mpz_set (w2, w1); mpz_mul_ui (w1, w1, p);  mpz_add_ui (w1, w1, 1);  mpz_add_ui (w2, w2, 1);		// w1 = p^3+1, w2 = p^2+1
	mpz_mul_si (w2, w2, a[0]);  mpz_add (w1, w1, w2);								// w1 = (p^3+1) + (p^2+1)a1
	mpz_set_ui (w2, p+1);  mpz_mul_si (w2, w2, a[1]);  mpz_add (w1, w1, w2);				// w1 = (p^3+1) + (p^2+1)a1 + (p+1)a2
	mpz_sub (N[0], N[0], w1);
	a[2] = mpz_get_si (N[0]);
	sts = 1;
done:
	smalljac_charpoly_gops += jac_gops-n;
	mpz_clear (e[0]); mpz_clear (e[1]); mpz_clear (w1); mpz_clear(w2); mpz_clear(w3); mpz_clear(w4); mpz_clear(w5); mpz_clear(w6); mpz_clear (Min); mpz_clear (Max); mpz_clear (N[0]); mpz_clear (N[1]);
	return sts;
}

// Distinguish among k > 1 possible group orders
// Modifies the array N so that on return N[0] holds distinguished group order (if there is one)
// Returns number of candidate orders remaining (this should be 1 but could be 0)
int smalljac_distinguish_group_orders (mpz_t N[], int k, hc_poly c[1])
{
	mpz_t X, Y;
	unsigned long q[MAX_UI_PP_FACTORS], e[MAX_UI_PP_FACTORS], o[JAC_MAX_GENERATORS], E, *m, u, v;
	jac_t b[JAC_MAX_GENERATORS];	// note that max generators is 2g+1, not 2g, jac_sylow requires room for one extra element
	int r, s, t, i, j, j2;
	
	assert ( k > 1 );

	mpz_init (X); mpz_init (Y);
	m = malloc (k*sizeof(*m));
	
	assert ( mpz_sgn(N[0]) && mpz_sgn(N[1]) );
	mpz_gcd (X, N[0], N[1]);
	for ( i = 2 ; i < k ; i++ ) { assert (mpz_sgn(N[i])); mpz_gcd (X, X, N[i]); }
	if ( ! mpz_fits_ulong_p (X) ) { gmp_printf ("GCD %Zd too large! k=%d, N[0]=%Zd, N[1]=%Zd\n", X, k, N[0], N[1]); abort(); }
	E = mpz_get_ui (X);
	t = ui_factor (q, e, E);
	for ( s = 0 ; s < t && k > 1 ; s++ ) {
		mpz_set_ui (Y, q[s]);
		for ( j = i = 0 ; i < k ; i++ ) { m[i] = mpz_remove (X, N[i], Y);  if ( m[i] > m[j] ) j = i; }		// find order with largest q[s]-valuation
		for ( j2 = j ? 0 : 1, i = 0 ; i < k ; i++ ) { if ( i == j ) continue; if ( m[i] > m[j2] ) j2 = i; }	// find order with second largest q[s]-valuation
		if ( j == j2 ) continue;														// no point if these are the same
		// rather than using the largest possible size, use q[s]*second largest, since this is enough to distinguish
		for ( v = q[s], i = 0 ; i < m[j2] ; i++ ) v *= q[s];
//printf ("Checking %ld-Sylow\n", q[s]);
		r = jac_sylow (b, o, q[s], E, v, 0, c);
		for ( u = o[0], i = 1 ; i < r ; i++ ) u *= o[i];
//printf ("Rank %d, size is %lu vs %lu\n", r, u, v);
		if ( u < v ) continue;
		// we found the largest possilbe p-Sylow subgroup, so we know we can eliminate at least one case (possibly more)
		for ( i = 0 ; i < k ; i++ ) if ( m[i] < m[j] ) m[i] = 0;
		for ( i = j = 0 ; i < k ; i++ ) if ( m[i] ) { mpz_set (N[j], N[i]); m[j] = m[i]; j++; }
		k = j;
	}
	return k;
}

// Computes Lpoly coefficients for curves over Q that are twists C' of  either C: y^2=x^6-5x^4-5x^2+1 orC: y^2=x^6+1 or , both of which have Jacobians
// that are Q-isogenous to the square of a CM elliptic curve (wich CM discriminants -2 and -3, respectively)
// Relies on Fite-Sutherland 2014, which implies that, up to a possible sign ambiguity on a1',
// the L-poly coefficients of C' are determined by those of C plus the pair (d_L, d_K), the residue degrees of p in L and K respectively.
// Since we don't necessarily know L and K, we don't try to determine d_L and d_K, rather we try all the possible cases (of which there are only 8).
int smalljac_genus2_charpoly_cmsquare (long a[], int D, hc_poly c[1])
{
	static mpz_t X, Y, N[8];
	static int init;
	unsigned long q[MAX_UI_PP_FACTORS], e[MAX_UI_PP_FACTORS], o[JAC_MAX_GENERATORS], u, v, E;
	jac_t g, h, g0, b[JAC_MAX_GENERATORS];	// note that max generators is 2g+1, not 2g, jac_sylow requires room for one extra element
	ff_t w;
	long a1[8], a2[8];
	long p, x, y, z;
	int m[8];
	int i, j, j2, k, r, s, t, n;
	
	p = (long) _ff_p;
	assert ( p > 3 && p < (1L<<40) );		// for convenience, we asssume 8*p^(3/2) < 2^63
	
	if ( ! init ) { mpz_init (X);  mpz_init (Y); for ( i = 0 ; i < 8 ; i++ ) mpz_init (N[i]); init = 1; }
	
	if ( D == -2 ) {
		if ( smalljac_cm2square_Lpoly (a) != 2 ) { err_printf ("smalljac_cm2square_Lpoly failed in smalljac_genus2_charpoly_cmsquare at p=%ld\n", _ff_p); return 0; }
	} else if ( D == -3 ) {
		_ff_set_one(w);
		if ( smalljac_x6pa_Lpoly (a, w) != 2 ) { err_printf ("smalljac_x6pa_Lpoly failed in smalljac_genus2_charpoly_cmsquare at p=%ld\n", _ff_p); return 0; }
	} else {
		err_printf ("smalljac_genus2_charpoly_cmsquare only supports CM discriminants -2 and -3\n"); abort(); 
	}
	
	if ( (D==-2 && (p&7)>3) || (D==-3 && ! _ff_p1mod3) ) {		// supersingular case
		assert (!a[0] && !(a[1]%p));
		a1[0] = 0;  a2[0] = 2*p;				// (2,2)
		a1[1] = 0;  a2[1] = -2*p;				// (4,2)
		a1[2] = 0;  a2[2] = p;				// (12,6)
		if ( D == -2 ) {
			a1[3] = 0;  a2[3] = 0;			// (8,4)
		} else {
			a1[3] = 0;  a2[3] = -p;			// (6,6)
		}
		k = 4;
	} else {
		assert (a[0] && !(a[0]&1));					// in the non-supersingular case a[0] is nonzero and even, this guarantees the 8 possiblities below are well-defined and distinct
		a1[0] = a[0]; 			a2[0] = a[1];			// (1,1)
		a1[1] = -a[0];			a2[1] = a[1];			// (2,1)
		a1[2] = 0;			a2[2] = 4*p-a[1];		// (2,2)
		a1[3] = -a[0]/2;		a2[3] = a[1]-3*p;		// (3,3)
		a1[4] = 0;			a2[4] = a[1]-4*p;		// (4,2)  we don't currently use the fact that case (2,2) and (4,2) differ only in sign
		a1[5] = a[0]/2;			a2[5] = a[1]-3*p;		// (6,3)
		if ( D == -2 ) {
			y = 2*(4*p-(a[0]*a[0])/4);
			x = (long) sqrt(y);
			if ( x*x != y ) { err_printf ("2*(4*p-a1^2) is not a perfect square for y^2=x^6-5x^4-5x^2+1 mod %ld\n", p);  abort(); }
			a1[6] = x;			a2[6] = 6*p-a[1];		// (8,4)
			a1[7] = -x;			a2[7] = 6*p-a[1];		// (8,4)			
		} else {
			y = 3*(4*p-(a[0]*a[0])/4);
			x = (long) sqrt(y);
			if ( x*x != y ) { err_printf ("3*(4*p-a1^2) is not a perfect square for y^2=x^6+1 mod %ld\n", p); abort(); }
			a1[6] = x;			a2[6] = 7*p-a[1];		// (6,6)
			a1[7] = -x;			a2[7] = 7*p-a[1];		// (6,6)
		}
		k = 8;
	}
	
	// eliminate pairs that are incompatible with 2-torsion -- this appears to never leave more than 4 possibilities (and sometimes only 1)
	j = (1 << hc_poly_two_rank (c)) - 1;
	x = p&0xF;  x = (x*x)&0xF;  x = (x+1)&0xF;	// x = p^2+1 mod 16
	y = (p+1)&0xF;						// y = p+1 mod 16
	for ( i = 0 ; i < k ; i++ ) {
		z = (x + y*a1[i] + a2[i]) &0xF;			// z = p^2 + a1*p + a2 + a1 + 1 mod 16
		m[i] = 1;
		if ( (z&j) ) m[i] = 0;					// group order must be 0 mod 2^two-rank
		if ( ! j && !(z&1) ) m[i] = 0;			// if two-rank is zero then the group order must be odd
	}
	for ( i = j = 0 ; i < k ; i++ ) if ( m[i] ) { a1[j] = a1[i];  a2[j] = a2[i]; m[j] = m[i]; j++; }
	k = j;
	if ( ! k ) {err_printf ("smalljac_genus2_charpoly_cmsquare failed with p=%ld, a1=%ld, a2=%ld because none of the cases applied!\n", p, a[0], a[1]); abort(); }
	if ( k == 1 ) { a[0] = a1[0]; a[1] = a2[0]; return 2; }

	// now we need to distinguish the reminaing possible cases
	n = jac_gops;
	for ( int r = 0 ; r < SMALLJAC_RETRIES ; r++ ) {
		do { _jac_random(g, c[0]); } while ( _jac_is_identity(g) );
		_jac_set (g0, g);
		jac_exp_ui (&h, &g, p, c);
		jac_exp_ui (&g, &h, p, c);
		_jac_mult (g, g, g0, c[0]);		// h = (p^2+1)*g
		_jac_invert (g);
		for ( i = j = 0 ; i < k ; i++ ) {
			if ( ! m[i] ) continue;
			jac_exp_si (&h, &g0, (p+1)*a1[i]+a2[i], c);
			if ( _jac_cmp (g,h) != 1 ) m[i] = 0;
		}
		for ( i = j = 0 ; i < k ; i++ ) if ( m[i] ) { a1[j] = a1[i];  a2[j] = a2[i]; m[j] = m[i]; j++; }
		k = j;
		if ( ! k ) {err_printf ("smalljac_genus2_charpoly_cmsquare failed with p=%ld, a1=%ld, a2=%ld because none of the cases applied!\n", p, a[0], a[1]); abort(); }
		if ( k > 1 ) continue;
		a[0] = a1[0];  a[1] = a2[0];
		smalljac_charpoly_gops += jac_gops-n;
		return 2;
	}

	// if we get here we assume that the group exponent divides 2 of the remaining candidate orders, now we compute sylow subgroups (using dlogs) to distinguish
	// the good news is that this likely means that there is a small p for which knowing the cardinality of the p-Sylow subgroup will determine which case we are in
	mpz_set_ui (X, p);  mpz_mul_ui (X, X, p);  mpz_add_ui (X, X, 1);
	for ( i = 0 ; i < k ; i++ ) { mpz_set_ui (Y, p+1); mpz_mul_si (Y, Y, a1[i]);  mpz_add_si (Y, Y, a2[i]); mpz_add (N[i], X, Y); }
	mpz_gcd (X, N[0], N[1]);
	for ( i = 2 ; i < k ; i++ ) mpz_gcd (X, X, N[i]);
	assert ( mpz_fits_ulong_p (X) );
	E = mpz_get_ui (X);
	t = ui_factor (q, e, E);
	for ( s = 0 ; s < t && k > 1 ; s++ ) {
		mpz_set_ui (Y, q[s]);
		for ( j = i = 0 ; i < k ; i++ ) { m[i] = mpz_remove (X, N[i], Y);  if ( m[i] > m[j] ) j = i; }		// find order with largest q[s]-valuation
		for ( j2 = j ? 0 : 1, i = 0 ; i < k ; i++ ) { if ( i == j ) continue; if ( m[i] > m[j2] ) j2 = i; }	// find order with second largest q[s]-valuation
		if ( j == j2 ) continue;														// no point if these are the same
		// rather than using the largest possible size, use q[s]*second largest, since this is enough to distinguish
		for ( v = q[s], i = 0 ; i < m[j2] ; i++ ) v *= q[s];
//printf ("Checking %ld-Sylow\n", q[s]);
		r = jac_sylow (b, o, q[s], E, v, 0, c);
		for ( u = o[0], i = 1 ; i < r ; i++ ) u *= o[i];
//printf ("Rank %d, size is %lu vs %lu\n", r, u, v);
		if ( u < v ) continue;
		// we found the largest possilbe p-Sylow subgroup, so we know we can eliminate at least one case
		for ( i = 0 ; i < k ; i++ ) if ( m[i] < m[j] ) m[i] = 0;
		for ( i = j = 0 ; i < k ; i++ ) if ( m[i] ) { a1[j] = a1[i];  a2[j] = a2[i]; m[j] = m[i]; j++; }
		k = j;
	}

	if ( ! k ) {err_printf ("smalljac_genus2_charpoly_cmsquare failed with p=%ld, a1=%ld, a2=%ld because none of the cases applied!\n", p, a[0], a[1]); abort(); }
	if ( k > 1 ) {
		// this should never happen!
		err_printf ("Unable to resolve ambiguity in smalljac_genus2_charpoly_cmsquare with p=%ld\n", p);
		ff_poly_print (c->f, c->d);
		for ( i = 0 ; i < k ; i++ ) printf ("a1=%ld, a2=%ld, N=%ld\n", a1[i], a2[i], (p*p+1) + (p+1)*a1[i] + a2[i]);
		return 0;
	}
	a[0] = a1[0];  a[1] = a2[0];
	smalljac_charpoly_gops += jac_gops-n;
	return 2;
}


/*
	Attempt to compute P(T) given the value of |a2| mod p and the parity of a1 (in genus 2)
	Uses O(p^(1/4)) group operations.

	Not used.
*/
int smalljac_genus2_charpoly_from_a2_modp (long a[2], hc_poly c[1])
{
	static mpz_t E[2], M[2];
	static int init;
	jac_t *new_a, r, g, h, w[9], wi[9];
	unsigned long n;
	double x;
	long p, e, e1, e2, Min, Max, t, a1min, a1max, a1, a2;
	int i, j, k, wcnt;
	
	if ( ! init ) { mpz_init(E[0]); mpz_init(E[1]); mpz_init(M[0]); mpz_init(M[1]); init = 1; }
	
	p = (long) _ff_p;

	n = jac_gops;
	a[1] = i_mod (a[1], p);
	wcnt = ( a[1] ? 8 : 9 );
//printf ("a2_modp input a1=%ld, a2=%ld mod %ld\n", a[0], a[1], _ff_p);
	do { _jac_random(r, c[0]); } while ( _jac_is_identity(r) );
//printf ("random point r is "); _jac_print(r);

	// compute slightly constrained bounds on a1 given a1, but just take the largest interval that contains all the possiblities, rather than customizing the interval for a1 for each possible a2			
	x = sqrt(p);
	a1min = a1max = 0;
	for ( k = 0 ; k < wcnt ; k++ ) {
		a2 = a[1]+(k-2)*p;
		a1 = -(long)((a2/x+2*x)/2);
		if ( a1 < a1min ) a1min = a1;
		a1 = (long)((a2/x+2*x)/2);
		if ( a1 > a1max ) a1max = a1;
		a2 = (p-a[1])+(k-2)*p;
		a1 = -(long)((a2/x+2*x)/2);
		if ( a1 < a1min ) a1min = a1;
		a1 = (long)((a2/x+2*x)/2);
		if ( a1 > a1max ) a1max = a1;
	}
	// make sure a1 bounds have the right parity
	if ( (a1min-a[0])&1 ) a1min++;
	if ( (a1max-a[0])&1 ) a1max--;
	// make bounds symmetric
	if ( -a1min > a1max ) a1max = -a1min;
	if ( a1min > -a1max ) a1min = -a1max;

	e1 = e2 = 0;
	// crossover between O(p^(1/2)) algorithm and O(p^(1/4)) algorithm is around p=10^6
	if ( p < 1000000 ) {
	
		// the order of r is p^2+1+(p+1)*a0+/-a[1]+k*p for some k such that (a[i]+k*p) is in [-2p,6p]
		// there are at most 9 (usually 8) possiblities for k, and we know that |a0| <= 4*sqrt{p}
		jac_exp_ui (&g, &r, p, c);	
		jac_exp_ui (w+2, &g, p, c);
		_jac_mult (w[2], w[2], r, c[0]);

		jac_exp_ui (&h, &r, a[1], c);
		_jac_set (wi[2],w[2]);
		_jac_mult (w[2], w[2], h, c[0]);
		if ( a[1] ) {
			_jac_mult (wi[2], wi[2], g, c[0]);
			_jac_invert (h);
			_jac_mult (wi[2], wi[2], h, c[0]);
		}
		for ( k = 3 ; k < wcnt ; k++ ) { _jac_mult (w[k], w[k-1], g, c[0]); if ( a[1] ) _jac_mult (wi[k], wi[k-1], g, c[0]); }
		_jac_invert (g);
		for ( k = 1 ; k >= 0 ; k-- ) { _jac_mult (w[k], w[k+1], g, c[0]); if ( a[1] ) _jac_mult (wi[k], wi[k+1], g, c[0]); }
		_jac_invert(r);
		_jac_mult(g,g,r,c[0]);

		// now w[k] holds r^(p^2+1+a[1]+k*p), wi[k] holds r^(p^2+p-a[1]+k*p) and g holds r^{-(p+1)}
		// we want to compute g^a0 for all |a0| <= 4*sqrt(p) and check whether g^a0 = w[k] for some k
		if ( ! (a[0]&1) ) {
			for ( k = 0 ; k < wcnt ; k++ ) {
				if ( _jac_is_identity(w[k]) ) { e=(p*p+1)+a[1]+(k-2)*p; if ( e != e1 ) { e2=e1; e1=e; } if ( e2 ) break; }
				if ( a[1] && _jac_is_identity(wi[k]) ) { e=(p*p+1)+(p-a[1])+(k-2)*p;  if ( e != e1 ) { e2=e1; e1=e; } if ( e2 ) break; }
			}
		}
		if ( ! e2 ) {
			if ( (a[0]&1) ) { i = 1; } else { _jac_square (g, g, c[0]); i = 2; }
			_jac_set (h,g);  
			for ( ; i <= a1max ; i+=2 ) {
				for ( k = 0 ; k < wcnt ; k++ ) {
					if ( (j=_jac_cmp(h,w[k])) ) { e=(p*p+1)+a[1]+(k-2)*p+i*j*(p+1); if ( e != e1 ) { e2=e1; e1=e; } if ( e2 ) break; }
					if ( a[1] && (j=_jac_cmp(h,wi[k])) ) { e=(p*p+1)+(p-a[1])+(k-2)*p+i*j*(p+1); if ( e != e1 ) { e2=e1; e1=e; } if ( e2 ) break; }
				}
				if ( e2 ) break;
				_jac_mult (h, h, g, c[0]);
			}
		}
	} else {
		new_a = &r;
		for ( k = 0 ; k < wcnt ; k++ ) {
			a2 = a[1]+(k-2)*p;
			mpz_set_ui (M[0], (p*p+1)+a2 + a1min*(p+1));
			mpz_set_ui (M[1], (p*p+1)+a2 + a1max*(p+1));
			i = jac_search (E, new_a, 2*(p+1), M[0], M[1], 2*wcnt, c);
			if ( i ) { e = mpz_get_ui(E[0]); if ( e != e1 ) { e2 = e1; e1 = e; }  if ( i > 1 ) { e2 = e1; e1 = mpz_get_ui(E[1]); } }
			if ( e2 ) break;
			new_a = 0;
			a2 = (p-a[1])+(k-2)*p;
			mpz_set_ui (M[0], (p*p+1)+a2 + a1min*(p+1));
			mpz_set_ui (M[1], (p*p+1)+a2 + a1max*(p+1));
			i = jac_search (E, new_a, 2*(p+1), M[0], M[1], 2*wcnt, c);
			if ( i ) { e = mpz_get_ui(E[0]); if ( e != e1 ) { e2 = e1; e1 = e; }  if ( i > 1 ) { e2 = e1; e1 = mpz_get_ui(E[1]); } }
			if ( e2 ) break;
		}
	}

	if ( ! e1 ) { err_printf ("smalljac_genus2_charpoly_from_a2_modp failed, a2 cannot be congruent to +/-%ld mod %ld\n", a[1], p); return 0; }
	if ( e2 ) {
//		err_printf ("Ambiguous result in smalljac_genus2_charpoly_from_a2_modp with e1=%ld, e2=%ld, p=%ld\n", e1, e2, p);
		e1 -= e2;
		if ( e1 < 0 ) e1 = -e1;
		if ( ! jac_fastorder_ui ((unsigned long*)&e2, &r, e1, c) ) { err_printf ("jac_fastorder_ui failed in smalljac_genus2_charpoly_from_a2_modp\n"); return 0; }
		Min = (unsigned long) floor(pow(sqrt(p)-1.0, 4.0));
		Max = (unsigned long) ceil(pow(sqrt(p)+1.0, 4.0));
		k = jac_order ((unsigned long*)&e1, Min, Max, JAC_INVALID_A1, e2, 0, 0, c);
		if ( k <= 0 ) { err_printf ("jac_order failed in smalljac_genus2_charpoly_from_a2_modp\n"); return 0; }
		if ( k > 1 ) { err_printf ("jac_order returned ambiguous result in smalljac_genus2_charpoly_from_a2_modp\n"); return 0; }
	}
	smalljac_charpoly_gops += jac_gops-n;
	// now try to determine a1 from Lp(1) and a2 mod p
	a[0] = JAC_INVALID_A1;
	for ( k = 0 ; k < wcnt ; k++ ) {
		t = e1-(p*p+1)-(a[1]+(k-2)*p);
		if ( ! (t%(p+1)) && (t/(p+1))*(t/(p+1)) <= 17*p ) { if ( a[0] != JAC_INVALID_A1 ) break; else a[0] = t/(p+1); }
		t = e1-(p*p+1)-((p-a[1])+(k-2)*p);
		if ( ! (t%(p+1)) && (t/(p+1))*(t/(p+1)) <= 17*p ) { if ( a[0] != JAC_INVALID_A1 ) break; else a[0] = t/(p+1); }
	}
	if ( k < wcnt ) { err_printf("Lpoly not determined by group order and a2 mod p, reverting to standard computation\n"); return 0; }
	a[1] = e1 - ((p*p+1)+a[0]*(p+1));
//printf ("success, a1=%ld a2=%ld mod %ld\n", a[0], a[1], p);
	return 1;
}
