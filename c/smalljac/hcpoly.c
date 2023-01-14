#include <stdio.h>
#include <gmp.h>
#include "mpzutil.h"
#include "hcpoly.h"

/*
    Copyright (c) 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// given a hyperelliptic curve y^2=F(x) with F  in Q[x], constructs an isomorphic curve y^2=f(x) with f  in Q[x]
// set den to the lcm of the denominators in the coefficients of F, times the leading coefficient of F
// the curve is not defined at primes dividing den
void hc_integral_model (mpz_t den, mpz_t f[], mpq_t F[], int d)
{
	mpz_t t1, t2, t3;
	long p;
	int i, j;
	
	mpq_get_den (den, F[0]);
	for ( i = 1 ; i <= d ; i++ ) if ( mpq_sgn(F[i]) && mpz_cmp_ui(mpq_denref(F[i]),1) != 0 ) mpz_lcm (den, den, mpq_denref(F[i]));

	// handle integer case quickly
	if ( mpz_cmp_ui(den,1) == 0 ) {
		for ( i = 0 ; i <= d ; i++ ) mpz_set (f[i], mpq_numref(F[i]));
		if ( !(d&1) || mpz_cmp_ui(f[d],1) == 0 ) return;
		mpz_init (t1);
		goto make_monic;
	}

	mpz_init (t1); mpz_init (t2);  mpz_init (t3);
	// first make leading coefficient integral by multiplying be the square of its denominator
	
	// explicitly handle powers of primes <= d which often appear in denominators
	mpz_set (t1, den);  mpz_set_ui (t3, 1);
	for ( i = 1 ; (p=ui_small_prime(i)) <= d ; i++ ) {
		mpz_set_ui (t2, p);
		j = mpz_remove (t1, t1, t2);
		mpz_ui_pow_ui (t2, p, j+(j&1));
		mpz_mul (t3, t3, t2);
	}
	mpz_mul (t1, t1, t1);  mpz_mul (t1, t1, t3);		// t1 is now a perfect square that is a multiple of den
	
	// multiply by t1 to clear all denominators
	for ( i = 0 ; i <= d ; i++ ) { mpz_divexact (t2, t1, mpq_denref(F[i])); mpz_mul (f[i], t2, mpq_numref(F[i])); }

	mpz_clear (t2);  mpz_clear (t3);

make_monic:
		
	// if d is odd make f monic -- incorporate lc(f) into den to note bad reduction (for d odd, y^2=cx^d+... necessarily has bad reduction at primes dividing c)
	if ( (d&1) ) {		
		mpz_set (t1, f[d]);
		mpz_lcm (den, den, t1);
		// replace f by t1^(d-1)*f(x/t1)
		for ( i = d-2 ; i>= 0 ; i-- ) { mpz_mul (f[i], f[i], t1); mpz_mul (t1, t1, f[d]); }
		mpz_set_ui (f[d], 1);
	}
	
	mpz_clear (t1);
}

// k[0] counts irred factors over Fp, k[i] counts factors of degree i for i in [1..d] (so k[0] = sum_{1<=i<=d}k[i])
int hc_two_rank_from_factorization (int d, int k[HECURVE_MAX_DEGREE+1])
{
	int g;
	
	// the odd degree case y^2=f(x) is easy: there is a 1-1 correspondence between divisors of degree <= g of f(x) (including the trivial divisor), and this is exactly half the total number 2^k[0] of divisors of f(x)
	if ( (d&1) ) return k[0]-1;
	g = (d-1)/2;
	// Important: we assume that for even degree curves y^2=f(x) the poly f has no roots over Fp (if it did, we would have moved this point to infinity to make the degree of f odd)
	assert (!k[1]);
	switch (g) {
	case 1: return 0;
	case 2:	// given that there are no linear factors, the number of 2-torsion points is just the number of divisors of degree 2 (0,1,or 3), plus 1 for the trivial divisor
		if ( ! k[2] ) return 0;
		if ( k[2] == 1 ) return 1;
		return 2;
	}
	// we could easily implement the genus 4 case but we don't currently need it.  The genus 3 case is actually a bit of a pain when there are no rational points at infinity, so we leave it out for the moment.
	printf ("two-rank computation for even degree curves of genus > 2 not implemented in hpoly.c\n");
	abort();
}

int hc_poly_compute_two_rank (hc_poly *c)
{
	int k[HECURVE_MAX_DEGREE+1];
	
	assert (c->n == 1);
	k[0] = ff_poly_factorization_pattern (k, c->f, c->d);
	c->two_rank = hc_two_rank_from_factorization (c->d, k);
	return c->two_rank;
}


/*
    For curves of the form y^2=f(x) with deg f even,  we translate the curve as required to ensure
    that there are no rational pts at infinity (per the suggestion of  Galbraith-Harrison-Morales in
    their ANTS VIII paper).  Currently this is only relevant in genus 2, but the code should work for even degree genus in general

    For curves of the form y^2=f(x) with deg f odd, we make the curve monic and kill the 2g coefficient.
    This assumes that we are working over a field whose characteristic does not divide deg f.
*/
void hc_poly_standardize (hc_poly *c)
{
	int k[HECURVE_MAX_DEGREE+1];
	ff_t s, t;
	int i;
	
	//printf ("hc_poly_standardize p=%ld, d=%d, f[d]=%ld, f[d-1]=%ld, n=%d\n", _ff_p, c->d, _ff_get_ui(c->f[c->d]), _ff_get_ui(c->f[c->d-1]), c->n);
	if ( c->n > 1 ) {  err_printf ("Internal error: hc_poly_standardize called with curve of degree %d defined over a number field of degree n=%d > 1, this is not currently supported\n", c->d, c->n); abort(); }
	
	if ( !(c->d&1) ) {
		if ( _ff_zero(c->f[c->d]) ) {
			c->d--; 	// if leading coefficient is 0 mod p (this can happen with a nonsingular curve y^2=deg 2g+2), revert to odd degree case (handled below)
		} else {
			// if f has a rational root, translate so (0,0) is on the curve and then move this point to infty to get an odd degree poly
			k[0] = ff_poly_factorization_pattern_and_root (k, c->f, c->d, &t);
			if ( k[1] ) {
				ff_poly_translate (c->f,0,c->f,c->d,t);
				ff_poly_reverse_inplace (c->f,c->d);
				if ( ! _ff_zero(c->f[c->d]) ) { err_printf ("error in smalljac_adjust_hecurve, translating root did not kill the constant term\n"); abort(); }
				if ( _ff_zero(c->f[c->d-1]) ) { err_printf ("singular curve in hc_poly_standardize at p=%ld\n", _ff_p); abort(); }
				c->d--;
				k[0]--; k[1]--;
			}
			// go ahead and set the two rank since we have the factorization
			c->two_rank = hc_two_rank_from_factorization (c->d, k);
		}
	}
	// if f now has odd degree we need to make it monic and kill the 2g coefficient
	if ( (c->d&1) ) {
		if ( ! _ff_one (c->f[c->d]) ) {
			// make sure leading coefficient is a residue before making f monic
			if ( ! ff_residue(c->f[c->d]) ) {
				ff_nonresidue(&t);
				_ff_mult(c->f[1],c->f[1],t);
				_ff_set(s,t);
				for ( i = 2 ; i <= c->d ; i++ ) { _ff_mult(s,s,t); _ff_mult(c->f[i],c->f[i],s); }
			}
			ff_poly_monic(c->f,0,c->f,c->d);		// note that deg 2g+1 coefficient cannot also be zero (ow curve is singular)
		}
		if ( ! _ff_zero (c->f[c->d-1]) && (_ff_p%c->d)) {
			// translate so 2g coefficient is 0
			ff_invert_small_int(&t,c->d);
			_ff_mult(s,t,c->f[c->d-1]); _ff_neg(t,s);
			ff_poly_translate (c->f,0,c->f,c->d,t);
		}
		return;
	}
	// for even degree curves y^2=f(x) with no rational root, we want to make the leading coefficient a non-residue to ensure no points at infinity
	// To do this we translate so that the constant coefficient is a non-residue and then reverse the coefficients
	if ( ! ff_residue (c->f[c->d]) ) return;
	_ff_set_one(t);
	do {
		ff_poly_eval (&s, c->f, c->d, &t);
		if ( ! ff_residue (s) ) break;
		_ff_inc(t);
	} while ( ! _ff_zero(t) );
	if ( _ff_zero(t) ) { err_printf ("Fatal error, can't translate sextic to get 0 pts at infinity, field F_%ld is too small!\n", _ff_p);  abort(); }
	// now we want to translate by t and reverse the poly to make the constant coeff the leading coeff
	ff_poly_translate (c->f,0,c->f,c->d,t);
	ff_poly_reverse_inplace (c->f,c->d);
	// make leading term equal to our default non-residue (which will be small and have a known inverse)
	_ff_invert(s,c->f[c->d]);
	ff_nonresidue(&t);
	_ff_mult(s,s,t);
	_ff_set(c->f[c->d],t);
	for ( i = 0 ; i < c->d ; i++ ) _ff_mult(c->f[i],c->f[i],s);
	if ( ! _ff_zero (c->f[c->d-1]) && (_ff_p%c->d) ) {
		// translate by -f[d-1]/(d*f[d]) so that f[d-1] is zero
		ff_nonresidue_inverse(&s);
		ff_invert_small_int(&t,c->d);
		_ff_mult(t,t,s); _ff_mult(s,t,c->f[c->d-1]);  _ff_neg(t,s);
		ff_poly_translate (c->f,0,c->f,c->d,t);
	}
}
