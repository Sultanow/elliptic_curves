#include <stdio.h>
#include "ff_poly.h"
#include "smalljac_internal.h"

/*
    Copyright (c) 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

/*
	This module contains code to efficiently compute L-polynomials (or just a1) for special curves,
	currently restricted to hyperelliptic curves of the form y^2=x^d+a*x and y^2=x^d+a.
	(support for special Fermat curves (including Picard curves) would be fairly easy to add.
*/


/*
	Computes Lpoly coefficient a_1 = - tr(E/Fp) where E: y^2 = x^3 + f0 is defined over Q.
	Such a curve has j-invariant 0 and CM by the imaginary quadratic order of discriminant -3.
	It follows that for p = 1 mod 3 we have 4p = t^2 + 3v^2, where t = tr(E/Fp).
	To distinguish the 6 possible values of t, we effectively invert Rubin-Silverberg
	"Choosing the correct elliptic curve in the CM method", Alg 3.5 (but we simplify things a bit)

	Return value is 1 for success, -1 if the curve is singular, and -2 for any internal errors (which should be impossible)
*/
int smalljac_x3pa_Lpoly (long a[], ff_t f0)
{
	ff_t w, x, y, ss, tt, uu, p;
	long s, t, u, v;
	
	if ( _ff_zero (f0) ) return -1;						// singular curve
	
	p = _ff_p;
	
	if ( (p%3) == 2 ) { a[0] = 0; return 1; }			// if p = 2 mod 3 then E is supersingular and a_1=0

	// get a solution to p=u^2+3v^2, which necessarily exists and is unique up to signs
	if ( ! ff_cornacchia (&u, &v, 3) ) { printf ("internal error, %ld = 1 mod 3 is not a prime of the form p = u^2 + 3v^2 !!!\n", p);  return -2; }
	
	// make u = 1 mod 3 and make v = 0,2 mod 3
	if ( u%3 == 2 ) u = -u;
	if ( v%3 == 1 ) v = -v;
	
	ff_exp_ui (&w, &f0, (p-1)/6);						// w = b^n is a 6th root of unity, which must be one of 1,-1,+/-2u/(3v-u),+/-2u/(3v+u)
	
	// handle each of the 6 possibilities for w
	if ( _ff_one(w) ) { a[0] = -2*u; return 1; }			// if w = 1 than a1 = -2u (minus the trace)
	if ( _ff_neg_one(w) ) { a[0] = 2*u; return 1; }		// if w = -1 than a1 = 2u
	
	s = 3*v-u;  t = 3*v+u;							// s = 3v-u, t = 3v+u
	_ff_set_i (uu, u);  _ff_set_i (ss, s);  _ff_set_i (tt, t);
	_ff_mult (x, ss, tt); ff_invert (x, x);
	_ff_mult (x, x, uu); _ff_x2 (x); 
	_ff_mult (y, x, ss); _ff_mult (x, x, tt);				// x = 2u/(3v-u), y= 2u/(3v+u);
	if ( _ff_equal (w, x) ) {							// if w = 2u/(3v-u) then a1 = 3v+u			
		a[0] = t;
		return 1;
	} else if ( _ff_equal (w, y) ) {						// if w = 2u/(3v+u) then a1 = 3v-u
		a[0] = s;
		return 1;
	}
	ff_negate (w);
	if ( _ff_equal (w, x) ) {							// if w = -2u/(3v-u) then a1 = -(3v+u)
		a[0] = -t;
		return 1;
	} else if ( _ff_equal (w, y) ) {						// if w = -2u/(3v+u) then a1 = -(3v-u)
		a[0] = -s;
		return 1;
	}
	printf ("Fell through to unreachable code in smalljac_x3pa_Lpoly_Q with p=%ld, interntal error, u=%ld, v=%ld, -w=%ld.\n", p, u, v, _ff_get_ui(w)); return -2;
}



/*
	Computes Lpoly coefficient a_1 = - tr(E/Fp) where E: y^2 = x^3 + ax is defined over Q.
	Such a curve has j-invariant 1728 and CM by the imaginary quadratic order of discriminant -4.
	It follows that for p=1 mod4 we have 4p = t^2 + 4v^2, and therfore p=u^2+v^2 where 2u = t = tr(E/Fp). 
	To distinguish the 4 possible values for u, we invert Rubin-Silverberg "Choosing the correct elliptic curve in the CM method", Alg 3.4.

	eturn value is 1 for success, -1 if the curve is singular, and -2 for any internal errors (which should be impossible)
*/
int smalljac_x3pax_Lpoly (long a[], ff_t f1)
{
	long s, u, v, p;
	ff_t w, x, y;

	if ( _ff_zero(f1) ) return -1;					// singular curve
	
	p = _ff_p;
	
	if ( (p&3) == 3 ) { a[0] = 0; return 1; }		// if p = 3 mod 4 then E is supersingular and a_1=0
	
	// get a solution to p=u^2+v^2
	if ( ! ff_cornacchia (&u, &v, 1) ) { printf ("internal error, %ld = 1 mod 4 is not a prime of the form u^2 + v^2 !!!\n", p);  return -2; }
	if ( v&1 ) { s = u; u = v; v = s; }				// make u odd
	if ( (u-1-v)&3 ) u = -u;						// make u-1=v mod 4

	// apply Alg 3.4 of Rubin-Silverberg in reverse.
	_ff_neg (w, f1);								// negate f1 because RS use y^2 = x^3 - ax
	ff_exp_ui (&w, &w, (p-1)/4);					// w is a 4th root of unity, which is either 1,-1,u/v,or,-u/v
	
	// handle each of the 4 possibilities for w
	if ( _ff_one (w) ) { a[0] = -2*u; return 1; }		// w=1 => a1 = -2u (negative trace)
	if ( _ff_neg_one(w) ) { a[0] = 2*u; return 1; }	// w=-1 => a1 = 2u
	_ff_set_i (x, u);  _ff_set_i (y, v);
	_ff_invert (y, y); _ff_mult (x, x, y);				// compute u/v in Fp
	if ( _ff_equal(w, x) ) { a[0] = 2*v; return 1; }		// w=u/v => a1 = 2v
	ff_negate(w);
	if ( _ff_equal(w, x) ) { a[0] = -2*v; return 1; }	// w=-u/v => a1 = -2v
	printf ("Fell through to unreachable code in smalljac_x3pa_Lpoly_Q with p=%ld, internal error, u=%ld, v=%ld, -w=%ld.\n", p, u, v, _ff_get_ui(w)); return -2;
}


/*
	Computes Lpoly coefficients for curves of the form y^2 = x^5 + ax over Q, doesn't compute group structure

	Based on results in E. Furukawa, M. Kawazoe and T. Takahashi,
	"Counting Points for Hyperelliptic Curves of Type y2 = x^5 +ax over Finite Prime Fields"
	Selected Areas in Cryptography (SAC 2003), LNCS 3006, pp. 26-41, Springer, 2004,
	which are nicely summarized in Theorem 1 of Kawazoe and Takahasi,
	"Pairing-friendly hyperelliptic curves with ordinary Jacobians of type y^2 = x^5+ax",
	in Pairings 2008, LNCS 5209, pp. 164-177, Springer, 2008

	Note: a[] is the array of L-poly coefficients being computed, the value of the "a" coefficient in the
        curve equation is stored in sc->f[0]
*/
int smalljac_x5pax_Lpoly (long a[], ff_t f1)
{
	long x, y, p;
	ff_t s, t, w;

	if ( _ff_zero(f1) ) return -1;					// singular curve
	
	p = _ff_p;
	
	// handle p=7 mod 8 (easiest case) first
	if ( (p&7) == 7 ) { a[0] = 0;  a[1] = 2*p;  return 2; }		// if p=7 mod 8 then L_p(T) = (pT^2+1)^2.
	
	// handle p=5 mod 8 (next easiest case)
	if ( (p&7) == 5 ) {
		ff_exp_ui (&s, &f1, (p-1)/4);
		if ( _ff_one (s) ) { a[0] = 0;  a[1] = 2*p;  return 2; } 	// if p=5 mod 8 and u^((p-1)/4) = 1 (u is a 4th power) then L_p(T) = (pT^2+1)^2.
		ff_negate (s);
		if ( _ff_one (s) ) { a[0] = 0;  a[1] = -2*p; } 			// if p=5 mod 8 and u^((p-1)/4) = -1 (u is a square but not a 4th power) then L_p(T) = (pT^2-1)^2.
		else { a[0] = a[1] = 0; }							// o.w. p = 5 mod 8 and u is not a QR mod p, so L_p(T)=pT^4+1
		return 2;
	}
	
	// now we need to compute x=1 mod 4 s.t. p = x^2+2y^2, which must exists for p=1,3 mod 8
	if ( ! ff_cornacchia (&x, &y, 2) ) { printf ("internal error in smalljac_x5pax_Lpoly, %ld is not a prime of the form x^2 + 2y^2\n", p);  return -2; }
	if ( (x&3) != 1 ) x = -x;	// make sure x is 1 mod 4

	// handle p=3 mod 8 (next easiest case)
	if  ( (p&7) == 3 ) {
		ff_exp_ui (&s, &f1, (p-1)/2);
		if ( _ff_one(s) ) {
			a[0] = 0;  a[1] = 2*p-4*x*x;					// p = 3 mod 8 and u^((p-1)/2) = 1 (u is a square)
		} else {
			a[0] = 0;  a[1] = 4*x*x - 2*p;					// p = 3 mod 8 and u^((p-1)/2) != 1 (u is not a square)
		}
		return 2;
	}
	
	// handle p=1 mod 8 (hardest case)
	ff_exp_ui (&s, &f1, (p-1)/8);
	_ff_square (t, s);
	if ( _ff_one (t) ) {									// if u is a 4th power
		a[1] = 4*x*x + 2*p;
		if ( _ff_one (s) ) {
			if ( (p&0xF) == 1 ) {
				a[0] = -4*x;							// p = 1 mod 16 and u^((p-1)/8) = 1 (u is an eight power)
			} else {
				a[0] = 4*x;							// p = 9 mod 16 and u^((p-1)/8) = 1 (u is an eigth power)
			}
		} else {
			if( (p&0xF) == 1 ) {
				a[0] = 4*x;							// p = 1 mod 16 and u^((p-1)/8) = -1 (u is a fourth power but not an eigth power)
			} else {
				a[0] = -4*x;							// p = 9 mod 16 and u^((p-1)/8) = -1 (u is a fourth power but not an eigth power)
			}
		}
		return 2;
	}
	
	// p=1 mod 8 and u is not a 4th power
	ff_negate (t);
	if ( _ff_one (t) ) { a[0] = 0;  a[1] = 4*x*x - 2*p;  return 2; }	// u is a square but not a fourth power
	ff_negate (t);
	
	// p=1 mod 8 and u is not a square (last and hardest case)
	
	// adjust the sign of y if necessary, we require that 2*(-1)^e*y = (u^e + u^(3e))*x mod p, where e = (p-1)/8
	ff_mult (t, t, s);
	_ff_addto (t, s);  _ff_set_i (w, x);  ff_mult (t, t, w);
	_ff_set_i (w, 2*y);
	if ( (p&0xF) != 1 ) ff_negate (w);
	if ( ! _ff_equal (t, w) ) y = -y;
	
	a[0] = -4*y;  a[1] = 8*y*y;
	return 2;
}


// Computes Lpoly coefficients for curves y^2=x^6+a defined over Q, doesn't compute group structure
// note, a[] is the array of L-poly coefficients being computed, the value of the "a" coefficient in the curve equation is sc->f[0]
int smalljac_x6pa_Lpoly (long a[], ff_t f0)
{
	long x, y, k, m, n, ap, a2, p;
	int i;
	ff_t s, t, u;
	ff_t f[7];

	p = _ff_p;
	if ( _ff_zero(f0) ) return -1;			// singular curve
	
	// The Jacobian of y^2=x^6+1 is Q-isogenous to the square of the elliptic curve y^2=x^3+1, which we handle separately
	if ( _ff_one(f0) ) { if ( ! smalljac_x3pa_Lpoly (a, f0) ) return 0;  a[1] = a[0]*a[0]+2*p;  a[0] *=2; return 2; }
	
	
	// If p is not 1 mod 3, the Hasse-witt matrix is 0, so a_1=a_2=0 mod p, and the curve is supersingular
	// by Thm 6.1 of Galbraith "Supersingular curves in cryptography", and the L-poly must be (pT^2+1)^2.
	if ( p%3 != 1 ) { a[0] = 0;  a[1] = 2*p;  return 2; }
	
	// p = 6m+1
	// the Hasse Witt matrix is diagonal, with entries binom(3m,m)f0^{2m} and binom(3m,m)f0^m
	// we compute (3m,m) by applying formula 9.1 (p. 453) of "Binomial coefficients and Jacobi Sums"
	// by Hudson and Williams.  Binom(3m,m) = 2x where x^2+3y^2 = p with x = 1 mod 3
	if ( ! ff_cornacchia (&x, &y, 3) ) { printf ("internal error in smalljac_x6pa_lpoly, %ld is not a prime of the form x^2 + 3y^2 ?!\n", p);  return -2; }
	if ( (x%3) != 1 ) x = -x;	// make sure x is 1 mod 3
	// binom(3m,m) = 2x

	m = (p-1)/6;
	ff_exp_ui (&s, &f0, m);
	_ff_square (t, s);
	_ff_set_i (u, x);
	_ff_x2 (u);					// u = binom(3m,m)
	_ff_mult(s,s,u);
	_ff_mult(t,t,u);
	_ff_add(u,s,t);
	ap = _ff_get_ui(u);				// trace = s+t = (binom(3m,3)f0^{2m} + binom(3m,m)f0^m mod p
	_ff_mult(u,s,t);
	a2 = _ff_get_ui(u);				// a2 = s*t mod p
	if ( ap > p/2 ) ap -= p;

	// We now know the L-poly mod p, we just need to nail down a_2, which we can do by using the fact that Lp(1)=0 mod 3 and using 2-torsion gives us Lp(1) mod 6
	m = 2*p + (ap*ap)/4;								// upper bound on a2
	if ( m > a2 ) m -= (m-a2)%p; else m -= p - ((a2-m)%p);	// set m to the largest possible value of a2, given a2 mod p
	n = (2-2*ap)%6;									// set n to Lp(1)-a2=p^2-p*ap-ap+1=(p^2+1)-(p+1)ap=2-2*ap modulo 6 (since p = 1 mod 6)
	// Use 2-torsion to determine k=Lp(1) mod 6 (if f0=1 we always have 2-torsion)
	if ( _ff_one(f0) ) {
		k=0;
	} else {
		_ff_set_one (f[6]);  _ff_set(f[0],f0);
		for ( i = 1 ; i < 6 ; i++ ) _ff_set_zero (f[i]);
		k = ( ff_poly_count_factors (f, 6) > 2 ? 0 : 3);		// NOTE: 3-3 split does not yield 2-torsion.  1-5 and 2-4 splits can't happen
	}
	
	// We now know Lp(1)=k mod  6 which uniquely determines a_2	
	while ( (m+n-k)%6 ) m -= p;
	a[0] = -ap;
	a[1] = m;
	return 2;
}


/*
	Compute the L-poly coefficients of the curve C_good: y^2=x^6 - 5x^4 - 5x^2 + 1.
	This curve is a twist of y^2=x^5+x, and it has the distinct advantage that its Jacobian is Q-isogenous to the square of an elliptic curve defined over Q (namely, y^2=x^3-5x^2-5x+1).
	This elliptic curve has CM by the order with discriminant -2, and it is thus *very* easy to compute the L-poly of C_good, even easier than for y^2=x^5+x, and this data can in turn
	be used (in smalljac_g23.c) to quickly compute the L-poly of any twist of y^2=x^5+x (most of which cannot be handled by smalljac_x5pax_Lpoly above)
*/
int smalljac_cm2square_Lpoly (long a[])
{
	long x, y, p;
	
	// primes congruent to 5 or 7 mod 8 are supersingular
	p = _ff_p;
	if ( (p&7) > 3 ) { a[0]=0; a[1]=2*p; return 2; }
	if ( ! ff_cornacchia (&x, &y, 2) ) { printf ("internal error in smalljac_cm2square_Lpoly, %ld is not a prime of the form x^2 + 2y^2 ?!\n", p);  return -2; }
	y = ((x-1)/2)  + (((p-1)/2) * ((p+5)/2)) / 4;
	if ( y&1 ) a[0] = 4*x; else a[0] = -4*x;
	a[1] = 4*x*x + 2*p;
	return 2;
}


int smalljac_x7pax_a1 (long a[], ff_t f1)
{
	long x, y, p;
	ff_t t, u, v, w;

	if ( _ff_zero (f1) ) return -1;			// singular curve
	
	p = _ff_p;
	if ( p < 144 ) { printf ("p=%ld must be greater than 144 to use smalljac_x7pax_a1\n", p); return -2; }

	// If p is not 1 mod 4 then the Hasse Witt matrix is 0
	if ( (p&3) != 1 ) { a[0] = 0;  return 1; }

	// We can write the trace of  the Hasse Witt matrix as [x^n](x^6+f1)^n +  [x^(3n)](x^6+f1)^n + [x^(5n)](x^6+f1)^n, where n=(p-1)/2 (note we factored an x out of x^7+f1*x)
	// This simplifies to binom(n,n/6)*(f1^(5n/6)+f1^(n/6)) + binom(n,n/2)*f1^(n/2), where binom(n,n/6) is zero unless 6|n.
	
	// Using Thm. 9.2.2 of Berndt-Evans-Williams "Gauss & Jacobi Sums", compute binom(2m,m) = 2 * (-1)^(m+1) * x, where p=4m+1=x^2+y^2 with x = - kron(2,p) mod 4  (x is called a_4 in BEW)
	if ( ! ff_cornacchia (&x, &y, 1) ) { printf ("internal error in smalljac_x7pax_a1, %ld is not a prime of the form x^2 + y^2 ?!\n", p);  return -2; }
	
	// choose the sign of x so that x = - kron(2,p) mod 4
	if ( !(x&1) ) x = y;
	if ( (p&7)==1 || (p&7) == 7 ) { if ( (x&3) == 1 ) x = -x; } else { if ( (x&3) == 3 ) x = -x; }
	_ff_set_i (t, x);					// t = a_4
	if ( (p&7) == 1 ) ff_negate(t);		// t = (-1)^(m+1)*a_4
	_ff_x2 (t);						// t = 2*(-1)^(m+1)*a_4
	if ( (p%3) == 1 ) {
		ff_exp_ui (&u, &f1, (p-1)/12);
		_ff_square (v, u); _ff_mult (v, v, u);
	} else {
		ff_exp_ui (&v, &f1, (p-1)/4);
	}
	_ff_mult (w, t, v);										// t = 2 * (-1)^(m+1) * a_4,  w = binom(n,n/2)*f1^(n/2)
	if ( (p%3) == 1 ) {
		// Using Thm 9.2.10 of BEW we compute binom(6m,m) = 2 * (-1)^(m+1) * x * rho^2,  where p=12m+1, x is as above, and rho^2 is -1 if a_4 is 0 mod 3 and +1 otherwise (note (-1)^(m+1) is the same as above)
		if ( !(x%3) ) ff_negate(t);
		_ff_mult (v, v, u); _ff_mult (v, v, u); _ff_addto (u, v);		// u = f1^m + f1^(5m)
		_ff_mult (t, t, u);
		_ff_addto (w, t);									// w =binom(n,n/6)*(f1^(5n/6)+f1^(n/6)) + binom(n,n/2)*f1n^(n/2)
	}
	x = _ff_get_ui (w);
	if ( x > (p>>1) ) x = x-p;									// determine integer value of the trace from its value mod p > 2*6*sqrt(p)
	a[0] = -x;												// a_1 = a[0] is minus the trace
	return 1;
}


int smalljac_x8pa_a1 (long a[], ff_t f0)
{
	long x, y, p;
	ff_t t, u, v, w;

	if ( _ff_zero (f0) ) return -1;			// singular curve
	
	p = _ff_p;
	if ( p < 144 ) { printf ("p=%ld must be greater than 144 to use smalljac_x8pa_a1\n", p); return -2; }

	// If p is not 1 mod 3 than the Hasse Witt matrix is 0
	if ( (p&3) != 1 ) { a[0] = 0;  return 1; }

	// We can write the trace of  the Hasse Witt matrix as [x^(2n)](x^8+f0)^n +  [x^(4n)](x^8+f0)^n + [x^(6n)](x^8+f0)^n, where n=(p-1)/2
	// This simplifies to binom(n,n/4)*(f0^(3n/4)+f0^(n/4)) + binom(n,n/2)*f0^(n/2), where binom(n,n/4) is zero unless 4|n.
	
	// Using Thm. 9.2.2 of Berndt-Evans-Williams "Gauss & Jacobi Sums", compute binom(2m,m) = 2 * (-1)^(m+1) * x, where p=4m+1=x^2+y^2 with x = - kron(2,p) mod 4  (x is called a_4 in BEW)
	// If we were really clever, for p=1 mod 8 we would use p=x^2+2y^2 (computed below) to figure out how to write p=u^2+v^2 and avoid making 2 calls to ff_cornacchia.  But we aren't that clever...
	if ( ! ff_cornacchia (&x, &y, 1) ) { printf ("internal error in smalljac_x8pa_a1, %ld is not a prime of the form x^2 + y^2 ?!\n", p);  return -2; }
	
	// make sure x = - kron(2,p) mod 4
	if ( !(x&1) ) x = y;
	if ( (p&7)==1 || (p&7) == 7 ) { if ( (x&3) == 1 ) x = -x; } else { if ( (x&3) == 3 ) x = -x; }
	if ( (p&7) == 1 ) x = -x;
	_ff_set_i (t, x);
	_ff_x2 (t);
	if ( (p&7) == 1 ) {
		ff_exp_ui (&u, &f0, (p-1)/8);
		_ff_square (v, u);
	} else {
		ff_exp_ui (&v, &f0, (p-1)/4);
	}
	_ff_mult (w, t, v);										// t = 2 * (-1)^(m+1) * x,  w = binom(n,n/2)*f0^(n/2)
	if ( (p&7) == 1 ) {
		// Using Thm 9.2.8 of BEW we compute binom(4m,m) = 2 * (-1)^(m+1) * x,  where p=8m+1=x^2+2y^2 with x=3 mod 4
		if ( ! ff_cornacchia (&x, &y, 2) ) { printf ("internal error in smalljac_x8pa_a1, %ld is not a prime of the form x^2 + 2y^2 ?!\n", p);  return -2; }
		if ( (x&3) != 3 ) x = -x;
		if ( (p&0xf) == 1 ) x = -x;
		_ff_set_i (t, x);
		_ff_x2 (t);
		_ff_mult (v, v, u);  _ff_addto (u, v);						// u = f1^m + f1^(3m)
		_ff_mult (t, t, u);
		_ff_addto (w, t);									// w = binom(n,n/4)*(f0^(3n/4)+f0^(n/4)) + binom(n,n/2)*f0^(n/2)
	}
	x = _ff_get_ui (w);
	if ( x > (p>>1) ) x = x-p;									// determine integer value of the trace from its value mod p > 2*6*sqrt(p)
	a[0] = -x;												// a_1 = a[0] is minus the trace
	return 1;
}

#define FK_CURVES		16

struct {
	long C[15];
	long f[9];
	int d;
} FKtwist_table[FK_CURVES] = {
	{{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, }, {0,0,0,0,0,0,0,0,0}, 0 },			// unused
#define FERMAT_OFFSET 1
	// fermat curves k=Q unless ow specified, M=k(i)
	{{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, }, {1,0,1,0,0,0,0,0,0}, 2 },
	{{ 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -12, 0, 0, 2, }, {1,0,1,0,0,0,0,0,0}, 2 },
	{{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, }, {1,0,1,0,0,0,0,0,0}, 2 },
	{{ 1, 0, 0, 0, -6, 0, 0, 0, -32, 0, 0, 36, 0, 32, -6, }, {1,2,2,2,3,6,8,4,1}, 8 },		// k=Q(sqrt(-5))
	{{ -4, 0, 0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, }, {16,0,-4,0,1,0,0,0,0}, 4 },
	{{ 1, 0, 0, 0, 14, 0, 0, 0, 192, 0, 0, -84, 0, -192, 14, }, {1,0,0,0,1,0,0,0,0}, 4 }, 
	{{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, }, {16,0,-4,0,1,0,0,0,0}, 4 },
	{{ 1, 0, 0, 0, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, }, {16,0,-4,0,1,0,0,0,0}, 4 },
	{{ 1, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 36, 0, 0, 2, }, {16,0,-4,0,1,0,0,0,0}, 4 },
	{{ 1, 0, 0, 0, -72, 0, 0, 0, 0, 0, 0, -144, 0, 0, -8, }, {16,0,-4,0,1,0,0,0,0}, 4 },
	{{ 7546, -9604, 4998, -980, 98, 3332, -2940, 1176, -84, 588, -252, 84, 56, 0, 3, }, {1,0,6,0,5,0,1,0,0}, 6},
	{{ 608, -1344, -768, 160, 24, -512, 480, 288, -16, 144, -48, -24, -16, 0, 1, }, {2,6,9,-2,0,0,1,0,0}, 6},
	{{ 1, 0, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, }, {441,0,372,0,238,0,-28,0,1}, 8},
#define KLEIN_OFFSET 14
	// klein curves
	{{ 1, 6, -3, 0, 1, 0, 3, 3, 6, -3, 3, -3, 6, 0, 1, }, {7,0,1,0,0,0,0,0,0}, 2 },
	{{ -1, 10, -27, 10, -1, 105, -63, 0, 21, -27, 135, -27, 0, 0, 243, }, {5,-3,0,1,0,0,0,0,0}, 3 },
};

int smalljac_FKtwist_lookup (mpz_t C[15])
{
	register int i, j;
	
	for ( i = 1 ; i < sizeof(FKtwist_table)/sizeof(FKtwist_table[0]) ; i++ ) {
		for ( j = 0 ; j < 15 ; j++ ) if ( mpz_cmp_si (C[j], FKtwist_table[i].C[j]) ) break;
		if ( j == 15 ) return i;
	}
	return 0;
}

int smalljac_g1_CM7t_Lpoly (long a[])
{
	long u, v;
	long d;
	
	d = _ff_p % 7;
	if ( d == 3 || d == 5 || d == 6 ) { a[0] = 0; return 1; }
	if ( ! ff_cornacchia (&u, &v, 7) ) { printf ("%ld is not a prime of the form x^2 + 7y^2 but p=%ld mod 7\n", _ff_p, _ff_p);  return -2; }
	d = u % 7;
	if ( d == 3 || d == 5 || d == 6 ) u = -u;
	a[0] = -2*u;
	return 1;
}

int smalljac_FKtwist_Lpoly (long a[], int id)
{
	ff_t f[9], g[9];
	int cnts[9];
	long a0,a1, a2, a3,p;
	ff_t w;
	int i,n,d,dg,s,t;
	
	if ( id <= 0 || id > FK_CURVES ) { printf ("Unknown curve id %d in smalljac_FKtwist_Lpoly\n", id); return -2; }
	// compute a0, s, t
	p = _ff_p;
	if ( p == 2 ) return -1;
	if ( id < KLEIN_OFFSET ) {
		if ( (p&2) ) t = 2; else t = 1;
		if ( t == 1 ) {
			_ff_set_one(w);
			n = smalljac_x3pax_Lpoly (&a0, w);
			if ( n < 1 ) return n;
		} else {
			a0 = 0;
		}
	} else {
		if ( p == 7 ) return -1;
		a0=(p%7);
		if ( a0==1 || a0==2 || a0==4 ) t = 1; else t = 2;
		if ( t == 1 ) {
			n = smalljac_g1_CM7t_Lpoly (&a0);
			if ( n < 1 ) return n;
		} else {
			a0 = 0;
		}
	}
	d = FKtwist_table[id].d;
	assert ( d > 1 );
	if ( d == 2 ) {
		s = t;
	} else {
		// This code does not check for ramification!  This would be easy to add...
		for ( i = 0 ; i <= d ; i++ ) _ff_set_i (f[i], FKtwist_table[id].f[i]);
		ff_poly_derivative (g, &dg, f, d);
		if ( ff_poly_gcd (g, &dg, f, d, g, dg) > 0 ) { printf ("Skipping ramified prime %ld\n", p); return -1; }		
		ff_poly_factorization_pattern (cnts, f, d);
		for ( i = 0 ; i <= d ; i++ ) if ( cnts[i] ) break;
		assert ( i <= d );
		s = i;
	}

	if ( t == 1 ) {
		switch (s) {
		case 1: a1 = 3*a0;  a2 = 3*(a0*a0+p);  a3 = a0*(a0*a0+6*p);  break;
		case 2: a1 = -a0;  a2 = -a0*a0+3*p;  a3 = a0*a0*a0-2*p*a0;  break;
		case 3: a1 = a2 = 0;  a3 = a0*(a0*a0-3*p); break;
		case 4: if ( id >= KLEIN_OFFSET ) { a1 = 0; a2 = a0*a0-p; a3=a0*a0*a0-2*a0*p; break; }
		default: printf ("Unhandled degree pair <%d,%d> for special Fermat curve with id %d at p=%ld\n", s, t, id, p); abort();
		}
	} else {
		a1 = a3 = 0;
		switch (s) {
		case 2: a2 = 3*p;  break;
		case 4: a2 = -p;  break;
		case 6: a2 = 0;  break;
		case 8: a2 = p;  break;
		default: printf ("Unhandled degree pair <%d,%d> for special Fermat curve with id %d at p=%ld\n", s, t, id, p); abort();
		}
	}

//printf ("p=%ld,s=%d,t=%d) a0=%ld, a1=%ld, a2=%ld, a3=%ld\n", p, s, t, a0, a1, a2, a3);
	a[0] = a1;  a[1] = a2;  a[2] = a3;
	return 3;
}
