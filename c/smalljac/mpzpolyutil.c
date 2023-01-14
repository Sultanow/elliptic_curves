#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>
#include "polyparse.h"
#include "mpzpolyutil.h"
#include "cstd.h"

/*
    Copyright 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

int mpz_rational_function_init (mpz_rational_function_t f, char *expr)
{
	register char *s;
	register int d,i;
	
	for ( s = expr ; isspace(*s) ; s++ );
	if ( *s != '[' && *s != '(' ) {
		if ( strchr (s, '/') ) { printf ("Did not find expected leading parenthesis in rational expression -- please parenthesize numerator\n");  return 0; }		
	}
	d = poly_parse_max_degree (s);
	if ( d < -1 ) return 0;
	mpq_init (f->y); f->den = 0;
	if ( d == -1 ) { f->num = 0; f->d_num = -1; return 1; }
	f->num = mem_alloc ((d+1)*sizeof(*f->num));
	for ( i = 0 ; i <= d ; i++ ) mpz_init (f->num[i]);
	f->d_num = mpz_poly_parse (f->num, d, expr);
	if ( f->d_num < -1 ) { mpz_rational_function_clear(f); return 0; }
	for ( ; *s && *s != ']' && *s != ')' ; s++ );
	if ( ! *s ) return 1;
	for ( s++ ; isspace(*s) ; s++ );
	if ( *s != '/' ) return 1;
	for ( s++ ; isspace(*s); s++ );
	if ( *s != '[' && *s != '(' ) { printf ("Did not find expected parenthesis in rational expression -- please parenthesize denominmator\n"); mpz_rational_function_clear(f); return 0; }
	d = poly_parse_max_degree (s);
	if ( d < 0 )  { mpz_rational_function_clear(f); return 0; }
	f->den = mem_alloc ((d+1)*sizeof(*f->den));
	for ( i = 0 ; i <= d ; i++ ) mpz_init (f->den[i]);
	f->d_den = mpz_poly_parse (f->den, d, s);
	if ( f->d_den < 0 ) { mpz_rational_function_clear(f); return 0; }
	return 1;
}

void mpz_rational_function_clear (mpz_rational_function_t f)
{
	int i;
	
	mpq_clear (f->y);
	for ( i = 0 ; i <= f->d_num ; i++ ) mpz_clear (f->num[i]);
	if ( f->num ) mem_free (f->num);
	if ( f->den ) {
		for ( i = 0 ; i <= f->d_den ; i++ ) mpz_clear (f->den[i]);
		mem_free (f->den);
	}
}

void mpq_ws_make_integral (mpq_t W[5])	// does not necessarily minimize the coefficients
{
	mpz_t d, x;
	mpq_t z;
	
	mpz_init (d);  mpz_init (x);  mpq_init (z);
	mpz_lcm (d,mpq_denref(W[0]),mpq_denref(W[1])); mpz_lcm (d,d,mpq_denref(W[2])); mpz_lcm (d,d,mpq_denref(W[3])); mpz_lcm (d,d,mpq_denref(W[4]));
	if ( mpz_cmp_ui(d,1) == 0 ) return;
	mpq_set_z (z,d);  mpq_mul (W[0], W[0], z);
	mpz_mul(x,d,d); mpq_set_z (z,x);  mpq_mul (W[1], W[1], z);
	mpz_mul(x,x,d); mpq_set_z (z,x);  mpq_mul (W[2], W[2], z);
	mpz_mul(x,x,d); mpq_set_z (z,x);  mpq_mul (W[3], W[3], z);
	mpz_mul(x,x,d); mpz_mul(x,x,d); mpq_set_z (z,x);  mpq_mul (W[4], W[4], z);
}


void mpz_short_ws (mpz_t A, mpz_t B, mpz_t a1, mpz_t a2, mpz_t a3, mpz_t a4, mpz_t a6)
{
	mpz_t b2, b4, b6, t1, t2;
	
	if ( ! mpz_sgn(a1) && ! mpz_sgn(a2) && ! mpz_sgn(a3) ) { mpz_set (A,a4); mpz_set (B,a6); return; }
	mpz_init(b2);  mpz_init(b4);  mpz_init(b6);  mpz_init(t1); mpz_init(t2);
	mpz_mul(t1,a1,a1);  mpz_mul_2exp (b2,a2,2);  mpz_add(b2,t1,b2);			// b2 = a1^2 + 4a2, t1 = a1^2 (fix typo in Silverman p. 42!!!!)
	mpz_add(b4,a4,a4);  mpz_mul (t2,a1,a3); mpz_add(b4,b4,t2);				// b4 = 2a4 + a1a3, t2 = a1a3
	mpz_mul(t1,a3,a3);  mpz_mul_2exp (b6,a6,2); mpz_add(b6,t1,b6);			// b6 = a3^2 + 4a6
//	mpz_mul(b8,t1,a6); mpz_mul(t1,a2,b6); mpz_add(b8,b8,t1);					// b8 = a1^2a6 + a2b6
//	mpz_add(t2,t2,a4); mpz_mul(t1,t2,a4); mpz_sub(b8,b8,t1);					// b8 = a1^2a6 + a2b6 - a4(a1a3+a4)
	mpz_mul(t2,b2,b2); mpz_mul_ui(t1,b4,24); mpz_sub(A,t1,t2);				// A = -c4 = 24b4-b2^2, t2 = b2^2
	mpz_mul_ui(t1,b4,36);  mpz_sub(t2,t2,t1); mpz_mul(t1,t2,b2);				// t1 = b2^3 - 36b2b4
	mpz_mul_ui(t2,b6,216); mpz_add(B,t1,t2);								// B = -c6 = b2^3 - 36b2b4 + 216b6
	mpz_mul_ui(A,A,27);  mpz_mul_ui(B,B,54);								// A = -27c4, B=-54c6
	mpz_clear (b2); mpz_clear (b4); mpz_clear (b6); mpz_clear (t1); mpz_clear (t2);
}


static long cm_jinv[13] = { 0L, 1728L, -3375L, 8000L, -32768L, 54000L, -884736L, 287496L, -12288000L, 16581375L, -884736000L, -147197952000L, -262537412640768000L };

int mpz_jinv_has_cm (mpz_t J)
{
	register int i;
	
	if ( mpz_sizeinbase (J, 2) > 61 ) return 0;
	for ( i = 0 ; i < 13 ; i++ ) if ( mpz_cmp_si(J, cm_jinv[i]) == 0 ) return 1;
	return 0;
}


void _mpq_coeff_clear (void *f, int i) { mpq_set_ui (((mpq_t*)f)[i], 0, 1); }
int _mpq_coeff_addto (void *f, int i, mpq_t c, void *unused) { mpq_add (((mpq_t*)f)[i], ((mpq_t*)f)[i], c); return 1; }
int _mpq_coeff_iszero (void *f, int i) { return ! mpq_sgn (((mpq_t*)f)[i]); }

int mpq_poly_parse (mpq_t f[], int maxd, char *expr) { return poly_parse (f, maxd, expr, _mpq_coeff_clear, _mpq_coeff_addto, _mpq_coeff_iszero, 0); }
int mpq_poly_parse_plane_quartic (mpq_t f[], char *expr) { return poly_parse_plane_quartic (f, expr, _mpq_coeff_clear, _mpq_coeff_addto, _mpq_coeff_iszero, 0); }

void _mpz_coeff_clear (void *f, int i) { mpz_set_ui (((mpz_t*)f)[i], 0); }
int _mpz_coeff_addto (void *f, int i, mpq_t c, void *unused) { if(  mpz_cmp_ui (mpq_denref(c),1) ) return 0;  mpz_add (((mpz_t*)f)[i], ((mpz_t*)f)[i], mpq_numref(c)); return 1; }
int _mpz_coeff_iszero (void *f, int i) { return ! mpz_sgn (((mpz_t*)f)[i]); }

int mpz_poly_parse (mpz_t f[], int maxd, char *expr) { return poly_parse (f, maxd, expr, _mpz_coeff_clear, _mpz_coeff_addto, _mpz_coeff_iszero, 0); }
int mpz_poly_parse_plane_quartic (mpz_t f[], char *expr) { return poly_parse_plane_quartic (f, expr, _mpz_coeff_clear, _mpz_coeff_addto, _mpz_coeff_iszero, 0); }

// discriminant and resultant code (such as it is) starts here

void mpz_print_matrix (mpz_t *M, int n)
{
	int i,j;
	
	for ( i = 0 ; i < n ; i++ ) { for ( j = 0 ; j < n ; j++ ) gmp_printf ("%4Zd  ", M[i*n+j]); puts (""); }
	puts ("");			
}

// compute the determinant of an integer n*n matrix using simple Bareiss algorithm
// no attempt is made to be clever about choosing pivots or to take advantage of sparseness
// M specifiies an array whose i*n+j entry holds the (i,j) entery of the matrix (indexed from zero).
// M is destroyed in the process
void mpz_determinant (mpz_t D, mpz_t *M, int n, mpz_t w[2])
{
	register int i, j, k, s;

	// the general code below works fine in these small cases, but its a little faster to hard wire them
	switch (n) {
	case 1: mpz_set(D,M[0]); return;
	case 2: mpz_mul(w[0],M[0],M[3]);  mpz_mul(w[1],M[1],M[2]);  mpz_sub(D,w[0],w[1]); return;
	case 3: mpz_mul(w[0],M[4],M[8]);  mpz_mul(w[1],M[5],M[7]);  mpz_sub(w[0],w[0],w[1]);  mpz_mul(M[0],M[0],w[0]);
		     mpz_mul(w[0],M[1],M[8]);  mpz_mul(w[1],M[2],M[7]);  mpz_sub(w[0],w[0],w[1]);  mpz_mul(M[3],M[3],w[0]);
		     mpz_mul(w[0],M[1],M[5]);  mpz_mul(w[1],M[2],M[4]);  mpz_sub(w[0],w[0],w[1]);  mpz_mul(M[6],M[6],w[0]);
		     mpz_sub(D,M[0],M[3]);  mpz_add(D,D,M[6]);  return;
	}
	s = 1;
	for ( k = 0 ; k < n-1 ; k++ ) {
		if ( ! mpz_sgn(M[k*n+k]) ) {			// if next diagonal entry is zero, look for a row-swap
			for ( i = k+1 ; i < n ; i++ ) if ( mpz_sgn(M[i*n+k]) ) break;
			if ( i == n ) { mpz_set_ui(D,0); return; }
			for ( j = k ; j < n ; j++ ) mpz_swap (M[i*n+j], M[k*n+j]);
			s *= -1;
		}
		for ( i = k+1 ; i < n ; i++ ) {
			for ( j = k+1 ; j < n ; j++ ) {
				mpz_mul(w[0],M[i*n+j],M[k*n+k]);
				mpz_mul(w[1],M[i*n+k],M[k*n+j]);
				mpz_sub (w[0],w[0],w[1]);
				// note that divexact doesn't seem to like zero numerator (which it really should handle...)
				if ( k >= 1 && mpz_sgn(w[0]) ) mpz_divexact(M[i*n+j],w[0],M[(k-1)*n+(k-1)]); else mpz_set (M[i*n+j],w[0]);
			}
		}
	}
	if ( s > 0 ) mpz_set (D, M[n*n-1]); else mpz_neg (D, M[n*n-1]);
}


// general discriminant code -- w must point to d*d+2 initialized mpz_t's, this code is thread safe
void mpz_poly_discriminant (mpz_t D, mpz_t f[], int d, mpz_t *w)
{
	mpz_t *M;
	register int i, k;
	
	assert ( mpz_sgn(f[d]) );
	if ( d < 1 ) { mpz_set_ui(D,0); return; }
	if ( d == 1 ) { mpz_set(D,f[1]); return; }
	if ( d == 2 ) { mpz_mul(D,f[1],f[1]);  mpz_mul(w[0],f[0],f[2]); mpz_mul_2exp(w[0],w[0],2); mpz_sub(D,D,w[0]); return; }

	// row reduce the Sylvester matrix to a d x d matrix and then compute the determinant
	M = w+2;
		
	// reverse the order of the rows so that smaller coefficients are on top (also removes the need to adjust the sign at then end)
	for ( i = 0 ; i < d ; i++ ) mpz_mul_ui(M[i],f[d-i],d-i);			// first row is the derivative of f, with terms stored from high to low degree
	for ( i = 0 ; i < d-1 ; i++ )
		{ mpz_mul_ui(w[0],f[d-i-1],d);  mpz_sub(M[d+i],M[i+1],w[0]); }	// second row is x*f'-d*f
	mpz_mul_ui(w[0],f[0],d);  mpz_neg(M[d+i],w[0]);
	for ( k = 2 ; k < d ; k++ ) {
		for ( i = 0 ; i < d-1 ; i++ ) {
			mpz_mul(w[0],f[d],M[(k-1)*d+i+1]);  mpz_mul(w[1],M[(k-1)*d],f[d-i-1]);
			mpz_sub (M[k*d+i],w[0],w[1]);
		}
		mpz_mul(w[0],M[(k-1)*d],f[0]);  mpz_neg (M[k*d+i],w[0]);
	}
	mpz_determinant (D, M, d, w);
	mpz_pow_ui (w[0], f[d], ((d-2)*(d-3))/2);
	mpz_divexact(D,D,w[0]);
}


void mpz_fast_poly_discriminant (mpz_t D, mpz_t f[], int d)
{
	static int init;
	static mpz_t w1, w2, w3, w4, f02, f12, f22, f32, f33;
	static mpz_t *w;
	static int maxd;
	int i;

	// hack to free up statically allocated resources
	if ( ! D ) {
		if ( ! init ) return;
		mpz_clear (w1);  mpz_clear (w2);  mpz_clear(w3); mpz_clear (w4); mpz_clear (f02); mpz_clear (f12); mpz_clear (f22); mpz_clear (f32); mpz_clear (f33);
		if ( maxd ) for ( i = 0 ; i < maxd*maxd+2 ; i++ ) mpz_clear (w[i]);
		free (w);
		init = maxd = 0;
		return;
	}
	
	// statically allocate -- this function needs to be as fast as we can make it, it may get called a lot (of course this is NOT thread safe)
	if ( ! init ) { mpz_init (w1);  mpz_init (w2);  mpz_init(w3); mpz_init (w4); mpz_init (f02); mpz_init (f12); mpz_init (f22); mpz_init (f32); mpz_init (f33); init = 1; }
	
	assert ( mpz_sgn(f[d]) );

	// if poly is monic and depressed
	if ( mpz_cmp_ui (f[d],1) == 0 && ! mpz_sgn(f[d-1]) ) {
		for ( i = d-1 ; i > 1 ; i-- ) if ( mpz_sgn(f[i]) ) break;
		if ( i == 1 ) {
			// optimized code for x^d+a*x+b case  D = s1*(d-1)^(d-1)*f1^d + s0*d^d*f0^(d-1), where s1 = 1 if (d-1)=0,1 mod 4 (s1=-1 ow) and s0 = 1 if d=0,1 mod 4 (s0=-1 ow)
			if ( mpz_sgn(f[1]) ) {
				mpz_ui_pow_ui (w1, (d-1), (d-1));
				mpz_pow_ui (D, f[1], d);
				mpz_mul (D, D, w1);
				if ( ((d-1)&3) >= 2 ) mpz_neg (D, D);
			} else {
				mpz_set_ui(D,0);
			}
			if ( mpz_sgn(f[0]) ) {
				mpz_ui_pow_ui (w1, d, d);
				mpz_pow_ui (w2, f[0], d-1);
				mpz_mul (w2, w2, w1);
				if ( (d&3) < 2) mpz_add (D, D, w2); else mpz_sub (D, D, w2);
			}
			return;
		} else if ( (d == 5)  ) {
			// optimized code for x^5+f3*x^3+f2*x^2+f1*x+f0
			mpz_mul(f02,f[0],f[0]);  mpz_mul(f12,f[1],f[1]);  mpz_mul(f22,f[2],f[2]);  mpz_mul(f32,f[3],f[3]);  mpz_mul(f33,f32,f[3]);
			// D := f0 * (f0*(f0*(3125*f0 - 3750*f2*f3) + f1*(2000*f1*f3 + 2250*f22 - 900*f33) + 825*f22*f32 + 108*f32*f33) + f2*(f1*(560*f1*f32 -1600*f12 - 630*f22*f3 - 72*f32^2) + f22*(108*f22 + 16*f33)))
			//        + f12 * (f1*(256*f12 - 128*f1*f32 + f3*(144*f22 + 16*f33)) - f22*(27*f22 + 4*f33));
			if ( mpz_sgn(f[0]) ) {
				mpz_mul_ui(w1,f[0],3125);  mpz_mul(w2,f[2],f[3]); mpz_mul_ui(w2,w2,3750); mpz_sub(w1,w1,w2); mpz_mul(w1,w1,f[0]);
				mpz_mul(w2,f[1],f[3]); mpz_mul_ui(w2,w2,2000); mpz_mul_ui(w3,f22,2250); mpz_add(w2,w2,w3); mpz_mul_ui(w3,f33,900); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[1]);mpz_add(w1,w1,w2);
				mpz_mul(w2,f22,f32); mpz_mul_ui(w2,w2,825); mpz_add(w1,w1,w2); mpz_mul(w2,f32,f33); mpz_mul_ui(w2,w2,108); mpz_add(w1,w1,w2); mpz_mul(w1,w1,f02);
				mpz_mul(w2,f[1],f32); mpz_mul_ui(w2,w2,560); mpz_mul_ui(w3,f12,1600); mpz_sub(w2,w2,w3); mpz_mul(w3,f22,f[3]); mpz_mul_ui(w3,w3,630); mpz_sub(w2,w2,w3);
				mpz_mul(w3,f32,f32); mpz_mul_ui(w3,w3,72); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[2]); mpz_mul(w2,w2,f[1]);
				mpz_mul_ui(w3,f22,108); mpz_mul_2exp(w4,f33,4); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f22); mpz_mul(w3,w3,f[2]); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f[0]); mpz_add(w1,w1,w2);
			} else {
				mpz_set_ui(w1,0);
			}
			if ( mpz_sgn(f[1]) ) {
				mpz_mul_2exp(w2,f12,8); mpz_mul(w3,f[1],f32); mpz_mul_2exp(w3,w3,7); mpz_sub(w2,w2,w3); mpz_mul_ui(w3,f22,144); mpz_mul_2exp(w4,f33,4); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f[3]); mpz_add(w2,w2,w3);
				mpz_mul(w2,w2,f[1]); mpz_mul_ui(w3,f22,27); mpz_mul_2exp(w4,f33,2); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f22); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f12);
			} else {
				mpz_set_ui(w2,0);
			}
			mpz_add (D,w1,w2);
			return;
		}
	}
	if ( d > maxd ) {
		w = realloc (w, (d*d+2)*sizeof(mpz_t));
		for ( i = maxd ; i < d*d+2 ; i++ ) mpz_init (w[i]);
		maxd = d;
	}
	// general case
	 mpz_poly_discriminant(D,f,d,w);
}

// thread safe, we don't worry particularly about efficiency here
void mpz_poly_resultant (mpz_t R, mpz_t f[], int d_f, mpz_t g[], int d_g)
{
	mpz_t *M, w[2];
	register int i, j, n;
	
	// just create the Sylvester matrix and compute its determinant, don't try to be clever (currently this code is used only for testing)
	mpz_init (w[0]); mpz_init (w[1]);
	n = d_f+d_g;
	M = mem_alloc (n*n*sizeof(mpz_t));
	for ( i = 0 ; i < n*n ; i++ ) mpz_init(M[i]);
	for ( i = 0 ; i < d_g ; i++ ) for ( j = i ; j <= i+d_f ; j++ ) mpz_set (M[i*n+j], f[i+d_f-j]);
	for ( ; i < n ; i++ ) for ( j = i-d_g ; j <= i ; j++ ) mpz_set (M[i*n+j], g[i-j]);
	mpz_determinant (R, M, n, w);
	mpz_clear (w[0]); mpz_clear (w[1]);
	for ( i = 0 ; i < n*n ; i++ ) mpz_clear(M[i]);
	mem_free (M);
}

// parses a hyperelliptic curve in the form [f(x),h(x)], also accepts f(x),h(x), or [f(x)] or f(x) (missing h(x) is 0), where f and h are both integra;l
// defining the hyperelliptic curve y^2 + h(x)y = f(x)
// maxd bounds the degree of f, (maxd+1)/2 bounds the degree of h
// returns 0/1 for failure/success
int mpz_hyperelliptic_curve_parse (mpz_t f[], int *pdf, mpz_t h[], int *pdh, int maxd, char *str)
{
	char *s;
	int d;
	
	*pdf = mpz_poly_parse (f, maxd, str);
	if ( *pdf < 0 ) return 0;	// note that f must by nonzero (o.w. curve would be singular)
	s = strchr (str, ',');
	if ( ! s ) { *pdh = -1; return ((*pdf)-1)/2; }
	*pdh = mpz_poly_parse (h, (maxd+1)/2, s+1);
	if ( *pdh < -1 ) return 0;
	d = 2*(*pdh) > (*pdf) ? 2*(*pdh) : *pdf;
	return (d-1)/2;
}

// computes the discriminant of the hyperelliptic curve y^2 + h(x)y = f(x)
// replacing y by (y-h/2) and multiply through by 4 gives the isomorphic curve y^2 = k(x) = 4*f(x)+h(x)^2 with genus = floor((deg(k)-1)/2)
// if k(x) has odd degree then the discriminant of the curve is lc(k)^2*disc(k) / 2^(4g+4)  (note that lc(k) = 4*lc(f), so this is equivalent to lc(f)^2*disc(k) / 2^(4g))
// if k(x) has even degree then the discriminant of the curve is disc(k) / 2^(4g+4)
// This is taken from the Magma code and agrees with Lockhart's paper in the odd degree case
void mpz_hyperelliptic_curve_discriminant (mpz_t D, mpz_t f[], int df, mpz_t h[], int dh)
{
	mpz_t *k;
	mpz_t *w;
	int dk, i, j, g;
	
	assert (df >= 0);
	dk = df;
	if ( 2*dh > dk ) dk = 2*dh;
	k = malloc ((dk+1)*sizeof(*k));
	for ( i = 0 ; i <= dk ; i++ ) mpz_init (k[i]);
	for ( i = 0 ; i <= dh ; i++ ) for ( j = 0 ; j <= dh ; j++ ) mpz_addmul (k[i+j], h[i], h[j]);
	for ( i = 0 ; i <= df ; i++ ) mpz_addmul_ui (k[i], f[i], 4);
	g = (dk-1)/2;
	w = malloc ((dk*dk+2)*sizeof(*w));
	for ( i = 0 ; i < dk*dk+2 ; i++ ) mpz_init (w[i]);
	mpz_poly_discriminant (D, k, dk, w);
	for ( i = 0 ; i < dk*dk+2 ; i++ ) mpz_clear(w[i]);
	if ( (dk&1) ) { mpz_mul (D, D, k[dk]);  mpz_mul (D, D, k[dk]); }
	mpz_div_2exp (D, D, 4*g+4);
	for ( i = 0 ; i <= dk ; i++ ) mpz_clear (k[i]);
	free (k);
}

// sets g to 4*f+h^2, remove factor of 4 if possible.  g and f may overlap
int mpz_hyperelliptic_curve_normalize (mpz_t g[], mpz_t f[], int df, mpz_t h[], int dh)
{
	register int i, j;
	
	if ( dh < 0 ) {
		if ( g != f ) for ( i = 0 ; i <= df ; i++ ) mpz_set (g[i], f[i]);
		return df;
	}
	for ( i = 0 ; i <= df ; i++ ) mpz_mul_2exp (g[i], f[i], 2);
	while ( i <= 2*dh ) mpz_set_ui(g[i++], 0);
	for ( i = 0 ; i <= dh ; i++ ) for ( j = 0 ; j <= dh ; j++ ) mpz_addmul (g[i+j],h[i],h[j]);
	if ( 2*dh > df ) df = 2*dh;
	for ( i = 0 ; i <= df ; i++ ) if ( ! mpz_divisible_2exp_p (g[i], 2) ) break;
	if ( i > df ) for ( i = 0 ; i <= df ; i++ ) mpz_div_2exp (g[i], g[i], 2);
	return df;
}

// dynamically allocate buffer, it could be big
void mpz_poly_print (mpz_t f[], int df)
{
	char *buf;
	register int i, n;
	
	for ( n = 4, i = 0 ; i <= df ; i++ ) n += 16 + mpz_sizeinbase (f[i], 10);	// 16 bytes accounts for sign, space, x^, and exponent < 2^31
	assert ( n > 0 );
	buf = malloc (n);
	i = mpz_poly_sprint (buf, f, df);
	assert ( i > 0 && i < n );
	puts (buf);
	free (buf);
}


int mpz_poly_sprint (char *s, mpz_t f[], int df)
{
	register char *t;
	register int i;
	
	while ( df >= 0 && ! f[df] ) df--;
	if ( df < 0 ) { strcpy (s, "[0]");  return strlen(s); }
	t = s;
	if ( df >= 2 ) {
		if ( mpz_cmp_ui (f[df],1) == 0 ) { t += gmp_sprintf (t, "[x^%d", df); } else if ( mpz_cmp_si (f[df],-1) == 0 ) { t += gmp_sprintf (t, "[-x^%d", df); } else { t += gmp_sprintf (t, "[%Zd*x^%d", f[df], df); }
	} else if ( df == 1 ) {
		if ( mpz_cmp_ui (f[df],1) == 0 ) { t += gmp_sprintf (t, "[x"); } else if ( mpz_cmp_si(f[df],-1) == 0 ) { t += gmp_sprintf (t, "[-x"); } else { t += gmp_sprintf (t, "[%Zd*x", f[df]); }
	} else {
		t += gmp_sprintf (t, "[%Zd", f[df]);
	}
	for ( i = df-1 ; i >= 0 ; i-- ) {
		if ( mpz_cmp_ui (f[i],1) > 0 ) {
			t += gmp_sprintf (t, " + %Zd", f[i]); if ( i ) *t++ = '*';
		} else if ( mpz_cmp_ui(f[i],1) == 0 ) {
			t += gmp_sprintf (t, " + ");
		} else if ( mpz_cmp_si(f[i],-1) == 0 ) {
			t += gmp_sprintf (t, " - ");
		} else if ( mpz_cmp_si(f[i],-1) < 0 ) {
			mpz_neg (f[i], f[i]);
			t += gmp_sprintf (t, " - %Zd", f[i]); if ( i ) *t++ = '*';
			mpz_neg (f[i], f[i]);
		}
		if ( mpz_sgn(f[i]) ) {
			if ( i >= 2 ) {
				t += gmp_sprintf (t, "x^%d", i);
			} else if ( i == 1 ) {
				t += gmp_sprintf (t, "x");
			} else {
				if ( mpz_cmpabs_ui (f[i],1)==0 ) *t++ = '1';
			}
		}
	}
	*t++ = ']';
	*t= '\0';
	return t-s;
}
