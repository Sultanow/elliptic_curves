#ifndef _JAC_INCLUDE_
#define _JAC_INCLUDE_

/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff_poly.h"
#include "hecurve.h"
#include "hcpoly.h"
#include "smalljac.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
	Basic group operations for Jacobians of hyperelliptic curves are defined here,
	plus some generic fastorder computation functions (computations of element
	orders from a known exponent).
*/	

#define JAC_MAX_CURVES			100
#define JAC_MAX_GAPS				64
#define JAC_BUFSIZE				4096
#define JAC_MAX_FACTORS			100
#define JAC_U_DEGREE				HECURVE_GENUS
#define JAC_V_DEGREE				HECURVE_GENUS-1
#define JAC_CONFIDENCE_LEVEL		20			// larger value used internally in jac_sylow - note that smalljac results are (except where noted) provably correct regardless of this setting
#define JAC_MAX_GENERATORS		9			// this should be at least 2g+1 (not 2g, we need room for one extra)
#define JAC_CYCLIC_TESTS			4			// number of times to attempt to prove p-Sylow subgroup is cyclic before computing it

unsigned long jac_gops;

typedef struct {
	ff_t u[JAC_U_DEGREE+1];
	ff_t v[JAC_V_DEGREE+1];
} jac_t;

struct _jac_vector_struct {
	jac_t a[JAC_MAX_CURVES];
};
typedef struct _jac_vector_struct jac_vec_t;

#define _jac_set(o,a)			_hecurve_set((o).u,(o).v,(a).u,(a).v)
#define __jac_mult(o,a,b,c,ctx)		hecurve_compose ((o).u, (o).v, (a).u, (a).v, (b).u, (b).v, (c).f, ctx)
#define __jac_square(o,a,c,ctx)	hecurve_square ((o).u, (o).v, (a).u, (a).v, (c).f, ctx)
#define _jac_mult(o,a,b,c)		{ __jac_mult (o,a,b,c,0);  jac_gops++; }
#define _jac_square(o,a,c)		{ __jac_square (o,a,c,0);  jac_gops++; }
#define _jac_invert(o)			hecurve_invert((o).u,(o).v)
#define _jac_random(o,c)			hecurve_random((o).u,(o).v,(c).f)

#define _jac_cmp(a,b)			hecurve_cmp ((a).u,(a).v,(b).u,(b).v)  	// 1 if equal, -1 if inverses, 0 o.w.
#define _jac_set_identity(a)		_hecurve_set_identity ((a).u,(a).v)
#define _jac_is_identity(a)		_hecurve_is_identity ((a).u,(a).v)
#define _jac_2tor(a)				_hecurve_2tor((a).u,(a).v)
#define _jac_bits(a,c)			hecurve_bits ((a).u,(a).v,(c).f)

#define _jac_print(a)			hecurve_print((a).u, (a).v)
#define _jac_sprint(s,a)			hecurve_sprint(s,(a).u,(a).v)

static inline void jac_parallel_set (jac_t a[], jac_t b[], int n) { register int i;  for ( i = 0 ; i < n ; i++ ) { _jac_set(a[i],b[i]); }}

void jac_parallel_mult (jac_t o[], jac_t a[], jac_t b[], hc_poly c[], int n);
void jac_parallel_mult_c (jac_t o[], jac_t a[], jac_t b[], hc_poly c[1], int n);		// multiple a's and b's but same c
void jac_parallel_mult_1 (jac_t o[], jac_t a[], jac_t b[1], hc_poly c[1], int n);
void jac_parallel_square (jac_t o[], jac_t a[], hc_poly c[], int n);

void jac_square_mult (jac_t a[1], jac_t b[1], hc_poly c[1]);
void jac_mult2 (jac_t o1[1], jac_t a1[1], jac_t b1[1], jac_t o2[1], jac_t a2[1], jac_t b2[1], hc_poly c[1]);	// o1 can't overlap a2 or b2

void jac_exp_ui (jac_t o[1], jac_t a[1], unsigned long e, hc_poly c[1]);
static inline void jac_exp_si (jac_t o[1], jac_t a[1], long e, hc_poly c[1])
	{ if ( e >= 0 ) jac_exp_ui (o, a, e, c); else { jac_exp_ui (o, a, -e, c); _jac_invert(o[0]); } }
void jac_exp2_ui (jac_t o[2], jac_t a[1], unsigned long e1, unsigned long e2, hc_poly c[1]);
void jac_exp_ui_powers (jac_t o[1], jac_t a[], unsigned long e, hc_poly c[1]);			// o cannot overlap a, computes g^e given a[i]=g^(2^i) for 2^i <= e
static inline void jac_exp_si_powers  (jac_t o[1], jac_t a[], long e, hc_poly c[1])
	{ if ( e >= 0 ) jac_exp_ui_powers (o, a, e, c); else { jac_exp_ui_powers (o, a, -e, c); _jac_invert(o[0]); } }
void jac_exp_mpz (jac_t o[1], jac_t a[1], mpz_t e, hc_poly c[1]);
int jac_verify_group_exponent (mpz_t e, hc_poly c[1]);
unsigned long jac_pp_order_ui (jac_t a[1], unsigned long p, unsigned long h, hc_poly c[1]);
int jac_fastorder_ui (unsigned long *po, jac_t a[1], unsigned long e, hc_poly c[1]);
int jac_factored_order_ui (unsigned long *po, jac_t a[1], unsigned long p[], unsigned long h[], int w, hc_poly c[1]);
unsigned long jac_fastorder_powersmooth (mpz_t o, jac_t a[1], unsigned long L, hc_poly c[1]);
int jac_sylow (jac_t a[JAC_MAX_GENERATORS], unsigned long ords[JAC_MAX_GENERATORS], unsigned long p, unsigned long E, unsigned long M, unsigned long limit, hc_poly c[1]);
int jac_structure (long m[], hc_poly c[1], unsigned long order, int fExponent);

#ifdef __cplusplus
}
#endif

#endif
