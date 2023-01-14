#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "hecurve.h"
#include "cstd.h"

/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// Genus 2 Jacobian arithmetic formulas taken from the Handbook of Elliptic Curve and Hyperelliptic Curve Cryptography (with several important bug fixes!!)

#if HECURVE_GENUS == 2

#define _deg_u(u)			(_ff_nonzero(u[2])?2:(_ff_nonzero(u[1])?1:0))
#define _deg_v(v)			(_ff_nonzero(v[1])?1:(_ff_nonzero(v[0])?0:-1))
#define _dbg_print_uv(s,u,v)	if ( dbg_level >= DEBUG_LEVEL ) { printf (s);  hecurve_print(u,v); }
#define _print_uv(s,u,v)		{ printf (s);  hecurve_print(u,v); }

void hecurve_g2_square_1 (ff_t u[3], ff_t v[2], ff_t u0, ff_t v0, ff_t f[6]);
void hecurve_g2_square_2_r0 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[6]);
void hecurve_g2_compose_1_2 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[6]);
int hecurve_g2_compose_special (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[3], ff_t f[6]);
int hecurve_g2_compose_d6 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[7], hecurve_ctx_t *ctx);
int hecurve_g2_square_d6 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[7], hecurve_ctx_t *ctx);

/*
	TODO: Clean up code to not use quite so many variables.  Useful for debugging at the moment.
*/
// note that (u,v) may be equal to (u1,v1) or (u2,v2) but we assume overlap is aligned.
int hecurve_g2_compose (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[7], hecurve_ctx_t *ctx)
{
	register ff_t r, inv0, inv1, w0, w1, w2, w3, w4, w5, s0, s1, t1, L0, L1, L2;

#if ! HECURVE_FAST
	_dbg_print_uv ("g2_compose input 1: ", u1, v1);
	_dbg_print_uv ("g2_compose input 2: ", u2, v2);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_compose: invalid input\n");  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_compose: invalid input\n");  exit (0); }
#endif
		
	if ( _ff_nonzero(f[6]) ) return hecurve_g2_compose_d6 (u, v, u1, v1, u2, v2, f, ctx);
	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value of r*s1
		_ff_set (r, ctx->r);
		_ff_set (s0, ctx->s0);
		_ff_set (s1, ctx->s1);
		_ff_set (inv1, ctx->inv1);
///		dbg_printf ("Restored r = %ld, s1 = %ld, s0 = %ld, inv1 = %ld, w1 = %ld\n",_ff_get_ui(r), _ff_get_ui(s1), _ff_get_ui(s0), _ff_get_ui(inv1), _ff_get_ui(w1));
		goto hecurve_g2_compose_2_2_inverted;
	}
	if ( _ff_nonzero(f[4]) ) { hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);  return 1; }		// revert to Cantor if f4 is non-zero - should only happen over F_5, of if f has degree 6
	if ( _ff_zero(u1[2]) || _ff_zero(u2[2]) ) { hecurve_g2_compose_special (u, v, u1, v1, u2, v2, f);  return 1; }

	// 1. compute r = Resultant(u1,u2) and save z1=inv1 and z3=inv0
	_ff_sub(inv1,u1[1],u2[1]);			// reduce inv1
	_ff_qsub(w0,u2[0],u1[0]);			// w0 = z2 (unreduced)
	_ff_multadd(inv0,u1[1],inv1,w0);		// z3 = inv0 = u1[1]*inv1+w0 = u1[1]*z1+z2
	_ff_mult(w1,inv1,inv1); 
	ff_mult(w1,w1,u1[0]);
	_ff_multadd(r,w0,inv0,w1);
///	dbg_printf ("r = %ld\n",_ff_get_ui(r));
//	if ( _ff_zero(r) ) { hecurve_compose_special (u, v, u1, v1, u2, v2, f);  return 1; }		this is checked below
	
	// 2. compute almost inverse of u2 mod u1, i.e. inv = (r/u2) mod u1
	// nothing to do, simply note that inv = inv1x + inv0, i.e. inv1 = z1 and inv0 = z3
///	dbg_printf ("inv = %ldx + %ld\n",_ff_get_ui(inv1), _ff_get_ui(inv0));
	
	// 3. compute s' = rs = ((v1-v2)inv) mod u1
	_ff_qsub(w0,v1[0],v2[0]);
	_ff_mult(w2,inv0,w0);
	_ff_qsub(w1,v1[1],v2[1]);
	_ff_mult(w3,inv1,w1); 
	_ff_incmult(t1,u1[1],w3);
	_ff_qadd (w4,inv0,inv1);
	_ff_qadd (w5,w0,w1);
	_ff_qnegsum(w0,w2,t1);
	_ff_multadd(s1,w4,w5,w0);
	_ff_mult(t1,u1[0],w3);
	_ff_qsub(s0,w2,t1);
	
	// Note that we need to save w3, w4, and w5, but w0, w1, and w2 are free
//	if ( _ff_zero(s1) ) goto hecurve_compose_s1_zero;								this is checked below

	// 4. compute s'' = x+s0/s1 = x + s'0/s'1 and s1
	_ff_mult(w2,r,s1);			// assume r and s1 both non-zero - this is the usual case
	if ( _ff_zero(w2) ) {
		if ( _ff_nonzero(r) ) goto hecurve_g2_compose_s1_zero;
		hecurve_g2_compose_special (u, v, u1, v1, u2, v2, f);  return 1;
	}
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->s0, s0);
		_ff_set (ctx->s1, s1);
		_ff_set (ctx->r, r);
		_ff_set(ctx->inv1, inv1);
		_ff_set (ctx->invert, w2);	// return to let caller invert
		ctx->state = 1;
		return 0;
	}
	_ff_invert(w1,w2);
	
hecurve_g2_compose_2_2_inverted:	// caller returns here after inverting
								// we assume w1, r, s0, s1, and inv1 are set, all others free, s0 is unreduced, the rest are reduced
	_ff_mult(w2,w1,r);
	_ff_square(w3,s1); 
	ff_mult(w3,w3,w1);
	_ff_mult(w4,r,w2); 
	_ff_square(w5,w4);
	_ff_mult(s0,s0,w2);
///	dbg_printf ("s\'\' = x + %ld, s1 = %ld\n",_ff_get_ui(s0), _ff_get_ui(w3));
	
	// 5. compute L = s''0*u2 = x^3 + L2*x^2 + L1*x + L0
	_ff_qadd(L2,u2[1],s0);		// ok to leave unreduced
	_ff_multadd(L1,u2[1],s0,u2[0]);
	_ff_mult(L0,u2[0],s0);
///	dbg_printf ("L = x^3 + %ldx^2 + %ldx + %ld\n",_ff_get_ui(L2), _ff_get_ui(L1), _ff_get_ui(L0));
	
	// 6. compute u = (s(L+2v2)-t)/u1 = x^2 + u[1]*x + u[0]		note h = 0
	_ff_qsub(w0,s0,u1[1]);
	_ff_qsub(w1,s0,inv1);
	_ff_mult(w2,w0,w1);			// w2 = (s0-u1[1])(s0-inv1)
	_ff_qadd(w1,u2[1],u2[1]);
	_ff_qaddto(w1,inv1);
	ff_mult(w1,w1,w5);
	_ff_qaddto (w2,w1);
	_ff_qneg(w1,u1[0]);
	_ff_qaddto(w2,w1);
	_ff_addto(w2,L1);
	_ff_qadd(w0,v2[1],v2[1]);
	_ff_multadd(u[0],w0,w4,w2);				// potentially destroys u1[0] and u2[0]
	// no mults to work with so we are better off reducing as we go to compute u[1]
	_ff_add(w0,s0,s0);
	_ff_subfrom(w0,inv1);
	_ff_sub(u[1],w0,w5);
	_ff_set_one(u[2]);
///	dbg_printf ("u = %ldx^2 + %ldx + %ld\n",_ff_get_ui(u[2]), _ff_get_ui(u[1]), _ff_get_ui(u[0]));
	
	// 7. compute v = (-(L+v2)mod u = v[1]*x + v[0]					note h = 0
	_ff_qsub(w1,L2,u[1]);
	_ff_qsub(w0,u[0],L1);
	_ff_multadd(w2,w1,u[1],w0);
	_ff_qneg(w0,v2[1]);
	_ff_multadd (v[1],w2,w3,w0);
	_ff_qneg(w0,L0);
	_ff_multadd(w2,w1,u[0],w0);
	_ff_qneg(w0,v2[0]);
	_ff_multadd(v[0],w2,w3,w0);
///	dbg_printf ("v = %ldx + %ld\n",_ff_get_ui(v[1]), _ff_get_ui(v[0]));
#if ! HECURVE_FAST
	_dbg_print_uv ("Compose_2 Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;

	// Special case is not optimized
hecurve_g2_compose_s1_zero:
///		dbg_printf ("Special case, s\'\'1 == 0\n");

	// 4'. compute s
	_ff_invert(inv0,r);
	ff_mult(s0,s0,inv0); 
///	dbg_printf ("s = %ld\n",_ff_get_ui(s0));
	
	// Part of 6' moved here incase u2[0] gets overwritten below
	_ff_mult(w2,u2[0],s0);
	_ff_addto(w2,v2[0]);

	// 5'. compute u = (t-s(L+2v2))/u1 = x + u[0]		note h = 0
	_ff_square(w1,s0);
	_ff_addto(w1,u2[1]);
	_ff_addto(w1,u1[1]);
	_ff_neg(u[0],w1);				// this potentially destroys u1[0] and u2[0], we assume u2[1] is ok
///	dbg_printf ("u = x + %ld\n",_ff_get_ui(u[0]));
	
	// 6. compute v = -(L+v2) mod u = v[0]				note h = 0
	_ff_sub (w1,u2[1],u[0]);				// BUG FIX - handbook has w1 = s0(u21+u[0]) (+ SHOULD BE -)
	ff_mult(w1,w1,s0);
	_ff_addto(w1,v2[1]);
	ff_mult(w1,w1,u[0]);
	_ff_sub(v[0],w1,w2);
	_ff_set_zero(v[1]);

	_ff_set_zero(u[2]);  _ff_set_one(u[1]);		// set these last to avoid overlap issues

///	dbg_printf ("v = %ld\n",_ff_get_ui(v[0]));
#if ! HECURVE_FAST
	_dbg_print_uv ("Compose_2 Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}


int hecurve_g2_square (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[7], hecurve_ctx_t *ctx)
{
	register ff_t w0, w1, w2, w3, w4, w5, inv0, inv1, L0, L1, L2, r, s0, s1, t0, t1;
	
#if ! HECURVE_FAST
	_dbg_print_uv ("Square_2: ", u1, v1);
#endif	
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_square: invalid input\n");  exit (0); }
#endif

	if ( _ff_nonzero(f[6]) ) return hecurve_g2_square_d6 (u, v, u1, v1, f, ctx);
	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value of r*s1
		_ff_set (r, ctx->r);
		_ff_set (s0, ctx->s0);
		_ff_set (s1, ctx->s1);
///		dbg_printf ("Restored r = %ld, s1 = %ld, s0 = %ld, w1 = %ld\n",_ff_get_ui(r), _ff_get_ui(s1), _ff_get_ui(s0), _ff_get_ui(w1));
		goto hecurve_g2_square_2_2_inverted;
	}
	if ( _ff_nonzero(f[4]) ) {hecurve_compose_cantor (u, v, u1, v1, u1, v1, f);  return 1; }		// revert to Cantor if f4 is non-zero - should only happen over F_5
	if ( _ff_zero(u1[2]) ) { hecurve_g2_square_1 (u, v, u1[0], v1[0], f);  return 1; }

	// 1. compute v' = 2v mod u = v1*x + v0			note h = 0, we use w5=~v1, inv0=~v0
	_ff_add(w5,v1[1],v1[1]);				// w5 gets negated below, must be reduced
	_ff_qadd(inv0,v1[0],v1[0]);
///	dbg_printf ("2v = %ldx + %ld\n",_ff_get_ui(w5), _ff_get_ui(inv0));
	
	// 2. compute resultant r = Resultant(2v,u)
	_ff_square(w0,v1[1]);
	_ff_square(w1,u1[1]);
	_ff_qadd(w2,w0,w0); 
	_ff_qadd(w2,w2,w2);
	_ff_mult(w3,u1[1],w5); 
	_ff_qsub(w4,inv0,w3);
	ff_mult(w4,w4,inv0);
	_ff_multadd(r,u1[0],w2,w4);
//	_ff_addto(r,w4);
///	dbg_printf ("w0 = %ld, w1 = %ld, w2 = %ld, w3 = %ld, w4 = %ld r = res(2v,u) = %ld\n",
///			  _ff_get_ui(w0), _ff_get_ui(w1), _ff_get_ui(w2), _ff_get_ui(w3), _ff_get_ui(w4), _ff_get_ui(r));
//	if ( _ff_zero(r) ) { hecurve_square_2_r0 (u, v, u1, v1, f);  return 1; }		handled below
	
	// 3. compute almost inverse invv = r*inv
	_ff_neg(inv1,w5);
	_ff_subfrom(inv0,w3);
///	dbg_printf ("inv = %ldx + %ld\n", _ff_get_ui(inv1), _ff_get_ui(inv0));

	// 4. compute t = ((f-v^2)/u) mod u = tt1*x+tt0	note h = 0
	// this code could be improved further when f[3] == 0
#if HECURVE_SPARSE
	_ff_set (w3,w1);	// this could be further optimized
#else
	_ff_add(w3,f[3],w1);
#endif
	_ff_add(w4,u1[0],u1[0]);
	_ff_qadd(t1,w1,w1);
	_ff_qaddto(t1,w3);
	_ff_qsubfrom(t1,w4);
	_ff_qadd(t0,w4,w4);
	_ff_qsubfrom(t0,w3);
	_ff_qneg(w5,w0);
	_ff_multadd(t0,t0,u1[1],w5);
#if ! HECURVE_SPARSE
	_ff_addto(t0,f[2]);
#endif
///	dbg_printf ("t = %ldx + %ld\n", _ff_get_ui(t1), _ff_get_ui(t0));  

	// 5. compute ss = (tt*invv) mod u
	_ff_mult(w0,t0,inv0);
	_ff_mult(w1,t1,inv1);
	_ff_add(s1,inv0,inv1);		// need s1 reduced to keep mult from getting too big
	_ff_add(w2,t0,t1);
	_ff_qneg(w5,w0);
	_ff_multadd(s1,s1,w2,w5);
	_ff_incmult(w2,u1[1],w1);
	_ff_qsubfrom(s1,w2);
	_ff_mult(w5,w1,u1[0]);
	_ff_qsub(s0,w0,w5);
///	dbg_printf ("s' = %ldx + %ld\n",_ff_get_ui(s1), _ff_get_ui(s0));
//	if ( _ff_zero(s1) ) goto hecurve_square_s1_zero;		handled below

	// 6. compute s'' = x + s0/s1 and s1
	_ff_mult(w0,r,s1);		// assume non-zero r and s1 is usual case
	if ( _ff_zero(w0) ) {
		if ( _ff_nonzero(r) ) goto hecurve_g2_square_s1_zero;
		hecurve_g2_square_2_r0 (u, v, u1, v1, f);
		return 1;
	}
	if ( ctx && ! ctx->state  ) {
		_ff_set (ctx->s0, s0);
		_ff_set (ctx->s1, s1);
		_ff_set (ctx->r, r);
		_ff_set (ctx->invert, w0);
		ctx->state = 1; 								// return to let caller invert
		return 0;
	}
	_ff_invert(w1,w0);

hecurve_g2_square_2_2_inverted:	// caller returns here after inverting
							// we assume w0, s0, s1, and r are set and everything else is free
	_ff_mult(w2,w1,r); 
	_ff_square(w3,s1);
	ff_mult(w3,w3,w1);
	_ff_mult(w4,w2,r);
	_ff_square(w5,w4);
	ff_mult(s0,s0,w2);
///	dbg_printf ("s\'\'0 = %ld, w0 = %ld, w1 = %ld, w2 = %ld, w3 = %ld, w4 = %ld, w5 = %ld\n",
///			  _ff_get_ui(s0), _ff_get_ui(w0), _ff_get_ui(w1), _ff_get_ui(w2), _ff_get_ui(w3), _ff_get_ui(w4), _ff_get_ui(w5));
	// note that w3, w4, and w5, need to be saved, but w0, w1 and w2 are free
		
	// 7. compute LL = sss*u = x^3 + LL2*x^2 + LL1*x + LL0
	_ff_qadd(L2,s0,u1[1]);
	_ff_multadd(L1,s0,u1[1],u1[0]);
	_ff_mult(L0,s0,u1[0]);
///	dbg_printf ("L' = s''u = x^3 + %ldx^2 + %ldx + %ld\n",_ff_get_ui(L2), _ff_get_ui(L1), _ff_get_ui(L0));
	
	// 8. compute u = s^2 + 2vs/u + (v^2-f)/u1^2		note h = 0
	_ff_square(w0,s0);
	_ff_qadd(w1,v1[1],v1[1]);
	ff_mult(w1,w1,w4);
	_ff_qaddto(w0,w1);
	_ff_qadd(w1,u1[1],u1[1]);
	_ff_multadd(u[0],w1,w5,w0);						// potentially destroys u1[0] and u2[0]
	_ff_add(w0,s0,s0);								// no mult to absorb so stay reduced
	_ff_sub(u[1],w0,w5);							// potentially destroys u1[1] and u2[1]
	_ff_set_one(u[2]);
///	dbg_printf ("u = x^2 + %ldx + %ld\n",_ff_get_ui(u[1]), _ff_get_ui(u[0]));
	
	// 9. compute v = -(L+v1) mod u = v[1]x + v[0]			note h = 0
	_ff_qsub(w1,L2,u[1]);
	_ff_qsub (w0,u[0],L1);
	_ff_multadd(w2,w1,u[1], w0);
	_ff_qneg(w0,v1[1]);		// v1[1] could overlap v[1]
	_ff_multadd(v[1],w2,w3,w0);
	_ff_qneg(w0,L0);
	_ff_multadd(w2,w1,u[0],w0);
	_ff_qneg(w0,v1[0]);		// v1[0] could overlap v[0]
	_ff_multadd(v[0],w2,w3,w0);
///	dbg_printf ("v = %ldx + %ld\n",_ff_get_ui(v[1]), _ff_get_ui(v[0]));
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_square: output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that input has been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;

hecurve_g2_square_s1_zero:
///	dbg_printf ("Special case, s'1 == 0\n");
	// 6'. compute s and precomputations
	_ff_invert(w1,r);
	ff_mult(s0,s0,w1);
	_ff_mult(w2,s0,u1[0]);
	_ff_addto(w2,v1[0]);
///	dbg_printf ("w1 = %ld, w2 = %ld, s0 = %ld\n",_ff_get_ui(w1), _ff_get_ui(w2), _ff_get_ui(s0));

	// 7'. compute u' = (f-v^2)/u^2 - 2vs/u - s^2			note h =0
	_ff_square(w3,s0);
	_ff_add(w4,u1[1],u1[1]);
	_ff_addto(w3,w4);					
	_ff_neg(u[0],w3);							// potentially destroys u1[0] and u2[0]
///	dbg_printf ("u = x + %ld\n",_ff_get_ui(u[0]));

	// 8'. compute v' = -(su+v) mod u
	_ff_sub(w1,u1[1],u[0]);
	ff_mult(w1,w1,s0);
	_ff_addto(w1,v1[1]);
	ff_mult(w1,w1,u[0]);
	_ff_sub(v[0],w1,w2);
	_ff_set_zero(v[1]);
///	dbg_printf ("v = %ld\n",_ff_get_ui(v[0]));

	_ff_set_one(u[1]);							// set these last in case of overlap
	_ff_set_zero(u[2]);
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v,f) ) {
		err_printf ("hecurve_square output: verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that input has been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}

// Note that (u,v) == (u1,v1) or (u2,v2) (or both) are possible.  We do assume that any overlap is aligned.
// Thus modifying u[1] shouldn't effect u1[0] or u2[0] but may destory u1[1] and/or u2[1]
//
// This algorithm handles all the unusual cases - it assumes that either u1 or u2 has degree < 2
// or that Resultant(u1,u2) = 0, otherwise hecurve_compose would have handled it
int hecurve_g2_compose_special (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[3], ff_t f[6])
{
	register ff_t t1, t2, t3;
	ff_t *swap;
	int d1, d2, swapd;

#if ! HECURVE_FAST
	gmp_dbg_printf ("Compose Special: (%ldx^2+%ldx+%ld, %ldx+%ld) * (%ldx^2+%ldx+%ld, %ldx+%ld)\n",
			    _ff_get_ui(u1[2]), _ff_get_ui(u1[1]), _ff_get_ui(u1[0]), _ff_get_ui(v1[1]), _ff_get_ui(v1[0]),
			    _ff_get_ui(u2[2]), _ff_get_ui(u2[1]), _ff_get_ui(u2[0]), _ff_get_ui(v2[1]), _ff_get_ui(v2[0]));
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u1, v1, f) )  { err_printf ("hecurve_compose_special: invalid input\n");  exit (0); }
	if ( ! hecurve_verify (u2, v2, f) ) { err_printf ("hecurve_compose_special: invalid input\n");  exit (0); }
#endif
	d1 = _deg_u(u1);
	d2 = _deg_u(u2);
	if ( d1 > d2 ) {
		swap = u2;  u2 = u1;  u1 = swap;
		swap = v2;  v2 = v1;  v1 = swap;
		swapd = d2;  d2 = d1;  d1 = swapd;
	}
	switch(d1){
	case 0:		// deg u1 == 0, i.e. u1 is the identity
///		dbg_printf ("Case 1, deg(u1) == 0\n");
		_hecurve_set (u, v, u2, v2);
		break;
	case 1: 	// deg u1 == 1
		// case 2
///		dbg_printf ("Case 2, deg(u1) == 1\n");
		if ( d2 == 1 ) { // deg u1 == dege u2 == 1
///			dbg_printf ("Case 2.A, deg(u1) == deg(u2) == 1\n");
			// case 2.A
			if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) ) {
				if ( _ff_zero(v1[0]) && _ff_zero(v2[0]) ) {
///					dbg_printf ("Case 2.A.i.a D1 == -D2 == D1\n");
					_hecurve_set_identity(u,v);
					break;
				}
				_ff_neg(t1,v2[0]);
				if ( _ff_equal(v1[0],t1) ) {
///					dbg_printf ("Case 2.A.i.b D1 == -D2\n");
					_hecurve_set_identity(u,v);
					break;
				}
				if ( _ff_equal(v1[0],v2[0]) ) {
///					dbg_printf ("Case 2.A.ii D1 == D2\n");
					hecurve_g2_square_1 (u, v, u1[0], v1[0], f);
					break;
				}
			}
///			dbg_printf ("Case 2.A.iii D1 != +- D2\n");
			if ( _ff_equal(u1[0],u2[0]) ) { err_printf ("u1[0] == u2[0] in case 2.A.iii of hecurve_compose!\n");  exit(0); }
			
			// u = u1*u2 and v = ((v2-v1)x+v2u1[0]-v1u2[0])/(u1[0]-u2[0])		note that v1 and v2 are constants (deg 0)
			_ff_sub(t1,u1[0],u2[0]);
			ff_invert(t1,t1);
			_ff_sub(t2,v2[0],v1[0]);
			_ff_mult(v[1],t2,t1);				// we assume writing [1] doesn't hurt [0]
			_ff_mult(t2,v2[0],u1[0]);
			_ff_mult(t3,v1[0],u2[0]);
			_ff_subfrom(t2,t3);
			_ff_mult(v[0],t2,t1);
			_ff_set_one(u[2]);					// update u last to handle overlap
			_ff_add(u[1],u1[0],u2[0]); 		// we assume writing [1] doesn't hurt [0]
			ff_mult(u[0],u1[0],u2[0]);		// mult should be safe
			break;
		} else { // deg u1 == 1, deg u2 == 2
///			dbg_printf ("Case 2.B, deg(u1) == 1, deg(u2) == 2\n");
			// compute u2(-u1[0]) = (-u1[0])(-u1[0] + u2[1]) + u2[0]
			_ff_neg(t1,u1[0]);				// save t1 = -u1[0] for later
			_ff_add(t2,t1,u2[1]);
			ff_mult(t2,t2,t1);
			_ff_addto(t2,u2[0]);
			if ( _ff_nonzero(t2) ) {
///			dbg_printf ("Case 2.B.i, u2(-u1[0]) != 0\n");
				hecurve_g2_compose_1_2 (u, v, u1, v1, u2, v2, f);
				break;
			}
///			dbg_printf ("Case 2.B.ii\n");
			_ff_add(t2,u1[0],u1[0]);
			if ( _ff_equal(u2[1],t2) ) {
				_ff_square(t2,u1[0]);
				if ( _ff_equal(u2[0],t2) ) {
///					dbg_printf ("Case 2.B.ii.a, u2[1] == 2u1[0] and u2[0] == u1[0]^2\n");
					// compute v2(-u1[0])
					_ff_mult(t2,v2[1],t1);
					_ff_addto(t2,v2[0]);
					if ( _ff_equal(t2,v1[0]) ) {
						ff_t u3[3], v3[2], v4[2];
										
///						dbg_printf ("Case 2.B.ii.a.I D2 = 2D1\n");
						// Double D2 then subtract D1
						// this is slightly inefficient but easier and is a rare case anyway
						hecurve_square (u3, v3, u2, v2, f, 0);  // double D2
						// Bug fix - need to avoid infinite recursion on order 3 elements, this occurs when 2D2 == -D2
						if ( _ff_equal(u3[0],u2[0]) && _ff_equal(u3[1],u2[1]) && _ff_equal(u3[2],u2[2]) ) {
							_ff_neg(t2,v2[0]);
							_ff_neg(t3,v2[1]);
							if ( _ff_equal(v3[0],t2) && _ff_equal(v3[1],t3) ) {
								hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);
								break;
//								err_printf ("Attempt to compute 2D+D for D with order 6 - unable to handle this case\n");
//								exit (0);
							}
						}
///						dbg_printf ("2D2 != -2D1, computing 2D2-D1\n");
						_ff_neg(v4[0],v1[0]);
						_ff_neg(v4[1],v1[1]);						// D4 = invert D1, note u4 = u1
						hecurve_compose (u, v, u1, v4, u3, v3, f, 0);		// compute D2+D2-D1
					} else {
///						dbg_printf ("Case 2.B.ii.a.II, D2 = -2P1\n");
						_hecurve_set (u, v, u1, v1);
						hecurve_invert (u, v);						// compute inverse of D1
					}
					break;
				}
			} else {
///				dbg_printf ("Case 2.B.ii.b, u2[1] != 2u1[0] or u2[0] != u1[0]^2\n");
				// compute -v2(-u1[0]) - this is a bug fix - p. 314/315 of handbook uses v2(-u1[0]) (or not)
				_ff_mult(t3,v2[1],t1);		// t1 = -u1[0]
				_ff_addto(t3,v2[0]);
				_ff_neg(t2,t3); 
				// we need to compute v2(u1[0]-u2[1]) either way
				_ff_sub(t3,u1[0],u2[1]);
				ff_mult(t3,t3,v2[1]);
				_ff_addto(t3,v2[0]);
///				dbg_printf ("t3 = v2[1]*(u1[0]-u2[1]) + v2[0] = %ld\n",_ff_get_ui(t3));
				if ( _ff_equal(t2,v1[0]) ) {
///					dbg_printf ("Case 2.B.ii.b.I, -P1 occurs in D2\n");
					_ff_add(u[0],u2[1],t1);						// t1 = -u1[0]
					_ff_set_zero(u[2]);  _ff_set_one(u[1]);				// u[1] must get set here, could destory u2[1]
					_ff_set_zero(v[1]);
					_ff_set(v[0],t3);
				} else {
					ff_t u3[3], v3[2], u4[3], v4[2];
					
///					dbg_printf ("Case 2.B.ii.b.II, Doubling D1 and composing with (x+u21-u10, v2(u10-u21)),\n");
					hecurve_g2_square_1 (u3, v3, u1[0], v1[0], f);
					_ff_set_zero(u4[2]);  _ff_set_one(u4[1]);
					_ff_sub(u4[0],u2[1],u1[0]);
					_ff_set_zero(v4[1]);
					_ff_set(v4[0],t3);
					hecurve_compose (u, v, u4, v4, u3, v3, f, 0);
				}
			}
		}
		break;
	case 2: // deg u1 == deg u2 == 2 (case 3)
///		dbg_printf ("Case 3, deg(u1) == deg(u2) == 2\n");
		if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) ) {
///			dbg_printf ("Case 3.A, u1 == u2\n");
			if ( _ff_equal(v1[1],v2[1]) && _ff_equal(v1[0],v2[0]) ) {
				// catch inverse case for order 2 elements
				if ( _ff_zero(v1[1]) && _ff_zero(v1[0]) ) {
///					dbg_printf ("Case 3.A.i*, v1 == -v2 == 0\n");
					_hecurve_set_identity (u, v);
					break;
				}
///				dbg_printf ("Case 3.A.ii, v1 == v2\n");
				hecurve_square (u, v, u1, v1, f, 0);				// handles case 3.A.ii.b (zero resultant) also
				break;
			}
			_ff_neg(t1,v2[0]);
			_ff_neg(t2,v2[1]);
			if ( _ff_equal(v1[0],t1) && _ff_equal(v1[1],t2) ) {
///				dbg_printf ("Case 3.A.i, v1 == -v2\n");
				_hecurve_set_identity (u, v);
				break;
			}
///			dbg_printf ("Case 3.A.iii v1 != +- v2\n");
			// case 3.A.iii
			_ff_sub(t1,v2[1],v1[1]);
#if ! HECURVE_FAST
			if ( _ff_zero(t1) ) { err_printf ("v2[1] == v1[1] in case 3.A.iii of hecurve_compose_special\n");  exit (0); }
#endif
			_ff_invert(t2,t1);
			_ff_sub(t1,v1[0],v2[0]);
			ff_mult(t1,t1,t2);
			_ff_neg(t2,t1);
			_ff_mult(t3,v1[1],t1);
			_ff_addto(t3,v1[0]);
			hecurve_g2_square_1 (u, v, t2, t3, f);
		} else { // u1 != u2
///			dbg_printf ("Case 3.B, u1 != u2\n");
			// case 3.B
			// We can assume Resultant(u1,u2) == 0, since otherwise we wouldn't have been called
			ff_t x2, y2, x3, y3;
					
///			dbg_printf ("Case 3.B.ii, res(u1,u2) == 0\n");
			// need to compute gcd(u1,u2) = x-x1 where x1 = (u1[1]-u2[1])/(u2[0]-u1[0]))
			_ff_sub(t1,u1[1],u2[1]);
#if ! HECURVE_FAST
			if ( _ff_zero(t1) ) { err_printf ("u1[1] == u2[1] in case 3.B.ii of hecurve_compose_special\n");  exit (0); }
#endif
			_ff_invert(t2,t1);
			_ff_sub(t1,u2[0],u1[0]);
			_ff_mult(t3,t1,t2);
///			dbg_printf ("gcd(u1,u2) = x-%ld\n", _ff_get_ui(t3));
			// compute v1(x1) and v2(x1) and compare
			_ff_mult(t1,v1[1],t3);
			_ff_addto(t1,v1[0]); 
			_ff_mult(t2,v2[1],t3);
			_ff_addto(t2,v2[0]); 
			// Need to extract coords of P2 and P3 in either case, so compute them here
			// This code could be cleaned up a bit - there is some double negation, but it is a special case...
			_ff_add(y2,u1[1],t3); 
			_ff_neg(x2,y2);  
			_ff_mult(y2,v1[1],x2);
			_ff_addto(y2,v1[0]);
			_ff_add(y3,u2[1],t3);
			_ff_neg(x3,y3); 
			_ff_mult(y3,v2[1],x3);
			_ff_addto(y3,v2[0]); 
///			dbg_printf ("P2 = (%ld,%ld), P3 = (%ld,%ld)\n", _ff_get_ui(x2), _ff_get_ui(y2), _ff_get_ui(x3), _ff_get_ui(y3));
			if ( _ff_equal(t1,t2) ) {
				ff_t u3[3], v3[2], u4[3], v4[2], u5[3], v5[2];
							
///					dbg_printf ("Case 3.B.ii.a, v1(x1) == v2(x1)\n");
				_ff_neg(t1,t3);
				_ff_mult(t2,v1[1],t3);
				_ff_add(t2,t2,v1[0]);
				hecurve_g2_square_1 (u4, v4, t1, t2, f);		// D' = 2(P1)
///				dbg_printf ("D' = (%ldx^2 + %ldx + %ld, %ldx + %ld)\n",_ff_get_ui(u4[2]), _ff_get_ui(u4[1]), _ff_get_ui(u4[0]), _ff_get_ui(v4[1]), _ff_get_ui(v4[0]));
				// could we use hecurve_compose_1_2 here?
				_ff_set_zero(u3[2]);  _ff_set_one(u3[1]);
				_ff_neg(u3[0],x2);
				_ff_set_zero(v3[1]);
				_ff_set(v3[0],y2);
				hecurve_compose (u5, v5, u3, v3, u4, v4, f, 0);		// D'' = D' + (P2)
///				dbg_printf ("D'' = (%ldx^2 + %ldx + %ld, %ldx + %ld)\n",_ff_get_ui(u5[2]), _ff_get_ui(u5[1]), _ff_get_ui(u5[0]), _ff_get_ui(v5[1]), _ff_get_ui(v5[0]));
				_ff_neg(u3[0],x3);
				_ff_set(v3[0],y3);
				hecurve_compose (u, v, u3, v3, u5, v5, f, 0);		// D = D'' + (P3)
			} else {
///				dbg_printf ("Case 3.B.ii.a, v1(x1) != v2(x1)\n");
				hecurve_make_2 (u, v, x2, y2, x3, y3);
			}
		}
		break;
	default:
		err_printf ("Invalid degree in hecurve_compose!\n");  exit (0);
	}
#if ! HECURVE_FAST
	_dbg_print_uv ("Compose Special Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose_special output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}


// Algorithm 14.20 p.318
// As above, handle aligned overlap of (u,v) with (u1,v1) and/or (u2,v2)
void hecurve_g2_compose_1_2 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[6])
{
	register ff_t r, inv, t1, t2, w1, L0, L1, s0;

///	dbg_printf ("Compose_1_2\n");	

	// 1. compute r = u2 mod u1
	_ff_sub(w1,u2[1],u1[0]);
	ff_mult(w1,w1,u1[0]);
	_ff_sub(r,u2[0],w1);
///	dbg_printf ("r = %ld\n",_ff_get_ui(r));
	
	// 2. compute inverse of u2 mod u1
	_ff_invert(inv, r);
///	dbg_printf ("inv = %ld\n",_ff_get_ui(inv));
	
	// 3. compute s = ((v1-v2)inv) mod u1
	_ff_mult(w1,v2[1],u1[0]);			// BUG FIX - handbook incorrectly uses -v2[1]u1[0] - should be positive
	_ff_addto(w1,v1[0]);		
	_ff_subfrom(w1,v2[0]); 
	_ff_mult(s0,w1,inv);
///	dbg_printf ("s = %ld\n", _ff_get_ui(s0));
	
	// 4. compute L = su2 = s0*x^2 + L1*x + L0
	_ff_mult(L1,s0,u2[1]);
	_ff_mult(L0,s0,u2[0]);
	// no reduction required for L coeffs
///	dbg_printf ("L = su2 = %ldx^2 + %ldx + %ld\n",_ff_get_ui(s0), _ff_get_ui(L1), _ff_get_ui(L0));
	
	// 5. compute t = (f - v2^2)/u2 = x^3 + t2*x2 + t1*x + t0   note h = 0
	_ff_neg(t2,u2[1]);
	_ff_square(t1,u2[1]);
	_ff_subfrom(t1,u2[0]);
	if ( _ff_nonzero(f[3]) ) _ff_addto(t1,f[3]);
///	dbg_printf ("t = x^3 + %ldx^2 + %ldx + t0\n",_ff_get_ui(t2), _ff_get_ui(t1));
	
	// 6. compute u = (t-s(L+2v2))/u1 = x^2 + u[1]*x + u[0]	 	note h = 0
	_ff_set_one(u[2]);
	_ff_square(r,s0);
	_ff_sub(w1,t2,r);
	_ff_sub(u[1],w1,u1[0]);			// ther should be no overlap - potentially destroys u1[1] and u2[1]
	_ff_add(w1,v2[1],v2[1]);
	_ff_addto(w1,L1);
	ff_mult(w1,w1,s0);
	_ff_sub(r,t1,w1);
	_ff_mult(t2,u1[0],u[1]);
	_ff_sub(u[0],r,t2);								// potentially destroys u1[0] and u2[0]
///	dbg_printf ("u = x^2 + %ldx + %ld\n",_ff_get_ui(u[1]), _ff_get_ui(u[0]));
	
	// 7. compute v = (-(L+v2)) mod u = v[1]x + v[0]			note h = 0
	_ff_mult(w1,s0,u[1]);
	_ff_subfrom(w1,L1);
	_ff_subfrom(w1,v2[1]);			// v2[1] could overlap v[1]
	_ff_set (v[1],w1);
	_ff_mult(w1,s0,u[0]);
	_ff_subfrom(w1,L0);
	_ff_subfrom(w1,v2[0]);			// v2[0] could overlap v[0]
	_ff_set(v[0],w1);
///	dbg_printf ("v = %ldx + %ld\n",_ff_get_ui(v[1]), _ff_get_ui(v[0]));
}


// This code would benefit from precomputing f' - todo later
// Note that overlap of u[0] and u0 (and v[0] and v0) is possible!
void hecurve_g2_square_1 (ff_t u[3], ff_t v[2], ff_t u0, ff_t v0, ff_t f[6])
{
	register ff_t x, y, t1, t2, t3;
	
	if ( _ff_zero(v0) ) { _hecurve_set_identity (u, v);  return; }
	// The code below corresponds to (14.9) on p. 314 of Handbook of E+HE Crypto
///	dbg_printf ("Square_1 (x + %ld,%ld)\n",_ff_get_ui(u0),_ff_get_ui(v0));
	// u = u1^2
	_ff_set (t3, u0);				// save t3=u0 since we need it later and may overwrite it here.
	_ff_set_one(u[2]);
	_ff_add(u[1],t3,t3);
	_ff_square(u[0],t3);
///	dbg_printf ("u = x^2 + %ldx + %ld\n",_ff_get_ui(u[1]), _ff_get_ui(u[0]));
	// v = (f'(-u1[0])x + f'(-u1[0])u1[0])/(2v1) + v1
	// Compute y = f'(-u1[0])  - note that this code assumes that 5* max coeff of f fits in an ff_t
	_ff_neg(x,t3);
	_ff_set(y,f[1]);
	if ( _ff_nonzero(f[2]) ) {
		_ff_add(t2,f[2],f[2]);
		ff_mult(t2,t2,x);
		_ff_addto(y,t2);
	}
	_ff_square(t1,x);
	if ( _ff_nonzero(f[3]) ) {
		_ff_add(t2,f[3],f[3]); 
		_ff_addto(t2,f[3]);
		ff_mult(t2,t2,t1);
		_ff_addto(y,t2);
	}
#if FF_WORDS == 1
	if ( _ff_p != 5 ) {
#endif
		ff_square(t1,t1);
		_ff_add(t2,t1,t1); 
		_ff_x2(t2); 
		_ff_addto(t2,t1);
		_ff_addto(y,t2); 
#if FF_WORDS == 1
	}
#endif
///	dbg_printf ("f'(-u0) = %ld\n",_ff_get_ui(y));
	_ff_add(t2,v0,v0);
	_ff_invert(t1,t2);
	_ff_mult(v[1],t1,y);
	_ff_mult(t1,v[1],t3);
	_ff_add(v[0],t1,v0);
///	dbg_printf ("v = %ldx + %ld\n",_ff_get_ui(v[1]), _ff_get_ui(v[0]));
}

/*
	This algorithm handles squaring a degree 2 element where Resultant(2v,u) = 0.
	In this case the divisor contains a point with order 2 and we simply want to square the other point.
	This corresponds to case 3.A.ii.b in the general algorithm, but may also be detected in hecurve_square
*/
void hecurve_g2_square_2_r0 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[6])
{
	register ff_t t1, t2;

	// Since h = 0, we simply need gcd(2v1, u1) which must be (x+v1[0]/v1[1]) hence x1 = -v1[0]/v1[1] is the x-coord of P1
	// We want to then double P2 = (-u1[1]-x1,v1(-u1[1]-x1))  see p. 315 of the handbook
	if ( _ff_zero(v1[1]) ) {
		if ( _ff_zero(v1[0]) ) { _hecurve_set_identity (u, v);  return; }
		// this should be impossible if Results(u,2v) = 0
		err_printf ("v1[1] = 0 and v1[0] != 0 in hecurve_square_2_r0!");
		exit (0);		
	}
	_ff_invert(t1,v1[1]);
	_ff_neg(t2,v1[0]);
	ff_mult(t1,t1,t2);
	_ff_addto(t1,u1[1]);
	_ff_neg(t2,t1);
	ff_mult(t2,t2,v1[1]);
	_ff_addto(t2,v1[0]);
	hecurve_g2_square_1 (u, v, t1, t2, f);
}

/*
	IMPORTANT: the degree 6 functions assume that the curve has no rational points at infinity,
	equivalently, f[6] is a nonresidue.  Note that every curve is birational to such a curve, so this
	assumption is not restrictive, but it does require the curve to be put in a suitable form.
*/

int hecurve_g2_compose_d6 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t u2[3], ff_t v2[2], ff_t f[7], hecurve_ctx_t *ctx)
{
	register ff_t z0, z1, w0, w1, w2, w3, s0, s1, r, d, t2, t3, t4, ell0, ell1, ell2;
	
	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value 
		_ff_set (r, ctx->r);
		_ff_set (s0, ctx->s0);
		_ff_set (s1, ctx->s1);
		_ff_set (w3, ctx->inv1);
///		dbg_printf ("Restored r = %ld, s1 = %ld, s0 = %ld, inv1 = %ld, w1 = %ld\n",_ff_get_ui(r), _ff_get_ui(s1), _ff_get_ui(s0), _ff_get_ui(inv1), _ff_get_ui(w1));
		goto hecurve_g2_compose_d6_inverted;
	}

	if ( _ff_zero(u1[2]) || _ff_zero(u2[2]) ) {
		if ( _hecurve_is_identity(u1,v1) ) { _hecurve_set(u,v,u2,v2); return 1; }
		if ( _hecurve_is_identity(u2,v2) ) { _hecurve_set(u,v,u1,v1); return 1; }
		hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);  return 1;
	}
	if ( _ff_equal(u1[0],u2[0]) && _ff_equal(u1[1],u2[1]) ) {
		if ( _ff_equal(v1[0],v2[0]) && _ff_equal(v1[1],v2[1]) ) { return hecurve_g2_square_d6 (u,v,u1,v1,f,0); }
		_ff_add(w0,v1[0],v2[0]); _ff_add(w1,v1[1],v2[1]);
		if ( _ff_zero(w0) && _ff_zero(w1) ) { _hecurve_set_identity(u,v); return 1; }
		hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);  return 1;
	}
	
	// compute inv = z1*x+z0 = r/u2 mod u1 (3M+1S)
	_ff_sub(z1,u1[1],u2[1]);				// z1 = u11-u21
	_ff_sub(w0,u2[0],u1[0]);				// w0=u20-u10
	_ff_multadd (z0,u1[1],z1,w0);			// z0 = u11*z1+w0
	_ff_square(w2,z1);					// w2=z1^2
	_ff_sum_2_mults (r,w0,w2,u1[0],z0);	// r = w0*z0+w2*u10
	if ( _ff_zero(r) ) { hecurve_compose_cantor (u, v, u1, v1, u2, v2, f);  return 1; }
	
	// compute s' = s1*x+s0 = rs = r((v1-v2)/u2) mod u1  (5M)
	_ff_sub(w0,v1[0],v2[0]);				// w0 = v10-v20
	_ff_sub(w1,v1[1],v2[1]);				// w1 = v11-v21
	_ff_mult(w2,z0,w0);					// w2 = z0*w0
	_ff_mult(w3,z1,w1);					// w3 = z1*w1
	_ff_addto(w0,w1);			
	_ff_addto(z0,z1);
	_ff_neg(w1,u1[1]);
	_ff_dec(w1);
	_ff_sum_2_mults(s1,w0,w3,w1,z0);
	_ff_subfrom(s1,w2);					// s1 = (z0+z1)*(w0+w1)-w3*(1+u11)-w2  (Karatsuba)
	_ff_mult(w0,u1[0],w3);
	_ff_sub(s0,w2,w0);					// s0=w2-u10*w3
	
	// compute s = s1*x+s0, 1/r and d=1/(s1^2-f6) (6M, 2S, 1C, I)	TODO: optimize mult by small f6
	_ff_square(t4,r);					// r2 = r^2
	_ff_neg(w0,s1);
	_ff_sum_2_mults (w3,t4,w0,s1,f[6]);	// w3 = r2*f6-s1^2 (note that we don't actually square, it isn't any faster
	_ff_mult (w2,r,w3);
	_ff_mult(r,t4,r);					// r = r^3
	
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->invert, w2);
		_ff_set (ctx->r, r);
		_ff_set (ctx->s0, s0);
		_ff_set (ctx->s1, s1);
		_ff_set(ctx->inv1, w3);
		ctx->state = 1;
		return 0;						// return to let caller invert
	}
	_ff_invert(w1,w2);					// w2 = 1/(r*w2)
	
hecurve_g2_compose_d6_inverted:	// caller returns here after inverting
	
	_ff_mult(d,r,w1);					// d = r^3*w1 = 1/(s1^2-f6)
	_ff_mult(r,w1,w3);					// 1/r = w1*w3 (stored in r)
	_ff_mult(s1,s1,r);					// s1 = r*s1/r
	_ff_mult(s0,s0,r);					// s0 = r*s0/r

	// compute ell = s*u2 - ell3*x^3 + ell2*x^2 + ell1*x + ell0 (but we don't need ell3=s1)  (3M)
	_ff_mult(ell0,s0,u2[0]);				// ell0=s0*u20
	_ff_mult(ell2,s1,u2[1]);				// ell2 = s1*u21 (partial computation)
	_ff_add(w0,s0,s1);
	_ff_add(w1,u2[0],u2[1]);
	_ff_mult(ell1,w0,w1);
	_ff_subfrom(ell1,ell2);
	_ff_subfrom(ell1,ell0);				// ell1 = (s0+s1)*(u20+u21)-ell2-ell0	(Karatsuba)
	_ff_addto(ell2,s0);					// ell2 = s1*u21+s0

	// compute top 3 coefficients of t = t4*x^2+t3*x^3+t2*x^2+O(x) = (f-v2^2)/u2 (1M+2C)		TODO: optimize f6 mults
	_ff_mult(w0,u2[1],f[6]);
	_ff_sub(t3,f[5],w0);					// t4=f6 (implicit) and t3=f5-u21*t4 (we could assume f[5] is zero here, but there is no reason to)
	_ff_sum_2_mults(w0,u2[0],u2[1],t3,f[6]);
	_ff_sub(t2,f[4],w0);					// t2 = f4-u20*f6-u21*t3
	
	// replace t by t-s(ell+2v2)  (3M,1S)
	_ff_square(w0,s1);
	_ff_sub(t4,f[6],w0);					// t4 = f6-s1^2
	_ff_add(w0,s0,ell2);
	_ff_mult(w1,w0,s1);
	_ff_subfrom(t3,w1);					// t3 -= s1*(s0+ell2)
	_ff_add(w0,v2[1],v2[1]);
	_ff_addto(w0,ell1);
	_ff_sum_2_mults(w1,s0,s1,w0,ell2);
	_ff_subfrom(t2,w1);					// t2 -= s0*ell2+s1*(ell1+2*v21)

	// compute u3 = x^2 + u31*x + u30 = (t-s(ell+2v2))/u1 made monic (5M)   Be careful with overlap, we must all u = u1
	_ff_mult(w0,t4,u1[1]);
	_ff_subfrom(t3,w0);					// t3 = u31 = t3-u32*u11 (u32=t4)
	_ff_sum_2_mults(w0,u1[0],u1[1],t3,t4);
	_ff_subfrom(t2,w0);					// t2 = u30 = t2-u31*u11-u32*u10 (u32=t4)
	_ff_set_one(u[2]);
	_ff_mult(u[1],t3,d);					// u31 = d*u31
	_ff_mult(u[0],t2,d);					// u30 = d*u30
	
	// compute v3 = -(ell+v2) mod u3 (3M)
	_ff_addto(ell1,v2[1]);				// ell1 += v21
	_ff_addto(ell0,v2[0]);				// ell0 += v20
	_ff_mult(w1,u[1],s1);				// w1 = u31*s1
	_ff_sub(w2,ell2,w1);					// w2 = ell2-w1
	_ff_mult(w3,u[0],w2);				// w3 = u30*w2
	_ff_add(w0,u[0],u[1]);
	_ff_addto(s1,w2);
	_ff_mult(w0,w0,s1);
	_ff_subfrom(w1,w0);
	_ff_addto(w1,ell1);
	_ff_addto(w1,w3);
	_ff_neg(v[1],w1);					// v31 = -ell1+(u30+u31)(s1+w2)-w1-w3
	_ff_sub(v[0],w3,ell0);				// v30 = we - ell0

	// total cost 29M+3S+3C+I

#if ! HECURVE_FAST
	_dbg_print_uv ("g2_compose_d6 Result: ", u, v);
#endif
#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_compose output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("      ");  hecurve_print(u2,v2);
		err_printf ("note that inputs have been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}


int hecurve_g2_square_d6 (ff_t u[3], ff_t v[2], ff_t u1[3], ff_t v1[2], ff_t f[7], hecurve_ctx_t *ctx)
{
	register ff_t d, w0, w1, w2, w3, w4, w5, ell0, ell1, ell2, r, s0, s1, t0, t1, z0, z1, g0, g1;
	
	if ( ctx && ctx->state == 1 ) {
		// get invert result and restore saved variables
		_ff_set (w1, ctx->invert);		// inverted value
		_ff_set (r, ctx->r);
		_ff_set (s0, ctx->s0);
		_ff_set (s1, ctx->s1);
		_ff_set (t0, ctx->s2);
		_ff_set (w3, ctx->inv1);
		goto hecurve_g2_square_d6_inverted;
	}
	if ( _ff_zero(u1[2]) || (_ff_zero(v1[0]) &&_ff_zero(v1[1])) ) { _hecurve_set_identity (u,v); return 1; }
	if ( _ff_nonzero(f[5]) ) { hecurve_compose_cantor (u, v, u1, v1, u1, v1, f);  return 1; }		// revert to Cantor if f5 is non-zero, shouldn't ever happen
	
	// compute z=z1*x+z0 and r s.t. z=r/(2v) mod u   [3M+2S]
	_ff_add(t0,v1[0],v1[0]);						// t0 = 2v0
	_ff_add(z1,v1[1],v1[1]);						// z1 = 2v1
	_ff_square(w0,v1[1]);						// w0 = v1^2
	_ff_square(w1,u1[1]);						// w1 = u1^2
	_ff_add(w2,w0,w0); _ff_x2(w2);				// w2 = t0^2 = 4w0
	_ff_mult(w3,u1[1],z1);						// w3 = u1*z1
	_ff_sub(z0,t0,w3);							// z0 = t0-w3
	_ff_sum_2_mults(r,u1[0],t0,z0,w2);				// r = u0*w2+t0*z0
	if ( _ff_zero(r) ) { hecurve_compose_cantor (u, v, u1, v1, u1, v1, f);  return 1; }
	ff_negate(z1);
//printf ("r=%ld, z=%ld*x+%ld\n", _ff_get_ui(r), _ff_get_ui(z1), _ff_get_ui(z0));

	// compute t'=t1*x+t0 s.t. t' = ((f-v^2)/u) mod u   [5M+2C]
	_ff_mult(w2,u1[0],f[6]);						// w2 = u0*f6
	_ff_add(w3,w2,w2);							// w3 = 2*u0*f6
	_ff_add(w5,w2,w3);							// w5 = 3*u0*f6
	_ff_mult(w4,w1,f[6]);						// w4 = u1^2*f6
	_ff_add(w2,w4,w4);							// w2 = 2u1^2*f6
	_ff_subfrom(w3,w2);
	_ff_subfrom(w3,w4);						// w3 = 2*u0*f6 - 3*u1^2*f6 (needed for g0 later)
	_ff_sub(t0,f[4],w5);
	_ff_addto(t0,w2);
	_ff_mult(t0,t0,u1[1]);
	_ff_x2(t0);
	_ff_sub(t1,f[3],t0);							// t1 = f3-2*u1*(f4-w5+w2)
	_ff_add(t0,w1,w1);
	_ff_sub(g0,u1[0],t0);
	_ff_neg(g1,u1[1]);
	_ff_add(s0,f[4],w4);
	_ff_neg(s1,u1[0]);
	_ff_x2(s1);
	_ff_sum_4_mults(t0,g0,g1,s0,s1,f[4],w1,f[3],w5);
	_ff_addto(t0,f[2]);
	_ff_subfrom(t0,w0);							// t0:=f2-w0+w5*(u0-2*w1)-f3*u1+w1*(f4+w4)-2*f4*u0;
//printf ("t' = %ld*x+%ld\n", _ff_get_ui(t1), _ff_get_ui(t0));

	// compute s' = s1*x+s0 = rs = t'*z' mod u   [5M]
	_ff_mult(w0,t0,z0);
	_ff_mult(w1,t1,z1);
	_ff_add_one(g0,u1[1]);
	ff_negate(g0);
	_ff_add(g1,t0,t1);
	_ff_addto(z0,z1);
	_ff_sum_2_mults(s1,z0,w1,g0,g1);
	_ff_subfrom(s1,w0);
	_ff_mult(t0,u1[0],w1);
	_ff_sub(s0,w0,t0);
//printf ("rs = %ld*x+%ld\n", _ff_get_ui(s1), _ff_get_ui(s0));

	// compute s=s1*x+s0, 1/r, and d=1/(f6-s1^2)  [6M+2S+1C+Inv]
	_ff_square(w4,r);
	_ff_square(w0,s1);
	_ff_mult(t0,w4,f[6]);
	_ff_subfrom(t0,w0);							// w1 = r^2*f6-s1^2
	_ff_mult (w2,r,t0);
	_ff_mult(r,r,w4);							// r now holds r^3
	if ( ctx && ! ctx->state ) {
		_ff_set (ctx->invert, w2);
		_ff_set (ctx->r, r);
		_ff_set (ctx->s0, s0);
		_ff_set (ctx->s1, s1);
		_ff_set (ctx->s2, t0);
		_ff_set(ctx->inv1, w3);
		ctx->state = 1;
		return 0;								// return to let caller invert
	}
	_ff_invert(w1,w2);							// w2 = w1/(r*w1)

hecurve_g2_square_d6_inverted:					// caller returns here after inverting	

	_ff_mult(d,r,w1);							// d = 1/(f6-s1^2)
	_ff_mult(r,w1,t0);							// r now holds 1/r
	_ff_mult(s0,r,s0);							// s0 = s0'/r
	_ff_mult(s1,r,s1);							// s1 = s1'/r
//printf ("d=%ld, 1/r=%ld s = %ld*x+%ld\n", _ff_get_ui(d), _ff_get_ui(r), _ff_get_ui(s1), _ff_get_ui(s0));

	// compute ell = s*u   [3M]
	_ff_mult(ell0,s0,u1[0]);						// ell0 = s0*u0
	_ff_mult(ell2,s1,u1[1]);
	_ff_add(w0,s0,s1);
	_ff_add(w1,u1[0],u1[1]);
	_ff_mult(ell1,w0,w1);
	_ff_subfrom(ell1,ell0);
	_ff_subfrom(ell1,ell2);						// ell1 = (s1+s0)*(u1+u0) - s0*u0 - s1*u1 (Karatsuba)
	_ff_addto(ell2,s0);							// ell2 = s1*u1+s0
//printf ("ell = %ld*x^3+%ld*x^2+%ld*x+%ld\n", _ff_get_ui(s1), _ff_get_ui(ell2), _ff_get_ui(ell1), _ff_get_ui(ell0));

	// compute g = g2*x^2+g1*x+g0 = (t-2vs)/u   [1M+1C]
	_ff_mult(w1,u1[1],f[6]);
	_ff_x2(w1);
	_ff_neg(g1,w1);							// g0 = -2u1*f6
	_ff_mult(w1,v1[1],s1);
	_ff_x2(w1);
	_ff_addto(w1,w3);
	_ff_sub(g0,f[4],w1);							// g0 = f4 - (2*v1*s1+2*u0*f6 - 3*u1^2*f6)
//printf ("g = -%ld*x^2+%ld*x+%ld\n", _ff_get_ui(f[6]), _ff_get_ui(g1), _ff_get_ui(g0));

	// compute u = x^2+u1*x+u0 = (g-s^2)/u made monic    [3M+1S]
	_ff_square(w0,s0);
	_ff_subfrom(g0,w0);
	_ff_mult(u[0],d,g0);							// u0 = d*(g0-s0^2)
	_ff_mult(w0,s0,s1);
	_ff_x2(w0);
	_ff_subfrom(g1,w0);
	_ff_mult(u[1],d,g1);							// u1 = d*(g1-2s0*s1)
	_ff_set_one(u[2]);
//printf ("u = x^2+%ld*x+%ld\n", _ff_get_ui(u[1]), _ff_get_ui(u[0]));

	// compute v = v1*x+v0 = -(ell+v) mod u    [3M]
	_ff_addto(ell1,v1[1]);						// ell1 now holds ell1+v1
	_ff_addto(ell0,v1[0]);						// ell0 now holds ell0+v0
	_ff_mult(w1,u[1],s1);						// w1 = u1*s1
	_ff_sub(w2,ell2,w1);							// w2 = ell2-w1
	_ff_mult(w3,w2,u[0]);						// w3 = u0*w2
	_ff_add(w0,u[0],u[1]);
	_ff_addto(s1,w2);
	_ff_mult(t0,w0,s1);
	_ff_subfrom(t0,ell1);
	_ff_subfrom(t0,w1);
	_ff_sub(v[1],t0,w3);							// v1 = -ell1+(u0+u1)*(s1+w2)-w1-w3
	_ff_sub(v[0],w3,ell0);
//printf ("v = %ld*x+%ld\n", _ff_get_ui(v[1]), _ff_get_ui(v[0]));

	// total operation count 26M+5S+4C

#if HECURVE_VERIFY
	if ( ! hecurve_verify (u, v, f) ) {
		err_printf ("hecurve_g2_square_d6: output verification failed for inputs:\n");
		err_printf ("      ");  hecurve_print(u1,v1);
		err_printf ("note that input has been modified if output overlapped.\n");
		exit (0);
	}
#endif
	return 1;
}

#endif
