#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "polyparse.h"
#include "mpzpolyutil.h"
#include "pointcount.h"
#include "smalljac_internal.h"

// note that we assume curves in weierstrass form are specified by a minimal model, we make no attempt to transform them

// for p=2,3 the tables below are indexed by n = a1+p*a2+p^2*a3+p^3*a4+p^4*a6 \in [0,p^5) with a_i reduced mod p.
char ws2goodtab[32] = {0,0,0,0,1,0,1,1,0,1,0,1,1,0,1,1,0,1,0,1,1,1,1,0,0,0,0,0,1,1,1,0,};
char ws2a1tab[32] = {0,-1,0,1,0,1,2,1,0,1,0,-1,2,1,0,1,0,1,0,-1,0,-1,-2,-1,0,-1,0,1,-2,-1,0,-1,};
char ws3goodtab[243] = {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,0,1,0,1,1,1,1,1,1,0,0,0,0,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,0,0,0,1,1,1,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,0,1,};
char ws3a1tab[243] = {0,-1,-1,-1,1,1,1,0,0,0,-1,2,2,1,1,1,3,0,0,2,-1,2,1,1,1,0,3,0,-1,-1,-1,1,1,1,0,0,0,2,-1,2,1,1,1,0,3,0,-1,2,2,1,1,1,3,0,0,2,2,2,-2,-2,-2,0,0,3,2,2,-1,1,1,1,0,0,3,2,2,-1,1,1,1,0,0,0,2,2,2,1,1,1,0,0,0,-1,-1,-1,1,-2,-2,-3,0,0,-1,-1,-1,-2,1,-2,0,-3,0,2,2,2,1,1,1,0,0,0,-1,-1,-1,-2,1,-2,0,-3,0,-1,-1,-1,1,-2,-2,-3,0,3,-1,-1,-1,1,1,1,3,3,-3,-1,-1,-1,-2,-2,1,0,0,-3,-1,-1,-1,-2,-2,1,0,0,0,-1,-1,-1,-2,-2,-2,0,0,0,2,-1,-1,-2,1,1,0,0,0,-1,2,-1,1,-2,1,0,0,0,-1,-1,-1,-2,-2,-2,0,0,0,-1,2,-1,1,-2,1,0,0,0,2,-1,-1,-2,1,1,0,0,-3,-1,-1,-1,1,1,1,-3,-3,0,-1,-1,2,1,1,-2,0,0,0,-1,-1,2,1,1,-2,0,0,};

// Computes either Lpoly coefficients or group structure for tiny values of p.  Returns -2 for error, -1 for bad reduction, number of entries in a[] otherwise (should never return 0)
// Curve must be defined over Q
int smalljac_tiny_Lpoly (long a[], smalljac_curve *sc, int p, unsigned long flags)
{
	unsigned long f[SMALLJAC_MAX_DEGREE+1];
	long P1,pts;

	assert ( sc->Qflag );
	if ( sc->special == SMALLJAC_SPECIAL_FK_TWIST || sc->special == SMALLJAC_SPECIAL_PICARD ) {
		if ( p == 2 ) return -1;	// TODO: this is a lie!
		return smalljac_internal_Lpoly_Q (a, sc, p, flags);
	}
	
	// handle special cases for Weierstrass equations at 2 and 3 first
	if ( p == 2 && (sc->flags&SMALLJAC_CURVE_FLAG_WS) ) {
		int ws2 = mpz_fdiv_ui (sc->H[1],2) + 2*mpz_fdiv_ui (sc->F[2],2) + 4*mpz_fdiv_ui (sc->H[0],2) + 8*mpz_fdiv_ui (sc->F[1],2) + 16*mpz_fdiv_ui (sc->F[0],2);
		a[0] = ws2a1tab[ws2];
		if ( (flags&SMALLJAC_GROUP) ) a[0] += 2 + ws2goodtab[ws2];		// Over F_2 we always have a cyclic group
		return ws2goodtab[ws2] ? 1 : -1;
	}
	if ( p == 3 && (sc->flags&SMALLJAC_CURVE_FLAG_WS) ) {
		int ws3 = mpz_fdiv_ui (sc->H[1],3) + 3*mpz_fdiv_ui (sc->F[2],3) + 9*mpz_fdiv_ui (sc->H[0],3) + 27*mpz_fdiv_ui (sc->F[1],3) + 81*mpz_fdiv_ui (sc->F[0],3);
		if ( (flags&SMALLJAC_GROUP) ) {
			// there are exactly 9 cases where the group is not cyclic (and this necessarily means good reduction), and all of them are iso to Z/2 x Z/2
			switch (ws3) {
			case  54: // [0,0,0,2,0]
			case  61: // [1,2,0,2,0]
			case  62: // [2,2,0,2,0]
			case 178: // [1,2,1,0,2]
			case 188: // [2,2,2,0,2]
			case 206: // [2,2,1,1,2]
			case 214: // [1,2,2,1,2]
			case 225: // [0,0,1,2,2]
			case 234: // [0,0,2,2,2]
			a[0] = a[1] = 2;  return 2;
			}
		}
		a[0] = ws3a1tab[ws3];
		if ( (flags&SMALLJAC_GROUP) ) a[0] += 3 + ws3goodtab[ws3];
		return ws3goodtab[ws3] ? 1 : -1;
	}
	
	// We assume that we have good reduction if we get here, sc->D should be checked by caller
	
	if ( p == 2 ) {
		ff2k_t f[2*SMALLJAC_MAX_GENUS+3], h[SMALLJAC_MAX_GENUS+2];
		int df, dh;
		
		ff2k_setup (1);
		df = ff2k_poly_set_mpz (f, sc->F, sc->dF);
		dh = ff2k_poly_set_mpz (h, sc->H, sc->dH);
		sc->pts = ff2k_hyperelliptic_pointcount (f, df, h, dh);
		a[0] = sc->pts-p-1;
		if ( sc->genus == 1 || (flags&SMALLJAC_A1_ONLY) ) return 1;
		ff2k_setup (2);	// note thate we don't need to reset f and h, F_2 elements automatically lift
		pts = ff2k_hyperelliptic_pointcount (f, df, h, dh);
		a[1] = (pts - p*p - 1 + a[0]*a[0]) / 2;
		if ( sc->genus == 2 ) return 2;
		ff2k_setup (3);	// note thate we don't need to reset f and h, F_2 elements automatically lift
		pts = ff2k_hyperelliptic_pointcount (f, df, h, dh);
		a[2] = (pts - (long)p*p*p - 1 - a[0]*a[0]*a[0] + 3*a[0]*a[1]) / 3;
		if ( sc->genus == 3 ) return 3;
		err_printf ("Unhandled genus %d in smalljac_Lpoly_tiny for p=%d\n", sc->genus, p);
		return -2;
	}

	ui_poly_set_mpz_mod_p (f, sc->f, sc->degree, p);
	if ( sc->special == SMALLJAC_SPECIAL_PICARD ) {
		sc->pts = pointcount_tiny_pd4 (f, sc->degree, p);
	} else {
		sc->pts = pointcount_tiny (f, sc->degree, p);
	}
	
	if ( sc->genus == 1 ) {
		if ( (flags&SMALLJAC_GROUP) ) {
			ff_setup_ui (p);
			hc_poly_set_mpz (sc->hc, sc->f, sc->degree);
			return ecurve_group_structure (a,sc->pts,0,sc->hc->f);					// we don't have any torsion info so specify d=0
		}
		a[0] = sc->pts-p-1;
		return 1;
	}
	a[0] = sc->pts-p-1;
	if ( (flags&SMALLJAC_A1_ONLY) ) return 1;
	pts = pointcount_tiny_d2 (f, sc->degree, p);
	a[1] = (pts - p*p - 1+a[0]*a[0]) / 2;
	if ( sc->genus == 2 ) {
		if ( (flags&SMALLJAC_GROUP) ) {
			P1 = smalljac_Lp1_ui (a, sc->genus, p);
			ff_setup_ui (p);
			hc_poly_set_mpz (sc->hc, sc->f, sc->degree);
			return jac_structure (a, sc->hc, P1, 0);
		} else {
			return 2;
		}
	}

	pts = pointcount_tiny_d3 (f, sc->degree, p);
	a[2] = (pts - (long)p*p*p - 1 - a[0]*a[0]*a[0] + 3*a[0]*a[1]) / 3;
	if ( sc->genus == 3 ) {
		if ( (flags&SMALLJAC_GROUP) ) {
			P1 = smalljac_Lp1_ui (a, sc->genus, p);
			ff_setup_ui (p);
			hc_poly_set_mpz (sc->hc, sc->f, sc->degree);
			return jac_structure (a, sc->hc, P1, flags&SMALLJAC_SGROUP);
		} else {
			return 3;
		}
	}
	err_printf ("Unhandled genus %d in smalljac_Lpoly_tiny for p=%d\n", sc->genus, p);
	return -2;
}

