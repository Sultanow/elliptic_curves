#ifndef _JACORDER_INCLUDE_
#define _JACORDER_INCLUDE_

/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <limits.h>
#include "jac.h"
#include "hcpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
	This module contains generic order computations for Jacobians using BSGS.
	The detail are described in [SutherlandThesis] and [KedlayaSutherland2007].
	
	Most of smalljac's life is spent in the function jac_parallel_search.
*/

#define JAC_INVALID_A1			LONG_MAX

struct a2tab_entry {
	double a1;
	double a2median;
	double a2mae;
};

int jac_order (unsigned long *pP1, unsigned long Min, unsigned long Max,  long a1, long d,  int fExponentOnly, int *constraints, hc_poly c[1]);
unsigned long jac_parallel_search (jac_t a[1],  unsigned long M, unsigned long W, unsigned long Min, unsigned long Max, int tiny, int parity, hc_poly c[1]);
int jac_search (mpz_t e[2], jac_t a[1],  unsigned long m, mpz_t Min, mpz_t Max, int repeats, hc_poly c[1]);

#ifdef __cplusplus
}
#endif

#endif
