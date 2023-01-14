#ifndef _SMALLJAC_G23_INCLUDE_
#define _SMALLJAC_G23_INCLUDE_

/*
Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "hcpoly.h"

#ifdef __cplusplus
extern "C" {
#endif

unsigned long smalljac_charpoly_gops;

int smalljac_genus2_charpoly (long a[3], long p, double sqrtp, long P1, long PN1);
int smalljac_genus3_charpoly (long a[3], long p, double sqrtp, long P1, long PN1, unsigned long pts);
int smalljac_genus2_charpoly_from_P1 (long a[2], long P1, long Min, long Max, hc_poly c[1]);
int smalljac_genus2_charpoly_from_Pmodp (long a[2], hc_poly c[1]);
int smalljac_genus3_charpoly_from_P1 (long o[3], long P1, long pts, long Min, long Max, hc_poly c[1]);
int smalljac_genus3_charpoly_from_Pmodp (long a[3], hc_poly c[1]);
int smalljac_genus2_charpoly_cmsquare (long a[], int D, hc_poly c[1]);
int smalljac_genus2_charpoly_from_a2_modp (long a[2], hc_poly c[1]);
int smalljac_genus2_charpoly_from_zero_Pmodp (long a[2], hc_poly c[1]);

#ifdef __cplusplus
}
#endif

#endif
