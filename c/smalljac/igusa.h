#ifndef _IGUSA_INCLUDE_
#define _IGUSA_INCLUDE_

/*
    Copyright (c) 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

int mpq_poly_igusa_inv (mpq_t I[3], mpz_t f[], int degree);			// returns 0 if curve is singular

// returns true if I is the Igusa triple (400000,-20000,-2000) for y^2=x^5+x, false otherwise
static inline int mpq_x5px_twist (mpq_t I[3])
{
	if ( mpq_sgn(I[0]) <= 0 || mpq_sgn(I[1]) >= 0 || mpq_sgn(I[2]) >= 0 ) return 0;
	return ( mpz_cmp_ui(mpq_numref(I[0]),400000) == 0 && mpz_cmp_ui(mpq_denref(I[0]),1) == 0  &&
		      mpz_cmpabs_ui(mpq_numref(I[1]),20000) == 0 && mpz_cmp_ui(mpq_denref(I[1]),1) == 0 &&
		      mpz_cmpabs_ui(mpq_numref(I[2]),2000) == 0 && mpz_cmp_ui(mpq_denref(I[2]),1) == 0 );
}

// returns true if I is the Igusa triple (51200000/3,480000,148000) for y^2=x^6+1, false otherwise
static inline int mpq_x6p1_twist (mpq_t I[3])
{
	return ( mpz_cmp_ui(mpq_numref(I[0]),51200000) == 0 && mpz_cmp_ui(mpq_denref(I[0]),3) == 0 &&
	              mpz_cmp_ui(mpq_numref(I[1]),480000) == 0 && mpz_cmp_ui(mpq_denref(I[1]),1) == 0 &&
		      mpz_cmp_ui(mpq_numref(I[2]),148000) == 0 && mpz_cmp_ui(mpq_denref(I[2]),1) == 0 );	
}

#ifdef __cplusplus
}
#endif

#endif
