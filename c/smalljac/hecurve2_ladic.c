#include <stdlib.h>
#include <stdio.h>
#include "ff_poly.h"
#include "g2tor3poly.h"
#include "hecurve.h"

/*
    Copyright 2009-2012 Andrew V. Sutherland

    This file is part of smalljac.

    smalljac is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    smalljac is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with smalljac.  If not, see <http://www.gnu.org/licenses/>.
*/

#define MAX_POWER			19

struct poly_4coeff *g2tor3[41] = { g2tor3d0, g2tor3d1, g2tor3d2, g2tor3d3, g2tor3d4, g2tor3d5, g2tor3d6, g2tor3d7, g2tor3d8, g2tor3d9,
g2tor3d10, g2tor3d11, g2tor3d12, g2tor3d13, g2tor3d14, g2tor3d15, g2tor3d16, g2tor3d17, g2tor3d18, g2tor3d19,	
g2tor3d20, g2tor3d21, g2tor3d22, g2tor3d23, g2tor3d24, g2tor3d25, g2tor3d26, g2tor3d27, g2tor3d28, g2tor3d29,
g2tor3d30, g2tor3d31, g2tor3d32, g2tor3d33, g2tor3d34, g2tor3d35, g2tor3d36, g2tor3d37, g2tor3d38, g2tor3d39,g2tor3d40};

// IMPORTANT: f is assumed to be of the form x^5+f3x^3+f2x^2+f1x+f0
void hecurve_g2_tor3_modpoly (ff_t g[POLY_G2TOR3_DEGREE+1], ff_t f[6])
{
	register int i, j, k;
	ff_t powers[4][20];
	ff_t t0,c;
	
	// compute coefficient powers
	for ( i = 0 ; i < 4 ; i++ ) {
		_ff_set_one(powers[i][0]);
		_ff_set (powers[i][1], f[i]);
		for ( j = 2 ; j <= MAX_POWER ; j++ ) _ff_mult(powers[i][j],powers[i][j-1],f[i]);
	}
	
	/*
		The 3-torsion modular poly provided by Gaudry & Schost assumes f(x) = x^5+f0x^3+f1x^2+f2x + f3, while we use x^5+f3x^3+f2x^2+f1x+f0.
		Thus the need for replacing k with 3-k below.
	*/
	for ( i = 0 ; i <= POLY_G2TOR3_DEGREE ; i++ ) {
		_ff_set_zero(c);
		for ( j = 0 ; g2tor3[i][j].c ; j++ ) {
			_ff_set_i (t0,g2tor3[i][j].c);
			for ( k = 0 ; k < 4 ; k++ ) ff_mult(t0,t0,powers[3-k][g2tor3[i][j].f[k]]);
			_ff_addto(c,t0);
		}
		_ff_set(g[i],c);
	}
}

int hecurve_g2_3tor(ff_t f[6])
{
	ff_t g[POLY_G2TOR3_DEGREE+1], h[POLY_G2TOR3_DEGREE], r[POLY_G2TOR3_DEGREE];
	int k, d_h, d_r;
	
	hecurve_g2_tor3_modpoly(g,f);
	
	// first make sure modular poly is square free
	ff_poly_derivative (h,&d_h,g,POLY_G2TOR3_DEGREE);
	ff_poly_gcd (r,&d_r,g,POLY_G2TOR3_DEGREE,h,d_h);
	if ( d_r > 0 ) return 0;
	k = ff_poly_count_roots(g,POLY_G2TOR3_DEGREE);
	if ( ! k ) return 1;
	if ( k ==2 ||k==5 ||k==8 ) return 3;
	return 0;
}
