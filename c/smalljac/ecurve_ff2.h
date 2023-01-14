#ifndef _ECURVE_FF2_INCLUDE_
#define _ECURVE_FF2_INCLUDE_

/*
    Copyright (c) 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#include "ff_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

// returns the number of points on the curve y^2=x^3+f1*x+f0 over Fp^2 or zero if the curve is singular
// in characteristic 3, curve may be in the form y^2=x^3+f2*x^2+f1*x+f0, o.w. it only looks at f0 and f1 
long ecurve_ff2_group_order (ff_t f[8]);

#ifdef __cplusplus
}
#endif

#endif
