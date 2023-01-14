#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include "mpzutil.h"
#include "ff_poly.h"
#include "ecurve_ff2.h"

/*
    Copyright (c) 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#define BSGS_MAX_STEPS          4096        // handles p up to 2^22, which means #Fp^2  up to 2^44
#define FF2_ECURVE_VERIFY       0           // nonzero to force point verification, used for debugging

/*
    BSGS point-counting for elliptic curves in short Weierstarss form y^2=x^3+f1*x+f0 over Fp^2
    We must handle the CM case f0=0 separately, since we use (0,0) to represent the identity in affine coords.
*/

#if FF2_ECURVE_VERIFY
static ff_t _ecurve_ff2_f[4];       // cached curve coefficients, set by ecurve_ff2_random_point, used by ecurve_ff2_verify_point
static inline void ecurve_ff2_verify_point (ff_t x[2], ff_t y[2], char *errstr)
{
    ff_t w[2], z[2];
    
    if  ( ff2_zero(x) && ff2_zero(y) ) return;
    ff2_square (z, y);  ff2_square (w, x);  ff2_add (w, w, _ecurve_ff2_f+2);  ff2_mult (w, w, x);  ff2_add (w, w, _ecurve_ff2_f);
    if ( ff2_equal (w,z) ) return;
    printf ("ecurve_ff2 point verification failed with p=%ld, z^2=%ld:\n    the point (%ld*z+%ld,%ld*z+%ld) is not on the curve y^2 = x^3 + (%ld*z+%ld)*x +(%ld*z+%ld)\n",
           _ff_p, _ff_get_ui (_ff_2g), _ff_get_ui (x[1]), _ff_get_ui (x[0]), _ff_get_ui (y[1]), _ff_get_ui (y[0]), _ff_get_ui(_ecurve_ff2_f[3]), _ff_get_ui(_ecurve_ff2_f[2]), _ff_get_ui(_ecurve_ff2_f[1]), _ff_get_ui(_ecurve_ff2_f[0]));
    if ( errstr ) printf ("    %s\n", errstr);
    abort();
}
#else
// this should compile away into nothingness...
static inline void ecurve_ff2_verify_point (ff_t x[2], ff_t y[2], char *errstr) {}
#endif


static inline int ecurve_ff2_id (ff_t x[2], ff_t y[2]) { return ff2_zero(x) && ff2_zero(y); }

static inline void ecurve_ff2_invert (ff_t x[2], ff_t y[2])
    { ff2_negate (y); }
    
static inline void ecurve_ff2_print_point (ff_t x[2], ff_t y[2])
    { printf ("(%ld*z+%ld,%ld*z+%ld)", _ff_get_ui(x[1]), _ff_get_ui(x[0]), _ff_get_ui(y[1]), _ff_get_ui(y[0])); }

// compute disc(x^3+f1*x+f0) = -4*f1^3 - 27*f0^2 in Fp^2
int ff2_poly_x3axb_disc (ff_t d[2], ff_t f[4])
{
    ff_t t1[2], t2[2], t;
    
    ff2_square (t1, f+2);  ff2_mult (t1, t1, f+2);  ff2_add (t1, t1, t1);  ff2_add (t1, t1, t1);                // t1 = 4*f1^3
    ff2_square (t2, f); _ff_set_ui (t, 27); ff2_scalar_mult (t2, t, t2);  ff2_add (d, t1, t2);                  // t2 = 27*f0^2
    ff2_negate (d);
    return ! ff2_zero (d);
}

// compute disc(x^3+f2*x^2+f1*x+f0) = f2^2*(f1^2-f0*f2) - f1^3 in F_9
int ff2_poly_cubic_disc_char3 (ff_t d[2], ff_t f[6])
    { ff_t t1[2], t2[2];  ff2_square (t1, f+2); ff2_mult (t2, f, f+4); ff2_sub (d, t1, t2);  ff2_square (t2, f+4); ff2_mult (d, d, t2); ff2_mult (t1, t1, f+2); ff2_sub (d, d, t1); return ! ff2_zero(d); }
    
// generates a random point on the curve y^2=x^3+f1*x+f0
void ecurve_ff2_random_point (ff_t x[2], ff_t y[2], ff_t f[4])
{
#if FF2_ECURVE_VERIFY
    ff2_set (_ecurve_ff2_f,f);  ff2_set (_ecurve_ff2_f+2,f+2);  // cache f0, and f1 for later use by ecurve_ff2_verify_point
#endif
    do {
        ff2_random(x);  ff2_square (y, x);  ff2_add (y, y, f+2);  ff2_mult (y, y, x);  ff2_add (y, y, f);
    } while ( ! ff2_sqrt (y, y) );
    ecurve_ff2_verify_point (x, y, "random_point output");
}

// computes (x3,y3) = 2*(x1,y1)
void ecurve_ff2_dbl (ff_t x3[2], ff_t y3[2], ff_t x1[2], ff_t y1[2], ff_t f1[2])
{
    ff_t t1[2], t2[2], t3[2];

    ecurve_ff2_verify_point (x1, y1, "dbl input");
    if ( ff2_zero(y1) ) { ff2_set_zero (x3);  ff2_set_zero (y3); return; }                  // doubling a 2-torsion point yields the identity
    ff2_square (t1, x1);  ff2_add (t2, t1, t1);  ff2_add (t2, t2, t1);  ff2_add (t2, t2, f1);       // t2 = 3x1^2+f1
    ff2_add (t1, y1, y1);  ff2_invert (t1, t1);
    ff2_mult (t1, t1, t2);                                                  // t1 = lambda = (3x1^2+f1)/(2y1)
    ff2_square (t3, t1);  ff2_sub (t3, t3, x1);  ff2_sub (t3, t3, x1);                      // t3 = x3 = lambda^2 - 2x1
    ff2_sub (t2, x1, t3);                                                       // t2 = x1-x3
    ff2_set (x3, t3);                                                       // x3 = t3
    ff2_mult (t3, t1, t2);  ff2_sub (y3, t3, y1);                                   // y3 = (x1-x3)lambda - y1
    ecurve_ff2_verify_point (x3, y3, "dbl output");
    // I+2M+2S+7A}
}  

// computes (x3,y3) = (x1,y1) + (x2,y2) on the curve y^2=x^3+f1*x+f0 over Fp^2.  Handles all cases, any combo of points may coincide.
void ecurve_ff2_add (ff_t x3[2], ff_t y3[], ff_t x1[2], ff_t y1[2], ff_t x2[2], ff_t y2[2], ff_t f1[2])
{
    ff_t t1[2], t2[2], t3[2];
    
    ecurve_ff2_verify_point (x1, y1, "add input 1");
    ecurve_ff2_verify_point (x2, y2, "add input 2");
    if ( ff2_zero(x1) && ff2_zero(y1) ) { ff2_set (x3, x2); ff2_set (y3, y2); return; }     // (x1,y1)=(0,0) is the identity
    if ( ff2_zero(x2) && ff2_zero(y2) ) { ff2_set (x3, x1); ff2_set (y3, y1); return; }     // (x2,y2)=(0,0) is the identity
    if ( ff2_equal (x1, x2) ) {
        if ( ff2_equal (y1, y2) ) { ecurve_ff2_dbl (x3, y3, x1, y1, f1); return; }          // double if points are equal
        ff2_set_zero (x3); ff2_set_zero (y3);                                       // points must be opposite, return the identity
        return;
    }
    ff2_sub (t1, x2, x1);  ff2_invert (t1, t1);  ff2_sub (t2, y2, y1);  ff2_mult (t1, t1, t2);      // t1 = lambda = (y2-y1)/(x2-x1)
    ff2_square (t3, t1);  ff2_sub (t3, t3, x1);  ff2_sub (t3, t3, x2);                      // t3 = x3 = lambda^2 - x1 - x2 (don't write x3 yet)
    ff2_sub (t2, x1, t3);  ff2_mult (t2, t2, t1);  ff2_sub (y3, t2, y1);                    // y3 = (x1-x3)lambda - y1
    ff2_set (x3,t3);                                                            // x3 = t3
    ecurve_ff2_verify_point (x3, y3, "add output");
    // I+3M+2S+6A
}

// standard 4-ary exponentiation (fixed 2-bit window)
void ecurve_ff2_scalar_mult (ff_t x3[2], ff_t y3[2], ff_t x1[2], ff_t y1[2], long e, ff_t f1[2])
{
    register int i;
    register unsigned long j, m;
    struct { ff_t x[2], y[2]; } b[4];
    ff_t x[2], y[2];
    int s;
    
    ecurve_ff2_verify_point (x1, y1, "scalar mult input");
    switch (e) {
    case -2: ecurve_ff2_dbl (x3, y3, x1, y1, f1);  ff2_negate (y3); return;
    case -1: ff2_set (x3, x1);  ff2_set (y3, y1); ff2_negate (y3); return;
    case 0:  ff2_set_zero (x3);  ff2_set_zero(y3);  return;
    case 1:  ff2_set (x3, x1);  ff2_set (y3, y1);  return;
    case 2:  ecurve_ff2_dbl (x3, y3, x1, y1, f1);  return;
    }
    if ( e < 0 ) { e = -e; s =-1; } else s = 1;
    
    i = _asm_highbit(e);
    if ( i&1 ) i--;
    m = 3UL<<i;
    ff2_set (b[1].x, x1);  ff2_set (b[1].y, y1);
    ecurve_ff2_dbl (b[2].x, b[2].y, x1, y1, f1);
    ecurve_ff2_add (b[3].x, b[3].y, b[2].x, b[2].y, x1, y1, f1);
    j = ((m&e)>>i);     // notej cannot be 0
    ff2_set (x, b[j].x);  ff2_set (y, b[j].y);
    for ( m>>=2,i-=2 ; m ; m>>=2, i-=2 ) {
        ecurve_ff2_dbl (x, y, x, y, f1);  ecurve_ff2_dbl (x, y, x, y, f1);
        j = (m&e)>>i;
        if ( j ) ecurve_ff2_add (x, y, x, y, b[j].x, b[j].y, f1);
    }
    ff2_set (x3, x);  ff2_set (y3, y);
    if ( s < 0 ) ff2_negate (y3);
    ecurve_ff2_verify_point (x1, y1, "scalar mult output");
}

/*
    Compute the order of affine point (x,y) given that e*(x,y)=1, using the classical algorithm.

    This algorithm is not particularly fast, but it doesn't need to be, as it is not invoked for most p.
    The average number of bits exponentiated per prime for a non-CM curve (fixed curve, varying p) is less than 1.

    Copied from hecurve1x.c
*/
long ecurve_ff2_fastorder (ff_t x[2], ff_t y[2], long e, ff_t f1[2])
{
    struct { ff_t x[2], y[2]; } b[MAX_UI_PP_FACTORS];
    unsigned long q, p[MAX_UI_PP_FACTORS], h[MAX_UI_PP_FACTORS];
    register long o;
    register int i, j, k, w;

    ecurve_ff2_verify_point (x, y, "fastorder input");
    if ( e < 0 ) e = -e;
    w = ui_factor(p,h,e);
    if ( w == 1 && h[0] == 1 ) return e;    // not hard when the exponent is prime
        
    for ( i = 0 ; i < w ; i++ ) {
        for ( q=p[i],j=1 ; j < h[i] ; j++ ) q *= p[i]; 
        o = e/q;
        ecurve_ff2_scalar_mult (b[i].x, b[i].y, x, y, o, f1);
    }
    o = 1;
    for ( i = 0 ; i < w ; i++ ) {
        for ( j = 0 ; j < h[i]-1 ; j++ ) {
            if ( ecurve_ff2_id (b[i].x, b[i].y) ) break;
            ecurve_ff2_scalar_mult (b[i].x, b[i].y, b[i].x, b[i].y, p[i], f1);
        }
        if ( ! ecurve_ff2_id (b[i].x, b[i].y) ) j++;
        for ( k = 0 ; k < j ; k++ ) o *= p[i];
    }
//printf ("ecurve_ff2_fastorder computed order of (%ld*z+%ld,%ld*z+%ld) is %ld from exponent %ld over F_%ld^2 (z^2=%ld)\n",
//     _ff_get_ui (x[1]), _ff_get_ui(x[0]), _ff_get_ui(y[1]), _ff_get_ui(y[0]),o, e, _ff_p, _ff_get_ui(_ff_2g));
    return o;
}


#define BSGS_TABSIZE            BSGS_MAX_STEPS          // don't make this too big, it takes time to initialize it.  A few collisions won't kill us.
#define BSGS_TABMASK            ((unsigned long)(BSGS_TABSIZE-1))

static struct {
    ff_t x[2];
    ff_t y[2];
} babys[BSGS_MAX_STEPS];

// this code is copied from hecurve1x.c
static unsigned short hashtab[BSGS_TABSIZE];
static struct tab_entry {
    ff_t x[2];
    short i;
    short next;
} entries[BSGS_MAX_STEPS+1];
short nexttabentry;

static inline void tab_clear() { memset(hashtab,0,sizeof(hashtab)); nexttabentry = 1; } // don't use entry 0

// we require inserts to have unique x values, return -1 if unique, otherwise return index value for existing entry
static inline int tab_insert(ff_t x[2], short i)
{
    register struct tab_entry *p;
    register int n,h;

    h = (x[0]^x[1])&BSGS_TABMASK;                   // brutally simple hash function, probably not worth making this more sophisticated
    n = hashtab[h];
    if ( n ) {
        p=entries+n;
        for(;;) {
            if ( ff2_equal(p->x,x) ) return p->i;
            if ( ! p->next ) break;
            p=entries+p->next;
        }
    }
    p = entries+nexttabentry;
    ff2_set(p->x,x);
    p->i = i;
    p->next = n;
    hashtab[h] = nexttabentry++;
    return -1;
}

static inline int tab_lookup(ff_t x[2])
{
    register struct tab_entry *p;
    register int n;

    n = hashtab[(x[0]^x[1])&BSGS_TABMASK];          // brutally simple hash function, probably not worth making this more sophisticated
    if ( ! n ) return -1;
    p=entries+n;
    for(;;) {
        if ( ff2_equal(p->x,x) ) return p->i;
        if ( ! p->next ) break;
        p=entries+p->next;
    }
    return -1;
}

// given P=(x,y) of order n in [low, high], computes either e=n (returning 0), e=unique multiple of n in [low,high] (returning 1)
int ecurve_ff2_bsgs (long *e, ff_t x[2], ff_t y[2], long low, long high, ff_t f1[2])
{
    long bsteps, gstepsize, gstepcenter, dstep, ustep, e1, e2;
    ff_t dx[2], dy[2], ux[2], uy[2], gx[2], gy[2], ngy[2];
    int i, j;

    ecurve_ff2_verify_point (x, y, "bsgs input");

//printf ("BSGS input point P="); ecurve_ff2_print_point (x,y); puts ("");
    
    // deal with 2-torsion now so we don't have to worry about it later
    if ( ff2_zero(y) ) { *e = ( ff2_zero(x) ? 1 : 2 ); return 0; }
    
    // by matching inverses, we expect to need half as many giant steps,
    // also, we simultaneously search up/down and expect to abort one of these half-way through
    // thus we should scale bsteps by 3/4*1/2=3/8
    bsteps = (long) sqrt((3*(high-low+1))/8);
    if ( bsteps >= BSGS_MAX_STEPS ) { printf ("Exceeded BSGS_MAX_STEPS = %d\n", BSGS_MAX_STEPS);  abort(); }
    tab_clear ();
    ff2_set (babys[1].x, x);  ff2_set (babys[1].y, y);  tab_insert (x, 1);
    for ( i = 2 ; i <= bsteps ; i++ ) {
        ecurve_ff2_add (babys[i].x, babys[i].y, babys[i-1].x, babys[i-1].y, x, y, f1);
        j = tab_insert (babys[i].x, i);
        // check for collision -- if this occurs it must be a collision with (x,-y) (note (x,y) is not 2-torsion)
        if ( j > 0 ) {
            ff2_add (babys[i].y, babys[i].y, babys[j].y);
            if ( ! ff2_zero (babys[i].y) ) { printf ("Baby-step collision between %d and %d didn't hit inverses.  This is impossible.\n", i, j); abort(); }
            *e = i+j;
            return 0;
        }
    }
    ecurve_ff2_dbl (gx, gy, babys[bsteps].x, babys[bsteps].y, f1);
    if ( ecurve_ff2_id (gx, gy) ) { *e = 2*bsteps;  return 0; }             // 2*bsteps must be the exact order, ow it would have a factor <= bsteps
    ecurve_ff2_add (gx, gy, gx, gy, x, y, f1);
    if ( ecurve_ff2_id (gx, gy) ) { *e = 2*bsteps+1;  return 0; }               // 2*bsteps+1 must be the exact order, ow it would have a factor <= bsteps
    ff2_neg (ngy, gy);
    gstepsize = 2*bsteps+1;
    gstepcenter = low + (high-low)/2;
    e1 = e2 = 0;
    ecurve_ff2_scalar_mult (ux, uy, x, y, gstepcenter, f1);
    if ( (j = tab_lookup (ux)) >= 0 ) e1 = ( ff2_equal (uy, babys[j].y) ? gstepcenter - j : gstepcenter + j );
    if ( ecurve_ff2_id (ux, uy) ) e1 = gstepcenter; 
    dstep = ustep = gstepcenter;
    ff2_set (dx, ux);  ff2_set (dy, uy);
    while ( dstep > low+bsteps || ustep < high-bsteps ) {
        if ( dstep > low+bsteps ) {
            ecurve_ff2_add (dx, dy, dx, dy, gx, ngy, f1);  dstep -= gstepsize;
            if ( ecurve_ff2_id (dx, dy) ) { e2 = dstep; if ( ! e1 ) { e1 = e2; e2 = 0; } else break; }
            else if ( (j = tab_lookup (dx)) >= 0 ) {
                e2 = ( ff2_equal (dy, babys[j].y) ? dstep - j : dstep + j );
                if ( e2 < low || e2 > high ) { e2 = 0; continue; }
                if ( ! e1 ) { e1 = e2; e2 = 0; } else break;
                dstep = low;                                        // done with downward steps
            }
        }
        if ( ustep < high-bsteps ) {
            ecurve_ff2_add (ux, uy, ux, uy, gx, gy, f1);  ustep += gstepsize;
            if ( ecurve_ff2_id (ux, uy) ) { e2 = ustep; if ( ! e1 ) { e1 = e2; e2 = 0; } else break; }
            else if ( (j = tab_lookup (ux)) >= 0 ) {
                e2 = ( ff2_equal (uy, babys[j].y) ? ustep - j : ustep + j );
                if ( e2 < low || e2 > high ) { e2 = 0; continue; }
                if ( ! e1 ) { e1 = e2; e2 = 0; } else break;
                ustep = high;                                           // done with upward steps
            }
        }
    }
    
//printf ("BSGS found exponents %ld and %ld\n", e1, e2);

    if ( ! e1 ) { printf ("BSGS search failed to find any exponent of the point (%ld*z+%ld,%ld*z+%ld) in the interval [%ld,%ld] over F_%ld^2 (z^2=%ld)\n",
                             _ff_get_ui (x[1]), _ff_get_ui(x[0]), _ff_get_ui(y[1]), _ff_get_ui(y[0]), low, high, _ff_p, _ff_get_ui(_ff_2g));  abort(); }
    if ( ! e2 ) { *e = e1;  return 1; }
    *e = ecurve_ff2_fastorder (x, y, e2-e1, f1);
    return 0;
}

// Algorithm to compute the trace of E/Fp^2 defined by y^2 = x^3 + ax (which has j-invariant 1728)
long ecurve_ff2_1728_trace (ff_t a[2])
{
    ff_t b[2], c[2], x, y, w;
    long u, v, s, p;
    
    if ( ff2_zero(a) ) { printf ("a must be nonzero in ecurve_ff2_1728_trace\n"); abort(); }
    
    p = (long)_ff_p;
    
    // handle supersingular case (p=3 mod 4) first
    if ( (p&3) == 3 ) {
        // compute b=a^((p^2-1)/4)
        ff2_bar (b, a);                         // b = a^p
        ff2_invert (c, a);
        ff2_mult (b, b, c);                     // b = a^(p-1)
        ff2_exp_ui (b, b, (p+1)/4);             // b = a^(p^2-1)/4)
        if ( ff2_pm_one (b) ) {
            u = ( ff2_one(b) ? -p : p );
        } else {
            u = 0;
        }
        return 2*u;
    }
    
    // compute a solution (u,v) to p=x^2+y^2
    if ( ! ff_cornacchia (&u, &v, 1) ) { printf ("cornacchia failed, unable to write p=%ld = 1 mod 4 as the sum of two squares ?!!\n", _ff_p); abort(); }
    // use it  to obtain a solution (u,v) to p^2=x^2+y^2
    v = 2*u*v;  u = p-2*u*u;
    if ( v&1 ) { s = u; u = v; v = s; }             // make u odd
    if ( (u-1-v)&3 ) u = -u;                        // make u-1=v mod 4
    
    // The algorithm below is based on reversing Algorithm 3.4 of Rubin-Silverberg "Choosing the correct elliptic curve in the CM method"
    // note that they write the curve as y^2 = x^3 - ax, but it makes no difference becase (-1)^((p^2-1)/4) = 1
    
    // compute w=(-a)^((p^2-1)/4), which is a 4th root of unity and must lie in Fp, since p=1 mod 4
    ff2_norm (&w, a);                       // w = a^(p+1)
    ff_exp_ui (&w, &w, (p-1)/4);                // w = a^((p^2-1)/4) = (-a)^((p^2-1)/4)
    
    // w is one of +/-1,+/-u/v, handle each of the 4 possibilities
    if ( _ff_one (w) ) return 2*u;              // w=1 =>trace is 2u
    if ( _ff_neg_one(w) ) return -2*u;          // w=-1 => trace is -2u
    _ff_set_i (x, u);  _ff_set_i (y, v);
    _ff_invert (y, y); _ff_mult (x, x, y);          // compute u/v in Fp
    if ( _ff_equal(w, x) ) return -2*v;         // w=u/v => trace is -2v
    ff_negate(w);
    if ( _ff_equal(w, x) ) return 2*v;          // w=-u/v => trace is 2v
    
    err_printf ("fell through to unreachable code in ecurve_ff2_1728_trace\n"); abort();
}


#define FF2_ECURVE_COUNT_P      7           // By Theorem 2 of CS2010, for p^2 > 49 (p > 7) we can use ecurve_ff2_group_order

static inline int ff2_idx(ff_t x[2]) { return _ff_get_ui(x[1])*_ff_p+_ff_get_ui(x[0]); }

// naive brute force point couting (this is very slow code that only gets used for p=3,5,7)
long ecurve_ff2_count_points (ff_t f[8])
{
    char sqrs[2*FF2_ECURVE_COUNT_P*FF2_ECURVE_COUNT_P];
    ff_t x[2], y[2];
    long pts;
    
    if ( _ff_p > FF2_ECURVE_COUNT_P ) { err_printf ("p must be <= %d to in ecurve_ff2_count_points\n", FF2_ECURVE_COUNT_P); abort(); }

    memset (sqrs, 0, sizeof(sqrs));
    ff2_set_zero (x);
    do {
        ff2_square(y,x);
        sqrs[ff2_idx(y)] = 1;
        _ff_inc (x[0]);  if ( _ff_zero (x[0]) ) _ff_inc(x[1]);
    } while ( ! ff2_zero (x) );
    pts = 0;
    ff2_set_zero (x);
    if ( _ff_p == 3 ) {
        do {
            ff2_add (y, x, f+4);  ff2_mult (y, y, x); ff2_add (y, y, f+2); ff2_mult (y, y, x); ff2_add (y, y, f);
            if ( ff2_zero (y) ) pts++;
            else if ( sqrs[ff2_idx(y)] ) pts += 2;
            _ff_inc (x[0]);  if ( _ff_zero (x[0]) ) _ff_inc(x[1]);
        } while ( ! ff2_zero (x) );
    } else {
        do {
            ff2_square (y, x); ff2_add (y, y, f+2); ff2_mult (y, y, x); ff2_add (y, y, f);
            if ( ff2_zero (y) ) pts++;
            else if ( sqrs[ff2_idx(y)] ) pts += 2;
            _ff_inc (x[0]);  if ( _ff_zero (x[0]) ) _ff_inc(x[1]);
        } while ( ! ff2_zero (x) );
    }
    return pts+1;
}

/*
    Las Vegas algorithm to compute the group order of E/Fp^2 defined by y^2 =f(x) where f is a monic cubic (f2 need not be zero)
    Implements the algorithm described in Cremona-Sutherland "On a Theorem of Mestre and Schoof"

    The coefficients f0, f1, and f2 are stored at f, f+2, f+4

    returns -1 if the curve is singular.
*/
long ecurve_ff2_group_order (ff_t F[8])
{
    ff_t x[2], y[2], g[4], *f, *h, _f[4];
    long a, m, a1, m1, e, low, high, p, q;
    int i, tor2, unique;
    
// printf ("f2_ecurve_group_order: "); ff2_poly_print(F,3);

    assert (ff2_one(F+6));
    
    // handle characteristic 3 separately
    if ( _ff_p == 3 ) {  ff2_poly_cubic_disc_char3 (x, F); return ( ff2_zero(x) ? -1 : ecurve_ff2_count_points(F) ); }
    
    // depress F if needed and set f to point to the appropiate poly (F or _f)
    if ( ! ff2_zero (F+4) ) { ff2_poly_depress_monic_cubic (_f, F); f = _f; } else { f = F; }
    
    // We need to handle the special case that f0=0 separately, since we assume we don't want (0,0) as an affine point
    // But in this case j=1728 and counting points is very easy
    if ( ff2_zero (f) ) { if ( ff2_zero(f+2) ) return -1; else return (long)(_ff_p*_ff_p)+1-ecurve_ff2_1728_trace (f+2); }

    // check that curve is non-singular (and get the discriminant of f(x) which we use for an easy 2-torsion check)
    ff2_poly_x3axb_disc (x, f);
    if ( ff2_zero(x) ) return -1;
    tor2 = ( ff2_sqrt (x, x) ? 0 : 1 );         // if  the discriminant is not a square then we must have 2-torsion (but otherwise don't assume anything)
    
    // note we require p^2 > 49 in order to apply Theorem 2 of Cremona-Sutherland, revert to naive point counting for p <= 7
    if ( _ff_p <= FF2_ECURVE_COUNT_P ) return ecurve_ff2_count_points (f);
    if ( _ff_p > (1UL<<31) ) { printf ("p=%ld is too large to count points in Fp^2 using ecurve_ff2_group_order\n", _ff_p); abort(); }
    p = _ff_p;                                                          // put p in a (signed) long for convenience
    q = p*p;
    
    low = q+1-2*p;  high = q+1+2*p;                                     // compute Hasse bounds
    h=f;                                                                // h will alternate between f and g (twist coefficients)
    a = 0;  m = 1;
    if ( tor2 ) { low = (low+1)/2;  high = high/2;  m = 2; }                        // optimize for (easily detected) 2-torsion but nothing else
    for ( i = 0 ; ; i++ ) {
        ecurve_ff2_random_point (x, y, h);
        if ( tor2 ) ecurve_ff2_dbl (x, y, x, y, h+2);                           // note that tor2 applies to both the curve and its twist
        unique = ecurve_ff2_bsgs (&e, x, y, low, high, h+2);
        if ( tor2 ) e *= 2;
        if ( unique ) return ( (i&1) ? 2*(q+1) - e : e );                           // if bsgs finds a unique multiple of the order of (x,y), we are done (this almost always happens!)
        // we know e divides q+1+/-t, thus t = +/-(q+1) mod e (sign is + when i is even, - ow)
        m1 = e;  a1 = (q+1)%e;  if ( (i&1) ) a1 = e-a1;
        if ( ! i_aseq_intersection (&m, &a, m, a, m1, a1) || ( a > 2*p && a-m < -2*p ) ) {
            printf ("No possible value of t found a=%ld, m=%ld, with f1=%ld*z+%ld, f0=%ld*z+%ld, p=%ld, z^2=%ld\n",
                        a, m, _ff_get_ui (f[3]), _ff_get_ui (f[2]), _ff_get_ui (f[1]), _ff_get_ui (f[0]), _ff_p, _ff_get_ui(_ff_2g));
            abort();
        }
        // if a or a-m is the only integer congruent to a mod m in the Hasse interval, then we know the trace and are done
        if ( a+m > 2*p && a-m < -2*p ) return q+1-a;
        if ( a > 2*p && a-2*m < -2*p ) return q+1-a+m;
        if ( ! i ) {
            // compute twist coefficients g0 and g1 (we often won't get this far, which is why we wait to compute them)
            ff2_nonresidue (x);                                         // get a nonresidue a in Fp2 (stored in x)
            ff2_square (y, x);  ff2_mult (x, x, y);
            ff2_mult (g, f, x);  ff2_mult (g+2, f+2, y);                        // y^2 = x^3 + a^2f1x + a^3f0 is the quadratic twist 
        }
        h = ( (i&1) ? f : g );                                              // alternate between input curve and its twist based on parity of i
    }
}
