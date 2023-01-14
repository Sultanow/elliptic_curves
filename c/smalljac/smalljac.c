#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <gmp.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "mpzpolyutil.h"
#include "polyparse.h"
#include "jac.h"
#include "jacorder.h"
#include "ecurve_ff2.h"
#include "hecurve.h"
#include "smalljac.h"
#include "smalljactab.h"
#include "smalljac_g23.h"
#include "smalljac_internal.h"
#include "pointcount.h"
#include "hcpoly.h"
#include "igusa.h"
#include "nfpoly.h"
#include "cstd.h"

/*
    Copyright (c) 2007-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

/*
	Main smalljac module. 
	
	Support for curves of genus 1 or 2 over Q is essentially complete.  Still work to be done
	on supporting number fields, but at least quadratic extensions in genus 1 are supported,
	and higher degree fields in genus 1 and 2 are supported at degree 1 primes.
	
	Lots of work still to do in genus 3, we currently only support point-counting to compute
	Frobenius traces of hyperelliptic and some special plane quartics (e.g. Picard curves)
*/

static inline void smalljac_curve_init_f (smalljac_curve *sc) 
	{ while ( sc->f_inits <= sc->degree ) mpz_init (sc->f[sc->f_inits++]); }
static inline void smalljac_curve_init_FH (smalljac_curve *sc) 
	{ while ( sc->F_inits <= sc->degree ) mpz_init (sc->F[sc->F_inits++]);  while ( sc->H_inits <= (sc->degree+1)/2) mpz_init (sc->H[sc->H_inits++]); }
static inline void smalljac_curve_init_Deltas(smalljac_curve *sc)
	{ while ( sc->Delta_inits <= sc->degree ) mpz_init (sc->Deltas[sc->Delta_inits++]); }


static int _smalljac_initted;
void smalljac_init (void)
{
	if ( _smalljac_initted ) return;
	smalljac_table_alloc (SMALLJAC_TABBITS);
	pointcount_init (SMALLJAC_INIT_COUNT_P);
	_smalljac_initted = 1;
}

long _smalljac_curve_points (smalljac_curve_t c, unsigned long p, int n, int affine)
{
	smalljac_curve *sc = c;
	unsigned long f[SMALLJAC_MAX_DEGREE+1];
	long D[SMALLJAC_MAX_DEGREE+1];
	long q, ipts;
	int d, i;

	if ( ! sc->Qflag ) return SMALLJAC_NOT_OVER_Q;
	if ( n < 1 || n > 3 || p > SMALLJAC_INIT_COUNT_P ) return SMALLJAC_INVALID_PP;
	if ( sc->type != SMALLJAC_CURVE_ELLIPTIC && sc->type != SMALLJAC_CURVE_HYPERELLIPTIC ) return SMALLJAC_UNSUPPORTED_CURVE;
	if ( p == 2 ) {
		ff2k_t f[SMALLJAC_MAX_DEGREE+1];
		ff2k_t h[SMALLJAC_MAX_DEGREE/2+1];
		int df, dh;
		
		ff2k_setup (n);
		df = ff2k_poly_set_mpz (f, sc->F, sc->dF);
		dh = ff2k_poly_set_mpz (h, sc->H, sc->dH);
		if ( df <= 2 && dh <= 1 ) return SMALLJAC_UNSUPPORTED_CURVE;	// don't deal with genus 0 cases in char 2
		return ff2k_hyperelliptic_pointcount (f, df, h, dh) - (affine?1:0)*ff2k_hyperelliptic_pointcount_at_infinity(df,dh,n);
	}
	ui_poly_set_mpz_mod_p (f, sc->f, sc->degree, p);
	for ( d = sc->degree ; d >= 0 && f[d]==0 ; d-- );
	if ( d <= 2 ) {	// genus 0 cases
		for ( q = p, i = 1 ; i < n ; q *= p, i++ );
		switch (d) {
		case -1: return q+(affine?0:1);									// y^2 = 0, double line
		case 0:  if ( !(n&1) || legendre(f[0],p) == 1 ) return 2*q+(affine?0:1);		// union of two rational lines that intersect at infinity
			     else return (affine?0:1);									// union of two irrational lines, intersection at infinity is rational
		case 1: return q+(affine?0:1);									// irreducible conic with 1 pt at infinity
		case 2:
			D[0] = (f[1]*f[1]-4*f[2]*f[0]) % p;								// we assume p < SMALLJAC_INIT_COUNT_P < 2^31 so this fits in a long
			if ( ! D[0] ) {												// y^2=(a*x+b)^2, union of two lines with affine intersection (-b/a:0:1) (so 2 pts at infinity)
				if ( legendre(f[2],p) == 1 ) return 2*q-1+(affine?0:2);			// lines are rational
				else return 1;											// lines are irrational, affine intersection is rational
			} else {													// irreducible conic with 0 or 2 pts at infinity
				if ( !(n&1) ) ipts = 2; else ipts = 1+ legendre(f[2], p);
				return q+1-(affine?ipts:0);
			}
		default:
			printf ("Unhandled case in _smalljac_curve_points"); assert(0);
		}
	}
	if ( (d&1) ) { ipts = 1; } else { if ( !(n&1) ) ipts = 2; else ipts = 1+ legendre(f[d],p); }
	if ( p <= d ) {
		switch (n) {
		case 1: return pointcount_tiny (f, d, p)-(affine?ipts:0);
		case 2: return pointcount_tiny_d2 (f, d, p)-(affine?ipts:0);
		case 3: return pointcount_tiny_d3 (f, d, p)-(affine?ipts:0);
		default:  return SMALLJAC_INVALID_PP;
		}
	}
	if ( n == 3 && p > POINTCOUNT_MAX_TINYP ) return SMALLJAC_INVALID_PP;
	if ( n == 2 ) {
		switch (d) {
		case 3: return pointcount_g1_d2 (f,p)-(affine?ipts:0);
		case 4: return pointcount_g1d4_d2 (f,p)-(affine?ipts:0);
		case 5: return pointcount_g2_d2 (f,p)-(affine?ipts:0);
		case 6: return pointcount_g2d6_d2 (f,p)-(affine?ipts:0);
		default: return SMALLJAC_UNSUPPORTED_CURVE;
		}
	}
	pointcount_precompute_long (D, (long *)f, d);
	for ( i = 0 ; i <= d ; i++ ) { D[i] %= (long)p; if ( D[i] < 0 ) D[i] += p; }
	switch (d) {
	case 3: return pointcount_g1 ((unsigned long *)D,p)-(affine?ipts:0);
	case 4: return pointcount_g1d4 ((unsigned long *)D,p,f[4])-(affine?ipts:0);
	case 5: return pointcount_g2 ((unsigned long *)D,p)-(affine?ipts:0);
	case 6: return pointcount_g2d6 ((unsigned long *)D,p,f[6])-(affine?ipts:0);
	default: return SMALLJAC_UNSUPPORTED_CURVE;
	}
}

long smalljac_curve_points (smalljac_curve_t c, unsigned long p, int n) { return _smalljac_curve_points (c, p, n, 0); }
long smalljac_curve_affine_points (smalljac_curve_t c, unsigned long p, int n) { return _smalljac_curve_points (c, p, n, 1); }



// callback function used by smalljac_Lpoly, just copies n and a
int smalljac_Lpoly_callback (smalljac_curve_t curve, unsigned long p, int good, long a[], int n, void *arg)
	{ int i; if ( good ) { ((long*)arg)[0] = n;  for ( i = 0 ; i < n ; i++ ) ((long*)arg)[i+1] = a[i]; } else { ((long*)arg)[0] = 0; } return 1; }

/*
	For the sake of simplicity, smalljac_Lpoly is implemented via a call to smalljac_Lpolys with start=end=p.
	This is slightly slower, but performance is not a major consideration for an individual call to smalljac_Lpoly.

	Those interested in computing Lpolys for large sets of curves will almost certainly want to use the underlying
	functionality directly (e.g. pointcount.c, which supports multi-curve point-counting, for example).
*/
int smalljac_Lpoly (long a[], char *curve, unsigned long q, unsigned long flags)
{
	mpz_t P;
	smalljac_curve *sc;
	unsigned long p;
	long b[SMALLJAC_MAX_GENUS+1];
	long sts;
	int i, h, n, error;

	mpz_init (P);
	mpz_set_ui (P, q);
	h = mpz_pp_base (P, P);
	p = mpz_get_ui (P);
	mpz_clear (P);
		
	if ( ! h ) return SMALLJAC_INVALID_PP;
	if ( h > 1 ) {
		if ( (flags&SMALLJAC_GROUP) ) return SMALLJAC_INVALID_PP;			// group computation is currently supported only for prime fields
		if ( (flags&SMALLJAC_A1_ONLY) ) flags &= ~SMALLJAC_A1_ONLY; 		// turn off A1_ONLY flag if q is a prime power
	}

	sc = smalljac_curve_init (curve, &error);
	if ( ! sc ) return error;
	if ( ! sc->Qflag ) { smalljac_curve_clear (sc); return SMALLJAC_NOT_OVER_Q; }
	if ( p > smalljac_curve_max_p(sc) && ! sc->special ) { smalljac_curve_clear (sc); return SMALLJAC_INVALID_PP; }
	b[0] = 0;

	sts = smalljac_Lpolys (sc, q, q, flags, smalljac_Lpoly_callback, b);
	smalljac_curve_clear (sc);
	
	if ( sts <= 0 ) return sts;
	n = b[0];
	if ( ! n ) return 0;
	for ( i = 0 ; i < n ; i++ ) a[i] = b[i+1];
	
	if ( h > 1 )  if ( ! smalljac_Lpoly_extend (a, n, p, h) ) return SMALLJAC_INTERNAL_ERROR;
	return n;
}

// parse the first line of the input file to smalljac_Lpolys_from_file(s) to extract the curve string and norm range [start,end], and optionally the genus
// does not try to validate, assumes the format is valid and only complains if it can't get the info it needs
// modifies buf inplace to contain null-terminated curve string
// returns 1 for success, 0 for failure
int smalljac_parse_file_header (char buf[SMALLJAC_CURVE_STRING_LEN+256], long *pstart, long *pend, int *pgenus)
{
	register char *s;

	if ( buf[0] != '[' ) { fprintf (stderr, "Expected first character [ missing from file header\n"); return 0; }
	for ( s = buf+strlen(buf)-1 ; s > buf && *s != ']'  && *s != ')' ; s-- );
	if ( *s !=']' && *s != ')' ) return 0;
	s++;  *s++ = '\0';
	while ( isspace(*s) ) s++;
	if ( *s == 'g' ) { while ( *s && ! isdigit(*s) ) s++;  *pgenus = atoi(s);  while ( isdigit(*s) ) s++; while ( isspace(*s) ) s++; } else *pgenus = 0;
	*pstart = atol(s);
	while ( isdigit(*s) ) s++; while ( isspace(*s) ) s++;
	*pend = atol(s);
	return 1;
}

// parse data line of input file to smalljac_Lpolys_from_file(s), which should be in the form "q,a1,a2,...a_n\n" or "q,?",
// where the integer q is a prime power and the integers a_1,...a_n are L-poly coefficients, with either n=1 or n=g.  ? indicates bad reduction
// returns n for success, 0 for bad reduction, or -1 for format error
static int smalljac_parse_file_record (long *q, long a[SMALLJAC_MAX_GENUS], int n, char buf[])
{
	register char *s;
	int i;
	
	*q = atol (buf);
	if ( ! *q ) return -1;
	s = buf;
	for ( i = 0 ; i < n ; i++ ) {
		for ( ; *s && *s != ',' ; s++ );
		if ( !*s++ ) return -1;
		if ( *s == '?' ) return 0;
		a[i] = atol(s);
	}
	return n;
}

// simulate smalljac_Lpolys using precomputed data -- MAKES NO ATTEMPT TO VALIDATE DATA
// expected format is [curve] in the first line of the file (and possbily several subsequent lines) followed by lines of text q,a_1,a_2,...,a_n\n
// where q,a_1,...,a_n are integers, with a_1,...a_n representing L-poly coefficients for the curve at a prime of norm q
// lines begining q,? indicate primes of bad reduction
long smalljac_Lpolys_from_file (char *filename, unsigned long start, unsigned long end, unsigned long flags,
						int (*callback)(smalljac_curve_t curve, unsigned long q, int good, long a[], int n, void *arg), void *arg)
{
	smalljac_curve *sc;
	char buf[SMALLJAC_CURVE_STRING_LEN+256];
	long q, fstart, fend;
	FILE *fp;
	int quit, err, sts, g, n;
	
	fp = fopen (filename, "r");
	if ( ! fp ) return SMALLJAC_FILENOTFOUND;
	
	if ( start <= 0 ) start = 1;
	
	// get curve string, and file start and end norms from the first line of the file
	if ( ! fgets(buf,sizeof(buf),fp) || buf[0] != '[' ) { err_printf ("curve spec missing from first line of data file %s\n", filename); fclose(fp); return SMALLJAC_BADFILE; }
	if ( ! smalljac_parse_file_header (buf, &fstart, &fend, &g) ) { err_printf ("unable to find curve spec and prime range in the first line of the specified field %s\n", filename); fclose(fp); return SMALLJAC_BADFILE; }
	if ( start < fstart ) { info_printf ("Increasing requested start from %ld to %ld\n", start, fstart); start = fstart; }
	if ( end > fend )  { info_printf ("Decreasing requested end from %ld to %ld\n", end, fend); end = fend; }
	if ( g ) {
		sc = smalljac_curve_alloc();
		strcpy (sc->str, buf);
		sc->genus =g;
	} else {
		sc = smalljac_curve_init (buf, &err);
		if ( ! sc ) { err_printf ("unable to parse curve string in file %s\n%s", filename, buf); fclose(fp); return err; }
	}
	quit = 0;
	n =  ( (flags&SMALLJAC_A1_ONLY) ? 1 : sc->genus );
	while ( fgets(buf,sizeof(buf),fp) ) {
		if ( buf[0] == '[' ) continue;	// ignore repeated curve spec -- handy for merged lpdata files
		sts = smalljac_parse_file_record (&q, sc->a, n, buf);
		if ( sts < 0 ) { err_printf ("smalljac file format error in %s\n%s", filename, buf);  fclose (fp); return SMALLJAC_BADFILE; }
		if ( q < start ) continue;  if ( q > end ) break;
		if ( (flags & SMALLJAC_FILTER) )  if ( ! (*callback)  (sc, q, -1, 0, 0, arg) ) continue;
		if ( ! sts ) { if ( !(flags&SMALLJAC_GOOD_ONLY) ) if ( ! (*callback) (sc, q, 0, 0, 0, arg) ) { quit = 1; break; } continue; }
		sc->q = q; sc->n = n;
		if ( ! (*callback) (sc, sc->q, 1, sc->a, sc->n, arg) ) { quit = 1; break; }
	}
	if ( ! quit ) q = end;
	smalljac_curve_clear (sc);
	fclose (fp);
	return q;
}

// simulate smalljac_Lpolys using precomputed data stored in multiple files (typically created using the lpdata program) -- MAKES NO ATTEMPT TO VALIDATE DATA
// note: unlike smalljac_Lpoly or smalljac_Lpoly_from_file, there is no implicit ordering among degree-1 primes of the same norm
long smalljac_Lpolys_from_files (char *fileprefix, int jobs, unsigned long start, unsigned long end, unsigned long flags,
						 int (*callback)(smalljac_curve_t curve, unsigned long q, int good, long a[], int n, void *arg), void *arg)
{
	smalljac_curve *sc;
	char filename[1024];
	char buf[SMALLJAC_CURVE_STRING_LEN+256];
	char cstr[SMALLJAC_CURVE_STRING_LEN];
	long q[SMALLJAC_MAX_JOBS];
	long a[SMALLJAC_MAX_JOBS][SMALLJAC_MAX_GENUS];
	FILE *fp[SMALLJAC_MAX_JOBS];
	long fstart, fend, qend;
	int sts, quit, err, genus;
	register int i, j, n;
	
	if ( strlen(fileprefix) + 64 > sizeof(filename) ) { err_printf ("specified file prefix is too long for buffer\n"); return SMALLJAC_INTERNAL_ERROR; }
	if ( jobs > SMALLJAC_MAX_JOBS ) { err_printf ("jobs cannot exceed SMALLJAC_MAX_JOBS = %d\n", SMALLJAC_MAX_JOBS); return SMALLJAC_INTERNAL_ERROR; }
	
	if ( start <= 0 ) start = 1;
	
	err = 0; qend = 0; sc = 0; n = 0;
	for ( i = 0 ; i < jobs ; i++ ) fp[i] = 0;
	for ( i = 0 ; i < jobs ; i++ ) {
		sprintf (filename, "%s_%d_%d.txt", fileprefix, jobs, i);
		fp[i] = fopen (filename, "r");
		if ( ! fp[i] ) { err_printf ("Error opening file %s\n", buf); err = SMALLJAC_FILENOTFOUND; goto done; }
		if ( ! fgets(buf,sizeof(buf),fp[i]) || buf[0] != '[' ) { err_printf ("curve spec missing from first line of data file %s\n", filename); err = SMALLJAC_NODATA; goto done; }
		if ( ! smalljac_parse_file_header (buf, &fstart, &fend, &genus) ) { err_printf ("unable to find curve spec and prime range in the first line of the specified field %s\n", filename); err = SMALLJAC_BADFILE; goto done; }
		if ( start < fstart || end > fend ) { err_printf ("Requested data range [%ld,%ld] extends outside the data range [%ld,%ld] in the file %s\n", start, end, fstart, fend, filename);  err = SMALLJAC_NODATA; goto done; }
		if ( ! i ) {
			strcpy (cstr, buf);
			if ( genus ) {
				sc = smalljac_curve_alloc();
				strcpy (sc->str, buf);
				sc->genus = genus;
			} else {
				sc = smalljac_curve_init (cstr, &err);
				if ( ! sc ) { err_printf ("unable to parse curve string in file %s\n%s", filename, buf); fclose(fp[i]); return err; }
			}
			n = (flags&SMALLJAC_A1_ONLY) ? 1 : sc->genus;
		} else {
			if ( strcmp (buf, cstr) != 0 ) { err_printf ("inconsistent curve string in file %s\n%s", filename, buf);  err = SMALLJAC_BADFILE; goto done; }
		}
		if ( ! fgets(buf,sizeof(buf),fp[i]) ) { q[i] = 0; continue; }
		sts = smalljac_parse_file_record (q+i, a[i], n, buf);
		if ( sts < 0 ) { err_printf ("smalljac file format error in %s\n%s", filename, buf);  err = SMALLJAC_BADFILE; goto done; }
		if ( ! sts ) a[i][0] = JAC_INVALID_A1;		// indicate bad reduction with an invalid a1 value
	}

	quit = 0;
	for (;;) {
		// get data for the next record, ordered by q
		for ( i = 0 ; i < jobs && ! q[i] ; i++ );
		if ( i == jobs ) break;
		for ( j = 1 ; j < jobs ; j++ ) if ( q[j] && q[i] > q[j] ) i = j;
		sc->q = q[i];  sc->n = ( a[i][0] == JAC_INVALID_A1 ? 0 : n );
		
		// replace the data we just used
		for ( j = 0 ; j < sc->n ; j++ ) sc->a[j] = a[i][j];
		if ( ! fgets(buf,sizeof(buf),fp[i]) ) {
			q[i] = 0;
		} else {
			sts = smalljac_parse_file_record (q+i, a[i], n, buf);
			if ( sts < 0 ) { err_printf ("smalljac file format error in %s\n%s", filename, buf);  err = SMALLJAC_BADFILE; goto done; }
			if ( ! sts ) a[i][0] = JAC_INVALID_A1;		// indicate bad reduction with an invalid a1 value
		}
			
		// process the data
		if ( sc->q < start ) continue;  if ( sc->q > end ) { q[i] = 0; continue; }
		if ( (flags & SMALLJAC_FILTER) )  if ( ! (*callback)  (sc, sc->q, -1, 0, 0, arg) ) continue;
		if ( ! sc->n ) {
			if ( !(flags&SMALLJAC_GOOD_ONLY) ) if ( ! (*callback) (sc, sc->q, 0, 0, 0, arg) ) { quit = 1; break; }
			continue;
		}
		if ( ! (*callback) (sc, sc->q, 1, sc->a, sc->n, arg) ) { quit = 1; break; }
	}
	qend = ( quit ? q[i] : end );
done:
	if ( sc ) smalljac_curve_clear (sc);
	for ( i = 0 ; i < jobs ; i++ ) if ( fp[i] ) fclose (fp[i]);
	return (err ? err : qend );
}


// reduces the poly defined over a number field at the ith degree-e prime ideal above p
// you must setup Fp (via ff_setup)  and call smalljac_nf_reduce_setup first,
// which (arbitrarily) determines the ordering of the prime ideals.
// output is a poly f(x) over Fp^n, return value is 0 if y^2=f(x) is found to be singular
// (but this is typically not checked here, other than ensuring the poly has large enough degree)
static inline int smalljac_reduce_nf_poly (smalljac_curve *sc, long p, int e, int i)
{
	register int d;

	sc->hc->n = e;
	sc->hc->d = d = nf_poly_reduce (sc->hc->f, sc->nfp, p, e, i);
	if  ( p == 2 ) return ( d ? 1 : 0 );
	if ( d < 2*sc->genus+1 ) return 0;
	if ( e == 1 ) if ( !(d&1) || ! _ff_one(sc->hc->f[d]) || ! _ff_zero(sc->hc->f[d-1]) ) hc_poly_standardize (sc->hc);
	// don't check discriminant here, it will be checked in smalljac_internal_Lpoly
	return 1;	
}


long smalljac_Lpolys (smalljac_curve_t curve, unsigned long start, unsigned long end, unsigned long flags,
				   int (*callback)(smalljac_curve_t curve, unsigned long p, int good, long a[], int n, void *arg), void *arg)
{
	static mpz_t P, D;
	static int init;
	smalljac_curve *sc;
	prime_enum_ctx_t *ctx;
	unsigned long h[SMALLJAC_MAX_BAD_PRIMES], badp[SMALLJAC_MAX_BAD_PRIMES];
	register unsigned long p, d, pbitmask, pbits;
	long window;
	int e, i, k, filter, good, good_only,  error, badpi, badpk;

	if ( ! init ) { smalljac_init();  mpz_init (P);  mpz_init (D); init = 1; }
	if ( (flags&SMALLJAC_A1_ONLY) && (flags&SMALLJAC_GROUP) ) return SMALLJAC_INVALID_FLAGS;
	sc = (smalljac_curve *)curve;
	if ( end < start || end > smalljac_curve_max_p(sc) ) { printf ("start=%lu, end=%lu, maxp=%lu\n", start, end, smalljac_max_p(sc->genus)); return SMALLJAC_INVALID_INTERVAL; }

	if ( flags&SMALLJAC_GROUP && !(sc->degree&1) ) { err_printf ("Currently group computations are supported only for hyperelliptic curves of the form y^2=f(x) with deg f = 2g+1 odd\n"); return SMALLJAC_UNSUPPORTED_CURVE; }
		
	if ( sc->genus > SMALLJAC_GENUS && ! (flags&SMALLJAC_A1_ONLY) && ! sc->special ) { err_printf ("The SMALLJAC_A1_ONLY flag must be set for curves of genus %d (or change SMALLJAC_GENUS and recompile)\n", sc->genus); return SMALLJAC_UNSUPPORTED_CURVE; }
	
	good_only = (flags&SMALLJAC_GOOD_ONLY);
	filter = (flags&SMALLJAC_FILTER);
	
	pbitmask = (flags&SMALLJAC_SPLIT)>>SMALLJAC_SPLIT_SHIFT;
	pbits = (flags&SMALLJAC_HIGH)>>SMALLJAC_HIGH_SHIFT;
	
	if ( sc->genus != 1 && flags&SMALLJAC_PRIME_ORDER ) { err_printf ("SMALLJAC_PRIME_ORDER flag only supported in genus 1\n");  return SMALLJAC_INVALID_FLAGS; }
	/*
		The standard (and most optimized) case is that we have a curve that is defined over Q, possibly being considered over a larger number field, but then only at degree-1 primes.
		It suffices to compute the Lpoly in Fp (just once per p), and then make multiple callbacks to account for the number of degree-1 primes above p (exactly 1 in Q)
		
		The code below is a bit unpleasant (compare to the number field case, which is much cleaner) and should probably be re-worked at some point (much of the unpleasantness
		has to do with trying to make the typical case (good reduction at a largish prime) as fast as possible, while still handling all the special cases.
	*/
	if  ( sc->Qflag && (sc->nfd == 1 || (flags&SMALLJAC_DEGREE1_ONLY)) ) {
		// Precompute delta values for pointcounting, if needed (but only when curve is defined over Q)
		if ( (start <= smalljac_count_p(sc->genus) || sc->genus > 2) && ! (sc->flags&SMALLJAC_CURVE_FLAG_DELTA) ) {
			smalljac_curve_init_Deltas(sc);
			pointcount_precompute (sc->Deltas, sc->f, sc->degree);
			sc->flags |= SMALLJAC_CURVE_FLAG_DELTA;
		}
		
		// If D is small, factor it to save time on bad reduction checks (but note that D is currently only used for curves over Q)
		badpk = badpi = 0;
		i = mpz_sizeinbase(sc->D,2);
		if ( i < 64 && i < 3*ui_len(end-start) ) badpk = ui_factor(badp,h,mpz_get_ui(sc->D)); else badpk = 0;
		badpi = 0;

		error = 0;  window = 0;
		// use fast prime enumeration (based on a wheeled sieve) in all cases (now supports up to 2^40, as of Feb 2010)
		if ( sc->genus==1 && (flags & SMALLJAC_PRIME_ORDER) ) {
			if ( flags&SMALLJAC_LOW_ORDER ) window = -2*sqrt(end);
			else window = 4*sqrt(end);
		}
		k = 1;
		ctx = fast_prime_enum_start_w (start, end, window);
		while ( (p = fast_prime_enum(ctx)) ) {
			if ( (p&pbitmask) != pbits ) continue;
			sc->q = p;
			if ( filter && ! (*callback) (curve, sc->q, -1, 0, 0, arg) ) continue;
			if ( badpk ) {
				while ( badpi < badpk && badp[badpi] < p ) badpi++;
				if ( badpi < badpk && p==badp[badpi] ) good = 0; else good = 1;			
			} else {
				good = ! mpz_divisible_ui_p (sc->D,p);
			}
			if ( ! good ) {
				if ( ! good_only ) {
					if ( (sc->flags&SMALLJAC_CURVE_FLAG_WS) ) {
						// sc->f holds curve in short ws form
						mpz_mul (D, sc->f[0], sc->f[1]);
						mpz_mul_2exp (D, D, 1);
						mpz_neg (D, D);
						sc->a[0] = -mpz_kronecker_ui (D, p);		// at bad primes the Frobenius trace is (-2AB/p) for y^2=x^3+Ax+B
						sc->n = 1;
						if ( (sc->flags&SMALLJAC_GROUP) ) sc->a[0] += p;
						if ( ! (*callback) (curve, p, 0, sc->a, sc->n, arg) ) break;	
					} else {
						if ( ! (*callback) (curve, p, 0, 0, 0, arg) ) break;
					}
				}
				continue;
			}
			if ( sc->nfd > 1 ) {
				k = nf_poly_reduce_setup (sc->nfp, p, 1);
				if ( ! k ) continue;
			}
			// handle tiny primes separately
			if ( p <= smalljac_tiny_p(sc->genus) ) {
				sc->n = smalljac_tiny_Lpoly (sc->a, sc, p, flags);	// returns -1 for bad reduction but will still compute ap=-1,0,1 correctly for genus 1 curves in weierstrass form at 2 and 3
				if ( sc->n < -1 ) return SMALLJAC_INTERNAL_ERROR;
				if ( sc->n < 0 ) {
					if ( ! good_only ) if ( ! (*callback) (curve, p, 0, sc->a, (sc->genus == 1 ? 1 : 0), arg) ) break;
					continue;
				}
				if ( sc->genus==1 ) {
					if ( (flags&SMALLJAC_LOW_ORDER) && sc->a[0] >= -1 ) continue;
					if ( (flags&SMALLJAC_PRIME_ORDER) && ! ui_is_prime(p+1+sc->a[0]) ) continue;	// note a[0] is the negated trace of Frobenius
				}
				// callback once for each degree-1 prime
				for ( i = 0 ; i < k ; i++ ) if ( ! (*callback) (curve, sc->q, 1, sc->a, sc->n, arg) ) break;
				if ( i < k ) break;
				continue;
			}
			// if only prime order groups are requested in genus=1, use Stickelberger to quickly rule out cases that must have a point of order 2
			if ( (flags& SMALLJAC_PRIME_ORDER) && sc->genus==1 ) { d = mpz_fdiv_ui(sc->disc,p);  if ( ui_legendre(d,p) < 0 ) continue; }
			sc->n = smalljac_internal_Lpoly_Q (sc->a, sc, p, flags);
			if ( sc->n == 0 ) continue;																	// indicates this case is excluded by flag settings (e.g. non-prime group order)
			if ( sc->n == -1 ) { if ( ! good_only ) if ( ! (*callback) (curve, sc->q, 0, sc->a, sc->n, arg) ) break; continue; }	// currently this can never happen
			if ( sc->n < -1 ) { error = 1;  break; }
			// callback once for each degree-1 prime
			for ( i = 0 ; i < k ; i++ ) if ( ! (*callback) (curve, sc->q, 1, sc->a, sc->n, arg) ) break;
			if ( i < k ) break;
		}
		fast_prime_enum_end (ctx);
		if ( ! p ) p = end;
		if ( error ) { printf ("smalljac internal error at p=%lu\n", p);  return SMALLJAC_INTERNAL_ERROR; }
		return (long) p;
	}

	if ( (flags&SMALLJAC_PRIME_ORDER) ) { err_printf ("Flag SMALLJAC_PRIME_ORDER is only supported for curves defined over Q (or at degree-1 primes)\n"); return SMALLJAC_INVALID_FLAGS; }
	if ( (flags&SMALLJAC_GROUP) ) { err_printf ("Flag SMALLJAC_GROUP is only supported for curves defined over Q (or at degree-1 primes)\n"); return SMALLJAC_INVALID_FLAGS; }
	
	if ( sc->nfd == 2 ) {
		if ( sc->genus != 1 ) { err_printf ("Curves over quadratic fields are currently only supported in genus 1 or when SMALLJAC_DEGREE1_ONLY is set\n"); return SMALLJAC_INTERNAL_ERROR; }
		error = 0;
		ctx = fast_prime_enum_start (start, end, 2);
		while ( (p = fast_prime_enum_powers(ctx)) ) {
			if ( (p&pbitmask) != pbits ) continue;
			sc->q = p;
			if ( filter && ! (*callback) (curve, sc->q, -1, 0, 0, arg) ) continue;
			e = fast_prime_enum_exp (ctx);
			if ( e > 1 ) p = fast_prime_enum_base (ctx);
			k =nf_poly_reduce_setup (sc->nfp, p, e);
			for ( i = 0 ; i < k ; i++ ) {
				if ( ! smalljac_reduce_nf_poly (sc, p, e, i) ) {
					if ( ! good_only ) if ( ! (*callback) (curve, sc->q, 0, 0, 0, arg) ) break;
					continue;
				}
				sc->n = smalljac_internal_Lpoly_nf (sc->a, sc->hc, p, flags);
				if ( sc->n == 0 ) continue;
				if ( sc->n == -1 ) { if ( ! good_only ) if ( ! (*callback) (curve, sc->q, 0, 0, 0, arg) ) break; continue; }
				if ( sc->n == -1 ) { error = 1; break; }
				if ( ! (*callback) (curve, sc->q, 1, sc->a, sc->n, arg) ) break;
			}
			if ( i < k ) break;
		}
		fast_prime_enum_end (ctx);
		if ( ! p ) p = end;
		if ( error ) { printf ("smalljac internal error at p=%lu\n", p);  return SMALLJAC_INTERNAL_ERROR; }
		return (long) p;
	}
	
	if ( sc->nfd > 2 ) {
		if ( ! (flags&SMALLJAC_DEGREE1_ONLY) ) { err_printf ("SMALLJAC_DEGREE1_ONLY flag must be set for number fields of degree > 2"); return SMALLJAC_INVALID_FLAGS; }
		error = 0;
		ctx = fast_prime_enum_start (start, end, 0);
		while ( (p = fast_prime_enum(ctx)) ) {
			if ( (p&pbitmask) != pbits ) continue;
			sc->q = p;
			if ( filter && ! (*callback) (curve, sc->q, -1, 0, 0, arg) ) continue;
			k =nf_poly_reduce_setup (sc->nfp, p, 1);
			for ( i = 0 ; i < k ; i++ ) {
				if ( ! smalljac_reduce_nf_poly (sc, p, 1, i) ) {
					if ( ! good_only ) if ( ! (*callback) (curve, sc->q, 0, 0, 0, arg) ) break;
					continue;
				}
				sc->n = smalljac_internal_Lpoly_nf (sc->a, sc->hc, p, flags);
				if ( sc->n == 0 ) continue;
				if ( sc->n == -1 ) { if ( ! good_only ) if ( ! (*callback) (curve, sc->q, 0, 0, 0, arg) ) break; continue; }
				if ( sc->n == -1 ) { error = 1; break; }
				if ( ! (*callback) (curve, sc->q, 1, sc->a, sc->n, arg) ) break;
			}
			if ( i < k ) break;
		}
		fast_prime_enum_end (ctx);
		if ( ! p ) p = end;
		if ( error ) { printf ("smalljac internal error at p=%lu\n", p);  return SMALLJAC_INTERNAL_ERROR; }
		return (long) p;
	}
	err_printf ("Fell through to unhandled case in smalljac_Lpolys sc->Qflag=%d, sc->nfd=%d, flags=%lx!\n", sc->Qflag, sc->nfd, flags); abort();
}

/*
	The following convention applies to the return value n of all the smalljac_*_Lpoly_* functions:

		n > 0:  either the number of Lpoly coefficients or cyclic factors in a[]
		n = 0:  no data returned due to flag settings (e.g. only prime order groups requested)
		n = -1: curve is singular
		n < -1: internal error of some sort (an explanatory error message will usually be printed at the point where the error occured)

	Depending on the function, not all of these cases apply, e. g. over Q singular curves (for primes > 3) are filtered out above (so -1 is never returned)
	and over number fields none of the flag restrictions are currently supported (so 0 is never returned).

	With smalljac_v4, we have switched all these function to use long rather than unsigned long for the prime p and the pointcount pts.
	This is to avoid easy mistakes that can arise when mixing signed and unsigned integers -- in fact it would be better if  _ff_p was a long also, but this hasn't changed (yet).
*/


/*
	The implementation of smalljac_internal_Lpoly_nf is very rudimentary at the moment, 
	it currently only handles genus 1 curves c over Fp or Fp^2 (and totally ignores flags). 
	All other cases get sent to smalljac_internal_Lpoly_Q.
*/
int smalljac_internal_Lpoly_nf (long a[], hc_poly *c, long p, unsigned long flags)
{
	ff_t d;
	long N;

	if ( c->d <= 4 ) {
		if ( p == 2 ) {
			if ( c->d != 4 ) { err_printf("Expected 5 Weierstrass coefficients for genus 1 curve in char 2 in smalljac_internal_Lpoly_nf, coeffs=%d\n", c->d+1); abort(); }
			N = ff2k_WS_pointcount (c->f);
			if ( N < 0 ) return -1;
			a[0] = N - ((1L<<_ff2k_k)+1);
			return 1;
		}
		assert ( p==_ff_p );
		assert ( c->d == 3 );
		switch ( c->n ) {
		case 1:
			if ( p == 3 ) { N = ecurve_order_F3 (0, c->f); if ( N < 0 ) return -1;  a[0] = N-(p+1); return 1; }
			assert (c->d==3 && _ff_one(c->f[3]) && _ff_zero (c->f[2]));
			if ( _ff_zero(c->f[1]) ) return smalljac_x3pa_Lpoly (a, c->f[0]);
			if ( _ff_zero(c->f[0]) ) return smalljac_x3pax_Lpoly (a, c->f[1]);
			if ( ! ff_poly_x3axb_disc(&d, c->f) ) return -1;		// check for singular curve here, since ecurve_order does not
			N = ecurve_order (0, c->f);				// N cannot be -1
			a[0] = N - (p+1);
			return 1;
		case 2:
			N = ecurve_ff2_group_order (c->f);
			if ( N < 0 ) return -1;
			a[0] = N - (p*p+1);
			return 1;
		default:
			err_printf ("Unsupported number field of degree %d in smalljac_internal_Lpoly\n", c->n); abort();
		}
	}
	ff_poly_print (c->f, c->d);
	printf ("Fell through to unimplemented case in smalljac_internal_Lpoly_nf, c->d = %d\n", c->d);
	return -2;
}

/*
	smalljac_internal_Lpoly_mod assumes p > sc->degree and good reduction.
        It computes the coefficients of L_p(T) for the curve specified by y^2=f(x) (defined over Q)
	The array D0 holds precomputed delta values used for pointcounting which will be
	used to determine a1 if p < smalljac_count_p.
	
	If p >= SMALLJAC_PADIC_P, a p-adic computation will be used to compute the 
	coefficients of L_p(T) modulo p (or p^N in genus > 3), using Harvey's frobenius
	code, followed by a generic group computation to get the exact function.
	
	Otherwise, a generic computation is performed, possibly using the a1 value
	obtained from pointcounting.
	
	The return value is either 0 (error) or the number of coefficients computed.
	If SMALLJAC_A1_ONLY is set and pointcounting is performed (not necessarily
	the case) then this will be 1, otherwise it will be g. 

	Note that for large p it can be more efficient to compute a1 via a method that
	computes all the coefficients rather than pointcounting, in this case all the
	coefficients are returned (may as well...).
	
	Return value is 0 if p is excluded by flags (e.g. not a prime order group),
	-2 for errors, or n > 0 for the number of coefficients or group rank
	Assumes curve is non-singular
*/
smalljac_curve **sc_ptr;
int smalljac_internal_Lpoly_Q (long a[], smalljac_curve *sc, long p, unsigned long flags)
{
	ff_t t;
	int n, sts;
	assert ( sc->Qflag );
	sc->pts = 0;

	// handle special curves that we can deal with more quickly -- in theory this is any CM curve, but for now we just deal with y^2=x^3+a and y^2=x^3*ax
	switch (sc->genus) {
	case 1:
		if ( sc->special == SMALLJAC_SPECIAL_CM && !(flags&SMALLJAC_GROUP) ) {
			sts = -3;
			if ( sc->degree == 3 && mpz_cmp_ui (sc->f[3],1)==0 && ! mpz_sgn(sc->f[2]) ) {
				if ( ! mpz_sgn(sc->f[0]) ) { ff_setup_ui (p); _ff_set_mpz (t, sc->f[1]); sts = smalljac_x3pax_Lpoly (a, t); }
				if ( ! mpz_sgn(sc->f[1]) ) { ff_setup_ui (p); _ff_set_mpz (t, sc->f[0]); sts = smalljac_x3pa_Lpoly (a, t); }
				// kludge to handle prime order filtering
				if ( sts != -3 ) {
					if ( sts == 1 ) {
						sc->pts = p+1+a[0];	// note a[0] is negated trace of Frobenius
						if ( (flags&SMALLJAC_LOW_ORDER) && sc->pts >= p ) return 0;
						if ( (flags&SMALLJAC_PRIME_ORDER) && ! ui_is_prime(sc->pts) ) return 0;
					}
					return sts;
				}
			}
		}
		break;
	case 2:
		// invoke special code for curves of the form y^2 = x^6+a, y^2=x^5 + ax, and y^2 = x^5 + a, and also twists of these curves
		if ( sc->special && !(flags&SMALLJAC_GROUP) ) {
			 ff_setup_ui(p); 
			switch (sc->special) {
			case SMALLJAC_SPECIAL_CM2SQUARE: return smalljac_cm2square_Lpoly (a);	// very special curve y^2=x^6-5x^4-5x^2+1, has Jacobian Q-isogenous to the square of an elliptic curve over Q with CM discriminant -2.
			case SMALLJAC_SPECIAL_X6PA: _ff_set_mpz (t, sc->f[0]);  return smalljac_x6pa_Lpoly (a, t);
			case SMALLJAC_SPECIAL_X5PAX: _ff_set_mpz (t, sc->f[1]);  return smalljac_x5pax_Lpoly (a, t);
			case SMALLJAC_SPECIAL_X5PA:
				switch (p%5) { case 2: case 3: a[0] = a[1] = 0; return 2; case 4: a[0] = 0; a[1] = 2*p; return 2; }		// only handle p != 1 mod 5, we could use LLL to deal with p = 1 mod 5 but we won't bother for now
				break;
			case SMALLJAC_SPECIAL_X5PX_TWIST:
			case SMALLJAC_SPECIAL_X6P1_TWIST:
				if ( p > 1024 ) {
					sc_ptr = &sc;	// this fixes a bizarre bug in which the pointer sc gets changed in the call below (optimization bug?)
					hc_poly_set_mpz (sc->hc, sc->f, sc->degree);
					if ( smalljac_genus2_charpoly_cmsquare (a, (sc->special == SMALLJAC_SPECIAL_X5PX_TWIST ? -2 : -3) , sc->hc) == 2 ) return 2;
					if ( sc->str ) puts (sc->str);
					assert ( p <= smalljac_max_p (sc->genus) );
				}
				break;
			}
		}
		break;
	case 3:
		if (  sc->special == SMALLJAC_SPECIAL_PICARD ) {
			unsigned long Deltaf0[SMALLJAC_MAX_DEGREE+1];
			if ( ! (flags&SMALLJAC_A1_ONLY) ) { err_printf ("Currently SMALLJAC_A1_ONLY must be set for Picard curves\n"); return -2; }
			ui_poly_set_mpz_mod_p (Deltaf0, sc->Deltas, sc->degree, p);
			sc->pts = pointcount_pd4 (Deltaf0,p);
			a[0] = (long)sc->pts - (p+1);
			return 1;
		}
		if ( (flags&SMALLJAC_A1_ONLY) && sc->special == SMALLJAC_SPECIAL_X7PAX && p > 144 ) { ff_setup_ui (p); _ff_set_mpz (t, sc->f[1]);  return smalljac_x7pax_a1 (a, t); }
		if ( (flags&SMALLJAC_A1_ONLY) && sc->special == SMALLJAC_SPECIAL_X8PA && p > 144 ) { ff_setup_ui (p); _ff_set_mpz (t, sc->f[0]);  return smalljac_x8pa_a1 (a, t); }
		if ( sc->special == SMALLJAC_SPECIAL_FK_TWIST ) { ff_setup_ui (p);   sts = smalljac_FKtwist_Lpoly (a, sc->special_curve_id);  if ( sts > 0 ) sc->pts = p+1+a[0];  return sts;}
		break;
	}

	if ( sc->genus > 2 ) {
		if ( (flags&SMALLJAC_A1_ONLY) ) {
			sc->pts = smalljac_pointcount_modp (sc, p);
			if ( ! sc->pts ) return -2;
			a[0] = (long)sc->pts - (p+1); 
			return 1;
		}
//		if ( p >= SMALLJAC_PADIC_P ) return smalljac_padic_Lpoly (a, sc, p, flags);
//		err_printf ("Currently unable to handle non-special curves of genus > 2 (hypellfrob not linked in)\n");
//		return -2;
	}

	if ( p <= smalljac_count_p(sc->genus) || sc->genus > 2 ) {
		sc->pts = smalljac_pointcount_modp (sc, p);
		if ( ! sc->pts ) { err_printf ("smalljac_pointcount failed at p=%lu\n", p); return -2; }
		if ( (flags&SMALLJAC_A1_ONLY) ||((flags&SMALLJAC_A1_ZERO) && sc->pts != p+1) ) { a[0] = (long)sc->pts - (p+1);  return 1; }
		if ( sc->genus == 1 ) {										// in genus 1, we're basically done, but we might need to do a group structure computation
			if ( (flags&SMALLJAC_LOW_ORDER) && sc->pts >= p ) return 0;
			if ( (flags&SMALLJAC_PRIME_ORDER) && ! ui_is_prime(sc->pts) ) return 0;
			// this code is used in genus 1 only for very small p
			if ( ! (flags&SMALLJAC_GROUP) ) { a[0] = (long)sc->pts - (p+1);  return 1; }
			ff_setup_ui (p);
			hc_poly_set_mpz (sc->hc, sc->f, sc->degree);
			n = ecurve_group_structure (a,sc->pts,0,sc->hc->f);			// we don't have any torsion info so specify d=0
			return ( n ? n : -2 );
		}
	}
	assert ( p <= smalljac_max_p (sc->genus) );
	ff_setup_ui (p);
	sc_ptr = &sc;	// this fixes a bizarre bug in which the pointer sc gets changed in the call below (optimization bug?)
	hc_poly_set_mpz (sc->hc, sc->f, sc->degree);
	return smalljac_generic_Lpoly (a, sc->hc, sc->pts, flags);
}

// point-counting over Fp (based on KS08)
// curve is assumed to be non-singular and defined over Q
unsigned long smalljac_pointcount_modp (smalljac_curve *sc,  long p)
{
	unsigned long Deltaf0[SMALLJAC_MAX_DEGREE+1];

	assert (sc->Qflag);
	ui_poly_set_mpz_mod_p (Deltaf0, sc->Deltas, sc->degree, p);
	if ( p < SMALLJAC_BIG_COUNT_P ) {
		switch (sc->degree) {
		case 3: return pointcount_g1 (Deltaf0,p);
		case 4: return pointcount_g1d4 (Deltaf0,p,mpz_fdiv_ui(sc->f[4],p));
		case 5: return pointcount_g2 (Deltaf0,p);
		case 6: return pointcount_g2d6 (Deltaf0,p,mpz_fdiv_ui(sc->f[6],p));
		case 7: return pointcount_g3 (Deltaf0,p);
		case 8: return pointcount_g3d8 (Deltaf0,p,mpz_fdiv_ui(sc->f[8],p));
		case 9: return pointcount_g4 (Deltaf0,p);
		case 10: return pointcount_g4d10 (Deltaf0,p,mpz_fdiv_ui(sc->f[10],p));
		}
	} else {
		switch (sc->degree) {
		// we should never get here in genus 1
		case 5: return pointcount_big_g2 (Deltaf0,p);
		case 6: return pointcount_big_g2d6 (Deltaf0,p,mpz_fdiv_ui(sc->f[6],p));
		case 7: return pointcount_big_g3 (Deltaf0,p);
		case 8: return pointcount_big_g3d8 (Deltaf0,p,mpz_fdiv_ui(sc->f[8],p));
		case 9: return pointcount_big_g4 (Deltaf0,p);
		}
	}
	err_printf ("Fell throught to unimplemented case in smalljac_pointcoint_modp, degree %d\n", sc->degree);
	return 0;
}

// Computes L-poly coefficients or group structure (or both) using generic group algorithms
// Returns -1 for errors, 0 if the prime is excluded by the settings of flags, and otherwise the
// number of entries in a which will be g for L-poly coefficients, or the rank of the group.
// the curve is assumed to be non-singular
int smalljac_generic_Lpoly (long a[], hc_poly *c, long pts, unsigned long flags)
{
	hc_poly twist;
	unsigned long Min, Max;
	unsigned long e, te, P1, PN1;
	unsigned long m;
	long a1, d, tcon, p;
	int constraints[SMALLJAC_MAX_CONSTRAINTS+1];
	int i, j, k, tk, n, sts;
	double x;

	p = _ff_p;
	// Use faster genus 1 code in ecurve instead of general jac functions
	if ( c->g == 1 ) {
		if ( pts ) {
			if ( (flags&SMALLJAC_LOW_ORDER) && pts >= p ) return 0;
			if ( (flags&SMALLJAC_PRIME_ORDER) && ! ui_is_prime(pts) ) return 0;
			a[0] = (long)pts - (long)(p+1);  return 1;
		}
		if ( flags & SMALLJAC_PRIME_ORDER ) {		
			P1 = ecurve_prime_order (c->f, ( (flags&SMALLJAC_LOW_ORDER) ? 1 : 0 ));
			if ( ! P1 ) return 0;
			if ( (flags&SMALLJAC_GROUP) ) { a[0] = P1; a[1] = 1; return 1; }
		} else {
			P1 = ecurve_order (&d,c->f);
			if ( (flags&SMALLJAC_LOW_ORDER) && P1 >= p ) return 0;
			if ( (flags&SMALLJAC_GROUP) ) { n = ecurve_group_structure (a,P1,d,c->f);  return ( n > 0 ? n : -2 ); }
		}
		a[0] = (long)(P1) - (long)(p+1);
		return 1;
	}
	if ( c->n != 1 ) { err_printf ("Finite fields of degree %d not supported for curves of genus %d\n", c->n, c->g);  return -2; }

	x = sqrt((double)p);
	Min = (unsigned long) (floor(pow(x-1.0, 2.0*c->g)));
	Max = (unsigned long) (ceil(pow(x+1.0, 2.0*c->g)));
	if ( pts ) {
		a1 = (long)pts - (long)(p+1);
	} else {
		a1 = JAC_INVALID_A1;
	}
	k = jac_order (&P1, Min, Max, a1, 0, flags&SMALLJAC_SGROUP, 0, c);
	if ( k <= 0 ) return 0;
	if ( k > SMALLJAC_MAX_CONSTRAINTS ) { printf ("%7lu: Exceeded SMALLJAC_MAX_CONSTRAINTS", p);  return -2; }
	if ( c->g == 2 ) {
		if ( k > 1 ) { printf ("%lu: Ambiguous result in genus 2 not handled\n", p);  return -2; }		// should be impossible provided pointcounting is used for p < SMALLJAC_TINY_P
		if ( (flags&SMALLJAC_GROUP) ) { n = jac_structure (a, c, P1, flags&SMALLJAC_SGROUP);  return ( n > 0 ? n : -2 ); }
		if ( pts ) {
			a[0] = (long)pts - (long)(p+1);
			a[1] = (long)(P1) - ((long)p*p+1)-(long)(p+1)*a[0];
			return 2;
		} else {
			if ( ! smalljac_genus2_charpoly_from_P1 (a, P1, Min, Max, c) ) {
				// If we can't deduce the charpoly from P1 (very rare) go ahead and compute PN1 from scratch
				// This is a bit silly, since we know the value of PN1 mod 2(p+1) and could use this knowledge to speed things up,
				// but it happens so rarely that we don't bother
				hc_poly_twist (&twist, c);
				k = jac_order (&PN1, Min, Max, JAC_INVALID_A1, 0, 0, 0, &twist);
				if ( k <= 0 ) { printf ("%lu: Attempted twist order computation failed\n", p);  return -2; }
				if ( k > 1 ) { printf ("%lu: Ambiguous result in genus 2 not handled\n", p);  return -2; }		// should be impossible
				if ( ! smalljac_genus2_charpoly (a, p, x, P1, PN1) ) { err_printf ("%lu: smalljac_genus2_charpoly failed\n", p);  return -2; }
				return 2;
			}
			return 2;
		}
	}
	if ( c->g == 3 ) {
		if ( k == 1 ) {
			if ( (flags&SMALLJAC_GROUP) ) return jac_structure (a, c, P1, flags&SMALLJAC_SGROUP);
			return ( smalljac_genus3_charpoly_from_P1 (a, P1, pts, Min, Max, c) ? 3 : -2 );
		} else {
			unsigned long unique_P1, unique_PN1;
			
			hc_poly_twist (&twist, c);
			m = 2*(p+1);
			e = P1;
			P1 = _ui_ceil_ratio(Min,e)*e;
			for ( i = 0 ; i < k ; i++ ) {
				tcon = 2*(p*p*p+1) - P1;
				constraints[i] = (tcon < 0 ? (int) (m - ((-tcon)%m)) : (int) (tcon%m) );
				P1 += e;
			}
			constraints[k] = -1;
			tk = jac_order (&PN1, Min, Max, (pts ? -a1 : JAC_INVALID_A1), 0, 0, constraints, &twist);
			if ( tk <= 0 ) { printf ("%lu: Attempted twist order computation failed\n", p);  return -2; }
			P1 = _ui_ceil_ratio(Min,e)*e;
			te = PN1;
			n = 0;
			unique_P1 = unique_PN1 = 0;
			for ( i = 0 ; i < k ; i++ ) {
				PN1 = _ui_ceil_ratio(Min,te)*te;
				for ( j = 0 ; j < tk ; j++ ) {
					sts = smalljac_genus3_charpoly (a, p, x, P1, PN1, pts);
					if ( sts ) {
						if ( ! n ) { unique_P1 = P1;  unique_PN1 = PN1; }
						n++;
					}
					PN1 += te;
				}
				P1 += e;
			}
			if ( n != 1 ) { printf ("%lu: Unable to resolve ambiguity, primary subgroup %lu, twist subgroup %lu, pts = %lu, n = %d, giving up\n", p, e, te, pts, n);  return -2; }
			if ( (flags&SMALLJAC_GROUP) ) { n = jac_structure (a, c, unique_P1, 0); return ( n > 0 ? n : -2 ); }
			return ( smalljac_genus3_charpoly (a, p, x, unique_P1, unique_PN1, pts) ? 3 : -2 );
		}
	}
	err_printf ("Fell through to unhandled case in smalljac_generic_Lpoly with curve of genus %d\n", c->g);  return -2;
}

/*
// only compile in genus > 2

// Computes Lpoly coefficients or group structure using p-adic computations to determine L_p(T) mod p
// by calling David Harvey's frobenius() function.  If group structure is requested, uses #J(C/F) = L_p(1)
// to compute the group order.  Returns 0 if an error occurs, the number of entries in a otherwise (either g or the rank of the group)
// curve is assumed to be non-singular and defined over Q
int smalljac_padic_Lpoly (long a[], smalljac_curve *sc, long p, unsigned long flags)
{
	hc_poly c;
	long f[SMALLJAC_MAX_DEGREE+1], n;
	int i;
	
	assert (sc->Qflag);
	ff_setup_ui (p);
	hc_poly_set_mpz (sc->hc, sc->f, sc->degree);

	// this is slightly annoying, but we need to convert to standard integers (mod p) and ff uses Montogemery
	for ( i = 0 ; i <= sc->degree ; i++ ) f[i] = _ff_get_ui (c.f[i]);
	if ( ! padic_charpoly (a, f, sc->degree, p) ) { err_printf ("%lu: padic_charpoly failed\n", p);  return -2; }
	switch ( sc->genus) {
	case 2:
		if ( ! smalljac_genus2_charpoly_from_Pmodp (a, sc->hc) ) return -2;
		break;
	case 3:
		if ( ! smalljac_genus3_charpoly_from_Pmodp (a, sc->hc) ) return -2;
		break;
	default:
		err_printf ("Invalid genus %d in smalljac_padic_lpoly\n", sc->genus);
		exit (0);
	}
	if ( (flags&SMALLJAC_GROUP) ) {
		n = jac_structure (a, sc->hc, smalljac_Lp1_ui (a, sc->genus, p), flags&SMALLJAC_SGROUP);
		return ( n ? n : -2 );
	}
	return sc->genus;
}
*/
smalljac_curve *smalljac_curve_alloc (void)
{
	smalljac_curve *sc;
	
	sc = mem_alloc (sizeof(*sc));
	mpz_init (sc->disc);
	mpz_init (sc->D);
	return sc;
}

smalljac_curve_t smalljac_curve_init (char *str, int *err)
{
	smalljac_curve *sc;

	sc = smalljac_curve_alloc();
	if ( ! smalljac_curve_set_str (sc, str, err) ) { smalljac_curve_clear ((smalljac_curve_t)sc);  return (smalljac_curve_t)0; }
	return (smalljac_curve_t) sc;
}


void smalljac_curve_check_special (smalljac_curve *sc)
{
	mpq_t I[3];
	int i;

	if ( sc->nfd != 1 ) return;
	switch (sc->degree) {
	case 3:
	case 4:
		if ( sc->dF == 3 && mpz_cmp_ui(sc->F[3],1)==0 && sc->dH <= 1 ) { 	// only curves over Q in Weierstrass form are checked
			for ( i = sc->dH+1 ; i <= 1 ; i++ ) mpz_set_ui (sc->H[i], 0);
			if ( mpz_ws_has_cm (sc->H[1], sc->F[2], sc->H[0], sc->F[1], sc->F[0]) ) sc->special = SMALLJAC_SPECIAL_CM;
		}
		break;
	case 5:
	case 6:
		if ( sc->degree == 5 && mpz_cmp_ui(sc->f[5],1)==0 && ! mpz_sgn(sc->f[4]) && ! mpz_sgn(sc->f[3]) && ! mpz_sgn(sc->f[2]) ) {
			if ( ! mpz_sgn (sc->f[0]) ) sc->special = SMALLJAC_SPECIAL_X5PAX;	// y^2 = x^5 +ax
			if ( ! mpz_sgn (sc->f[1]) ) sc->special = SMALLJAC_SPECIAL_X5PA;	// y^2 = x^5 + a
		}
		if ( sc->degree == 6 && mpz_cmp_ui(sc->f[6],1)==0 ) {
			for ( i = 1 ; i < 6 ; i++ ) if ( mpz_sgn(sc->f[i]) ) break;
			if ( i == 6 ) sc->special = SMALLJAC_SPECIAL_X6PA;				// y^2 = x^6 + a
			// check for the very special curve y^2=x^6-5x^4-5x^2+1, whose Jacobian is Q-isogenous to the square of the elliptic curve y^2=x^3-5x^2-5x+1 with CM by the order with discriminant -2.
			if ( mpz_cmpabs_ui(sc->f[2],5)==0 && mpz_sgn(sc->f[2])<0 && !mpz_sgn(sc->f[3]) && mpz_cmpabs_ui(sc->f[4],5)==0 && mpz_sgn(sc->f[4])<0 && ! mpz_sgn(sc->f[5]) ) sc->special = SMALLJAC_SPECIAL_CM2SQUARE;
		}
		if ( ! sc->special ) {
			// check for twist of y^2=x^5+x with igusa invariants (400000,-20000,-2000)
			// and for twist of y^2=x^6+1 with igusa invariants (51200000/3,480000,148000)
			for ( i = 0 ; i < 3; i++ ) mpq_init (I[i]);
			mpq_poly_igusa_inv (I, sc->f, sc->degree);
			if ( mpq_x5px_twist (I) ) {
				sc->special = SMALLJAC_SPECIAL_X5PX_TWIST;
				info_printf ("Twist of y^2=x^5+x\n");
			} else if ( mpq_x6p1_twist (I) ) {
				sc->special = SMALLJAC_SPECIAL_X6P1_TWIST;
				info_printf ("Twist of y^2=x^6+1\n");
			}
			for ( i = 0 ; i < 3; i++ ) mpq_clear (I[i]);
		}
		break;
	case 7:
		if ( mpz_cmp_ui(sc->f[7],1)==0 && ! mpz_sgn(sc->f[6]) && ! mpz_sgn(sc->f[5]) && ! mpz_sgn(sc->f[4]) && ! mpz_sgn(sc->f[3]) && ! mpz_sgn(sc->f[2]) ) {
			if ( ! mpz_sgn (sc->f[0]) ) sc->special = SMALLJAC_SPECIAL_X7PAX;	// y^2 = x^7 +ax
		}
		break;
	case 8:
		if ( mpz_cmp_ui(sc->f[8],1)==0 ) {
			for ( i = 1 ; i < 8 ; i++ ) if ( mpz_sgn(sc->f[i]) ) break;
			if ( i == 8 ) sc->special = SMALLJAC_SPECIAL_X8PA;				// y^2 = x^8 + a
		}		
		break;
	}
}

// Note, this will not handle zero leading coefficients
int inline quick_poly_degree_check (char *str, char x)
{
	register char *s;
	register int i, d;
	
	d = 0;  s = str;
	while ( *str && *str != ',' ) {
		while ( *s && *str != ',' && *s != x ) s++;
		if ( ! *s ) break;
		for ( s++ ; isspace(*s) ; s++ );
		if ( *s == '^' ) {
			for ( s++ ; isspace(*s) ; s++ );
			i = atoi(s);
			if ( i > d ) d = i;
		} else {
			if ( 1 > d ) d = 1;
		}
	}
	return d;
}


// attempt to determine the type, degree, and genus of a specified curve expression, and whether it is defined over Q or a number field
// Three types of expressions are currently supported:
//    (1) Weierstrass coordinates [a1,a2,a3,a4,a6] of an elliptic curve
//    (2) univariate polynomial f(x) specifying an elliptic or hyperelliptic curve y^2=f(x) or a Picard curve y^3=f(x) if deg f = 4
//    (3) bivariate (or trivariate) polynomial specifying a plane quartic
int smalljac_analyze_curve_string (char *str, int *type, int *degree, int *genus, int *Qflag, int *ws, int *nf)
{
	register char *s;
	char x, z;
	int i, d;

	while ( isspace (*str) ) str++;
	*Qflag =1;  *ws = *nf = 0;
	if ( str[0]=='[' && str[1]==']' ) return 0;
	z = 0;
	*nf = 0;
	if ( (s=strchr (str,'/')) ) {
		for ( s++ ; isspace(*s) ; s++ );
		if ( *s == '(' ) {
			for ( s++ ; *s && ! isalpha(*s) ; s++ );
			if ( *s ) { z = *s; *nf = 1; }
		}
	}
	if ( z ) {  for ( s = str ; *s && *s != '/' ; s++ ) if ( *s == z ) *Qflag = 0; }
	if ( str[0] == '['  && strchr (str, ',') ) {
		// there are two possibilities here, we could have Weierstrass coords [a1,a2,a3,a4,a6] or two polys [f(x),h(x)] specifying the hyperelliptic curve y^2+h(x)*y = f(x)
		for ( i = 0, s = str ; (s = strchr(s+1, ',')) ; i++ );
		if ( i == 4 ) {
			*type = SMALLJAC_CURVE_ELLIPTIC;  *degree = 3;  *genus = 1;  *ws = 1;
			return 1;
		} else if ( i == 1 ) {
			int df, dh;
			for ( s = str ; *s && *s != '(' && ! (isalpha(*s) && *s != z) ; s++ );
			if ( ! isalpha(*s) ) return 0;
			x = *s;
			df = quick_poly_degree_check (str, x);
			dh = quick_poly_degree_check (strchr(str,',')+1, x);
			*type = SMALLJAC_CURVE_HYPERELLIPTIC;
			if ( 2*dh > df ) *degree = 2*dh; else *degree = df;
			*genus = (*degree-1)/2;
			return 1;
		} else {
			return 0;
		}
	}
	for ( s = str ; *s && *s != '(' && ! (isalpha(*s) && *s != z) ; s++ );
	if ( ! isalpha(*s) ) return 0;
	x = *s;
	for ( s = str ; *s && *s != '(' && ! (isalpha(*s) && *s != x && *s != z) ; s++ );
	if ( isalpha(*s) ) {
		// assume it's a plane quartic if more than one variable is specified
		*type = SMALLJAC_CURVE_PLANE_QUARTIC;  *degree = 4;  *genus = 3;
	} else {
		d = quick_poly_degree_check (str, x);
		if ( d < 3 ) return 0;
		if ( d == 3 ) { *type = SMALLJAC_CURVE_ELLIPTIC;  *degree = 3;  *genus = 1; }
		else if ( d == 4 ) { *type = SMALLJAC_CURVE_PICARD;  *degree = 4;  *genus = 3; }
		else { *type = SMALLJAC_CURVE_HYPERELLIPTIC;  *degree = d;  *genus = (d-1)/2; }
	}
	return 1;
}

void smalljac_curve_cleanup_nf (smalljac_curve *sc)
{
	if ( sc->nfp ) { nf_poly_clear (sc->nfp);  mem_free (sc->nfp); }
	sc->nfp = 0;
}

int smalljac_nf_curve_is_singular (smalljac_curve *sc)
{
	long p;

	if ( ! sc->type == SMALLJAC_CURVE_ELLIPTIC && ! sc->type == SMALLJAC_CURVE_HYPERELLIPTIC ) return 0;	// only deal with elliptic/hyperelliptic curves
	for ( register int i = 2 ; i < 1000 ; i++ ) {
		p = ui_small_prime (i);
		if ( nf_poly_reduce_setup (sc->nfp, p, 1) ) {
			if ( ! smalljac_reduce_nf_poly (sc, p, 1, 0) ) continue;
			if ( ff_poly_discriminant_nonzero (sc->hc->f, sc->hc->d) ) return 0;	// curve is definitely not singular
		}
	}
	return 1;    	// if we cannot get good reduction at any of the first 1000 primes, assume curve is singular
}

// Given a string specifying the curve, parses the curve coefficients and dermines the genus, sets up the field of definition (Q or a number field), and if defined over Q, determines the discriminant (used to detect bad primes)
// Also checks for various special cases
int smalljac_curve_set_str (smalljac_curve *sc, char *str, int *err)
{
	int i, ws, nf;
	char *s;
	
	if ( ! _smalljac_initted ) smalljac_init();
	if ( err ) *err = 0;
	smalljac_curve_cleanup_nf (sc);						// cleanup any left-over number field poly if sc is being reused
	if ( strlen(str)+1 > sizeof(sc->str) ) { if ( err ) *err = SMALLJAC_PARSE_ERROR;  goto error; }
	strcpy (sc->str, str);
	sc->flags = sc->special = sc->ws2 = sc->ws3 = 0;  sc->nfd = 1;
	if ( ! smalljac_analyze_curve_string (sc->str, &sc->type, &sc->degree, &sc->genus, &sc->Qflag, &ws, &nf) ) { if ( err ) *err = SMALLJAC_PARSE_ERROR;  goto error; }
	if ( sc->genus < 1 || sc->genus > SMALLJAC_MAX_GENUS || sc->degree > SMALLJAC_MAX_DEGREE ) { if ( err ) *err = SMALLJAC_UNSUPPORTED_CURVE;  goto error; }
	smalljac_curve_init_f (sc);
	mpz_set_ui (sc->D,1);
	
	if ( sc->type == SMALLJAC_CURVE_PICARD ) {
		if ( nf ) { printf ("plane quartics currently only supported over Q (but you can use prime filtering (e.g. in lpdata or moments) to effectively extend the field of definition)\n");  if ( err ) *err = SMALLJAC_PARSE_ERROR;  goto error; }
		mpz_poly_parse (sc->f, sc->degree, str);
		mpz_fast_poly_discriminant (sc->disc, sc->f, sc->degree);
		mpz_mul_ui (sc->D, sc->disc, 3);
		return 1;
	}

	// code to handle twists of Klein and Fermat plane quartics -- all special cases, only supported over Q
	if ( sc->type == SMALLJAC_CURVE_PLANE_QUARTIC ) {
		if ( nf ) { printf ("plane quartics currently only supported over Q (but you can use moments filtering to effectively extend the field of definition)\n");  if ( err ) *err = SMALLJAC_PARSE_ERROR;  goto error; }
		while ( sc->f_inits < 15 ) mpz_init (sc->f[sc->f_inits++]);
		if ( ! mpz_poly_parse_plane_quartic (sc->f, str) ) { if ( err ) *err = SMALLJAC_PARSE_ERROR;  goto error; }
		sc->special = SMALLJAC_SPECIAL_FK_TWIST;
		sc->special_curve_id = smalljac_FKtwist_lookup (sc->f);
		if ( ! sc->special_curve_id ) { printf ("{ "); for ( i = 0 ; i < 15 ; i++ ) gmp_printf ("%Zd, ", sc->f[i]); puts ("}\n"); if ( err ) *err = SMALLJAC_UNSUPPORTED_CURVE;  goto error; }
		// we make no attempt to compute the discriminant of the curve
		return 1;
	}
	
	assert ( sc->type == SMALLJAC_CURVE_ELLIPTIC || sc->type == SMALLJAC_CURVE_HYPERELLIPTIC );

	if ( nf ) {
		if ( ! sc->nfp ) sc->nfp = mem_alloc (sizeof(*sc->nfp));
		if ( nf_poly_init (sc->nfp, str) != sc->degree ) { if ( err ) *err = SMALLJAC_PARSE_ERROR; goto error; }
		sc->nfd = sc->nfp->n;
		for ( s = sc->str ; *s != '/' ; s++ ); 
		*s++ = '\0';
		while ( isspace (*s) ) s++;
		sc->nfstr = s;
		if ( smalljac_nf_curve_is_singular (sc) ) { if ( err ) *err = SMALLJAC_SINGULAR_CURVE; goto error; }
		if ( ! sc->nfp->Qflag && !ws && strchr(sc->str, ',') ) { err_printf ("Currently hypelliptic curves of the form y^2+h(x)y=f(x) must be defined over Q, put curve in the form y^2=f(x)\n");  if ( err ) *err = SMALLJAC_PARSE_ERROR; goto error; }
		// for curves not defined over Q we deal with all reduction issues on a per prime basis by reducing in nf_poly, so we are a done here
		if ( ! sc->nfp->Qflag ) return 1;
	}

	// if we get here we have an elliptic or hyperelliptic curve defined over Q
	smalljac_curve_init_FH (sc);	// the "official" model of the curve is Y^2 + H(X)*Y = F(X), we use an isomorphic model y^2 = f(x) at all primes p > 2g+2

	if ( ws ) {
		 mpz_set_ui (sc->F[3], 1);
		if ( gmp_sscanf (str, "[%Zd,%Zd,%Zd,%Zd,%Zd]", sc->H[1], sc->F[2], sc->H[0], sc->F[1], sc->F[0]) != 5 ) { if ( err ) *err = SMALLJAC_PARSE_ERROR; goto error; }
		sc->dF = 3; 
		sc->dH = mpz_sgn(sc->H[1]) ? 1 : ( mpz_sgn(sc->H[0]) ? 0 : -1 );
		mpz_short_ws (sc->f[1], sc->f[0], sc->H[1], sc->F[2], sc->H[0], sc->F[1], sc->F[0]);  mpz_set_ui (sc->f[3], 1); mpz_set_ui (sc->f[2], 0);
		sc->flags |= SMALLJAC_CURVE_FLAG_WS;
	} else {
		if ( mpz_hyperelliptic_curve_parse (sc->F, &sc->dF, sc->H, &sc->dH, sc->degree, sc->str) != sc->genus )  { if ( err ) *err = SMALLJAC_PARSE_ERROR; goto error; }
	}
	if ( mpz_hyperelliptic_curve_normalize (sc->f, sc->F, sc->dF, sc->H, sc->dH) != sc->degree )  { if ( err ) *err = SMALLJAC_INTERNAL_ERROR; goto error; }
	mpz_hyperelliptic_curve_discriminant (sc->D, sc->F, sc->dF, sc->H, sc->dH);
	if ( ! mpz_sgn (sc->D) )  {if ( err ) *err = SMALLJAC_SINGULAR_CURVE; goto error; }
	mpz_fast_poly_discriminant (sc->disc, sc->f, sc->degree);

	// check for curves of special form (e.g. y^2=x^6+a, y^2=x^5+ax, y^3=quartic(x), etc...) -- note this must happen last)
	smalljac_curve_check_special (sc);

	return 1;

error:
	if ( err ) assert (*err);
	smalljac_curve_cleanup_nf (sc);
	return 0;
}

int smalljac_curve_cm (smalljac_curve_t curve)
{
	smalljac_curve *sc = curve;
    assert (sc->genus == 1 && sc->nfd == 1 && sc->dH <= 1 && sc->dF == 3 && mpz_cmp_ui(sc->F[3],1)==0); // complain if curve is not over Q and in Weierstrass form
	return (sc->special&SMALLJAC_SPECIAL_CM) ? 1 : 0;
}
char *smalljac_curve_str (smalljac_curve_t curve) { return ((smalljac_curve*)curve)->str; }
char *smalljac_curve_nf (smalljac_curve_t curve) { return ((smalljac_curve*)curve)->nfstr; }
int smalljac_curve_nf_degree (smalljac_curve_t curve) { return ((smalljac_curve*)curve)->nfd; }
int smalljac_curve_genus (smalljac_curve_t curve) { return ((smalljac_curve*)curve)->genus; }

// creates a smalljac elliptic/hyperelliptic curve y^2=f(x) with integer coefficients defined over Q
int smalljac_curve_set_mpz (smalljac_curve *sc, mpz_t f[], int degree,  char *str)
{
	int i;
	
	if ( degree < 3 || degree > SMALLJAC_MAX_DEGREE || degree == 4 ) return 0;
	sc->flags = sc->special = sc->ws2 = sc->ws3 = 0;
	sc->Qflag = 1; sc->nfd = 1;
	sc->degree = degree;  sc->genus = (degree-1)/2;
	sc->type = ( sc->genus == 1 ? SMALLJAC_CURVE_ELLIPTIC : SMALLJAC_CURVE_HYPERELLIPTIC );
	smalljac_curve_init_f (sc);
	smalljac_curve_init_FH (sc);
	for ( i = 0 ; i <= sc->degree ; i++ ) mpz_set (sc->f[i], f[i]);
	mpz_fast_poly_discriminant (sc->disc, sc->f, sc->degree);
	mpz_mul_2exp (sc->D, sc->disc, 4*sc->genus);																	// get the right power of 2
	if ( (sc->degree&1) ) { mpz_mul (sc->D, sc->D, sc->f[sc->degree]); mpz_mul (sc->D, sc->D, sc->f[sc->degree]); }	// include lc(f)^2 in D if deg(f) is odd
	smalljac_curve_check_special (sc);
	if ( str ) strcpy (sc->str, str); else sc->str[0] = '\0';
	return 1;
}

// creates a smalljac elliptic/hyperelliptic curve with integer coefficients defined over Q
int smalljac_curve_set_i (smalljac_curve *sc, long f[], int degree, char *str)
{
	int i;
	
	if ( degree < 3 || degree > SMALLJAC_MAX_DEGREE || (! (degree&1)&&degree!=6)  ) return 0;
	sc->flags = sc->special = sc->ws2 = sc->ws3 = 0;
	sc->Qflag = 1; sc->nfd = 1;
	sc->degree = degree;  sc->genus = (degree-1)/2;
	sc->type = ( sc->genus == 1 ? SMALLJAC_CURVE_ELLIPTIC : SMALLJAC_CURVE_HYPERELLIPTIC );
	smalljac_curve_init_f (sc);
	for ( i = 0 ; i <= degree ; i++ ) mpz_set_i (sc->f[i], f[i]);
	mpz_fast_poly_discriminant (sc->disc, sc->f, degree);
	mpz_mul_2exp (sc->D, sc->disc, 4*sc->genus);																	// get the right power of 2
	if ( (sc->degree&1) ) { mpz_mul (sc->D, sc->D, sc->f[sc->degree]); mpz_mul (sc->D, sc->D, sc->f[sc->degree]); }	// include lc(f)^2 in D if deg(f) is odd
	smalljac_curve_check_special (sc);
	if ( str ) strcpy (sc->str, str); else sc->str[0] = '\0';
	return 1;
}

void smalljac_curve_clear (smalljac_curve_t curve)
{
	smalljac_curve *sc;
	int i;
	
	sc = (smalljac_curve *)curve;
	smalljac_curve_cleanup_nf (sc); 
	mpz_clear (sc->D);
	mpz_clear (sc->disc);
	for ( i = 0 ; i < sc->f_inits ; i++ ) mpz_clear (sc->f[i]);
	for ( i = 0 ; i < sc->Delta_inits ; i++ ) mpz_clear (sc->Deltas[i]);
	for ( i = 0 ; i < sc->F_inits ; i++ ) mpz_clear (sc->F[i]);
	for ( i = 0 ; i < sc->H_inits ; i++ ) mpz_clear (sc->H[i]);
	mem_free (sc);
}

// given Lpoly at p for a curve defined over q, compute the Lpoly coefficients at q=p^h
// currently only handles the trace
int smalljac_Lpoly_extend (long a[], int n, long p, int h)
{
	register long t, t1, tn, tnm1;
	register int i;

	if ( h < 1 || n < 1 ) return 0;
	if ( h == 1 ) return 1;
	switch (n) {
	case 1:
		tnm1 = t1 = -a[0];  tn = t1*t1-2*p; 
		for ( i = 2 ; i < h ; i++ ) { t = t1*tn - p*tnm1;  tnm1 = tn;  tn = t; }
		a[0] = -tn;
		break;
	default:
		err_printf ("not implemented error in smalljac_Lpoly_extend\n"); exit (0);
	}
	return 1;
}
