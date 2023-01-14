#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "smalljac.h"
#include "smalljac_internal.h"

/*
    Copyright (c) 2012 Pavel Panchekha and Andrew V. Sutherland
    See LICENSE file for license details.
*/

#define MOMENTS		11
#define MAXM			4096
#define MAX_THREADS	256	// Maximum number of threads.  Actual number will be set based on # cores available.
#define A2MAX		((2*SMALLJAC_MAX_GENUS*(2*SMALLJAC_MAX_GENUS-1))/2)
#define SPLIT_POLY_MAX_DEGREE	100

typedef int (*smalljac_filter_callback)(smalljac_curve_t curve, unsigned long q, void *arg);

struct smalljac_moments_ctx {
	unsigned long ap_total;
	unsigned long tot_count;
	unsigned long good_count;
	unsigned long bad_count;
	unsigned long zero_count;
	unsigned long a2_icounts[2*A2MAX+1];	// counts of integer a2 values,  entry i is incremented whenever a2 is exactly i-A2MAX

	double moments[MOMENTS];
	double a2moments[MOMENTS];
	double a3moments[MOMENTS];

	smalljac_filter_callback filter_callback;
	void *filter_arg;
};

void smalljac_moments_ctx_init (struct smalljac_moments_ctx *context)
{
  memset(context, 0, sizeof(*context));
}

static int repeat;

int smalljac_moments_callback (smalljac_curve_t curve, unsigned long p, int good, long a[], int n, void *arg)
{
	struct smalljac_moments_ctx *ctx;
	double x, y, z;
	int i;
	long k;

	ctx = (struct smalljac_moments_ctx*) arg;
	if ( good < 0 ) { // negative good value indicates filtering
		// 2 is always a bad prime
		if ( p==2 ) goto bad;

		// Leading coefficient might be divisible by p
		if ( smalljac_curve_lc_divisible_p (curve, p) ) goto bad;

		// Must have an actual callback
		if ( ! ctx->filter_callback ) {
			printf ("Don't know how to filter\n");
			abort();
		}

		// Store how many times to count this prime
		repeat = ctx->filter_callback(curve, p, ctx->filter_arg);
		return repeat;

	bad:
		ctx->bad_count++;
		return 0;
	}

	// Only count after filtering out primes we ignore
	ctx->tot_count += repeat;

	if ( ! good ) { ctx->bad_count++; return 1; }

	if ( ! ctx->filter_callback ) repeat = 1;
	x = sqrt(p);
	for ( ; repeat ; repeat-- ) { // repeat is a static variable
		ctx->good_count++;
		ctx->ap_total -= a[0];
		z = a[0]/x;
		if ( ! a[0] ) ctx->zero_count++;
		for ( y = z, i = 1 ; i < MOMENTS ; y *= z, i++ ) ctx->moments[i] += y;
		if ( n > 1 ) {
			z =(double)a[1]/p;
			for ( y = z, i = 1 ; i < MOMENTS ; y *= z, i++ ) ctx->a2moments[i] += y;
			if ( n > 2 ) {
				z =(double)a[2]/(p*x);
				for ( y = z, i = 1 ; i < MOMENTS ; y *= z, i++ ) ctx->a3moments[i] += y;
			}
			k = a[1]/(long)p;
			if ( k >= -2 && k <= 2 && k*p == a[1] ) ctx->a2_icounts[A2MAX+k]++;
		}
	}
	return 1;
}

int smalljac_moments (smalljac_curve_t curve, unsigned long start, unsigned long end, double moments[], int n, int m, char STgroup[16], smalljac_filter_callback filter_callback, void *arg)
{
	struct smalljac_moments_ctx ctx;
	int genus, r, i;
	unsigned long flags;
	double z1, z2[5], m1sq[3], m2[4];

	genus = smalljac_curve_genus(curve);
	if ( genus < n ) { err_printf ("Invalid n=%d > genus=%d in smalljac_moments\n", n, genus); return SMALLJAC_INVALID_FLAGS; }

	flags = 0;
	flags |= SMALLJAC_DEGREE1_ONLY;

	if ( n == 1 ) flags |= SMALLJAC_A1_ONLY;
	
	smalljac_moments_ctx_init(&ctx);

	if ( filter_callback ) {
		ctx.filter_callback = filter_callback;
		ctx.filter_arg = arg;
		flags |= SMALLJAC_FILTER;
	}
	
	//smalljac_moments_compute(curve, start, end, flags, &ctx);
	smalljac_parallel_Lpolys(curve, start, end, flags, smalljac_moments_callback, &ctx);

	if ( STgroup ) {
		if ( genus == 2 && n == 2 ) {
			// Compute signature
			z1 = (double)ctx.zero_count/ctx.good_count;
			for ( i = 0 ; i < 5 ; i++ ) {
				z2[i] = (double)ctx.a2_icounts[A2MAX-2+i]/ctx.good_count;
			}

			m1sq[1] = ctx.moments[2]/ctx.good_count;
			m1sq[2] = ctx.moments[4]/ctx.good_count;
			m2[1] = ctx.a2moments[1]/ctx.good_count;
			m2[2] = ctx.a2moments[2]/ctx.good_count;
			m2[3] = ctx.a2moments[3]/ctx.good_count;

			// Look up signature
			r = smalljac_lookup_g2_STgroup (STgroup, z1, z2, m1sq, m2);
			if ( ! r ) STgroup[0] = '\0';
		} else if ( genus == 1 ) {
			// Compute signature
			z1 = (double)ctx.zero_count/ctx.good_count;
			m1sq[1] = ctx.moments[2]/ctx.good_count;
			m1sq[2] = ctx.moments[4]/ctx.good_count;

			// Look up signature
			r = smalljac_lookup_g1_STgroup (STgroup, z1, m1sq);
			if ( ! r ) STgroup[0] = '\0';
		}
	}
	if ( ! moments || ! m || ! n ) return 0;
	moments[0] = 1.0;
	for ( int i = 1; i < m; i++ ) {
		moments[i] = ctx.moments[i] / ctx.good_count;
	}
	if ( 2 <= n) {
		moments[m] = 1.0;
		for ( int i = 1; i < m; i++) {
			moments[i + m] = ctx.a2moments[i] / ctx.good_count;
		}
	}
	if ( 3 <= n) {
		moments[2*m] = 1.0;
		for ( int i = 1; i < m; i++) {
			moments[i + 2*m] = ctx.a3moments[i] / ctx.good_count;
		}
	}
    return 0;
}
