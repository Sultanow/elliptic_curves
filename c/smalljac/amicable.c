#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "hecurve.h"
#include "smalljac_internal.h"
#include "smalljac.h"
#include "cstd.h"

/*
    Copyright (c) 2010-2013 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// Program to enumerate amicable pairs and aliquot cycles for elliptic curve E/Q

#define MAX_PATH_LENGTH	100				// Length of longest aliquot path that can be followed (and cycle that can be detected).
#define MAX_THREADS		64				// Maximum number of threads.  Actual number will be set to 1,2,4,... depending on # cores available.

// callback context structure used to keep statistics
struct callback_ctx {
	unsigned long prime_count;				// counts primes where E/Fp has good reduction and #E(Fp) is a prime < p
	unsigned long pair_count;				// counts amicable pairs (p1,p2) with p1 > p2
	unsigned long cycle_count;				// counts distinct aliquot cycles of length > 2
	unsigned long cpu_musecs;				// cpu time in microseconds, per thread
	unsigned long sum;					// this will wrap in a large computation, but it is handy as a checksum
	int id;
};

void log_pair (struct callback_ctx *ctx, long p1, long p2)
{
	char filename[256];
	FILE *fp;
	
	sprintf (filename, "amicable_%d.log", ctx->id);
	fp = fopen (filename, "a");
	fprintf (fp, "(%ld,%ld)  [%ld]\n", p1, p2, ctx->pair_count);
	fclose(fp);
}

/*
	Callback function invoked by smalljac for each prime p where E has good reduction and #E(F_p) has prime order less than p.
	This is achieved by setting flags = SMALLJAC_PRIME_ORDER | SMALLJAC_LOW_ORDER | SMALLJAC_GOOD_ONLY below.
*/
int callback (smalljac_curve_t curve, unsigned long p, int good, long a[], int n, void *arg)
{
	struct callback_ctx *ctx;
	ff_t f[4];
	long N, N2, Nlist[MAX_PATH_LENGTH];
	int i,j, k;

	ctx = (struct callback_ctx*) arg;
	N = (long)p + 1L + a[0];						// note that a[0] is the L-poly coefficient and has sign opposite the trace
//printf ("callback p=%ld, a[0]=%ld, N=%ld\n", p, a[0], N);
	// we assume that since SMALLJAC_PRIME_ORDER and SMALLJAC_LOW_ORDER are set, N must be a prime < p.  We don't verify this here.
	ctx->prime_count++;
	ctx->sum += N;
	ff_setup_ui(N);
	if ( ! smalljac_Qcurve_reduce(f,curve) ) { printf ("Bad reduction at %ld, couldn't check for amicable pair\n", N);  return 1; }
/*
	// Uncomment to only test for amicable pairs, which is slightly faster than following every Aliquot sequence to termination (but the difference is small, a few percent).
 	if ( ecurve_test_exponent (p, f) ) { printf ("[%ld,%ld]\n", p, N); ctx->pair_count++; }
	return 1;
*/
	N2 = ecurve_prime_order (f,0);
	if ( N2 && N2 != N ) {
		if ( N2 == p) {
			printf ("(%ld,%ld)\n", N, p); fflush(stdout); ctx->pair_count++;
log_pair (ctx, N, p);
		} else {
			Nlist[0] = p; Nlist[1] = N;
			for ( i = 2 ; i < MAX_PATH_LENGTH ; i++ ) {
				Nlist[i] = N2;
				ff_setup_ui(N2);
				if ( ! smalljac_Qcurve_reduce(f,curve) ) { printf ("Aliquot path terminated due to bad reduction at %ld\n", N);  return 1; }
				N2 = ecurve_prime_order (f,0);
				if ( ! N2 || N2 == _ff_p ) break;		// dont't count cycles of length 1
				for ( j = 0 ; j <= i ; j++  ) if ( N2 == Nlist[j] ) break;
				if ( j <= i ) {
					if ( j < i ) {
						// Only count cycle if we detected it from its largest prime, to avoid duplicates.  This also means we ignore any leading tail.
						if ( j ) break;
						for ( j = 1 ; j <= i ; j++ ) if ( Nlist[j] > p ) break;
						if ( j <= i  ) break;
						// Print cycle starting with smallest prime, consistent with the normalization in Silverman-Stange
						for ( j = k = 0 ; j <= i ; j++ ) if ( Nlist[j] < Nlist[k] ) k = j;
						printf ("(%ld", Nlist[k]);
						for ( j = (k==i?0:k+1) ; j != k ; j = (j==i?0:j+1) ) printf (",%ld", Nlist[j]);
						printf (")\n"); fflush(stdout); 
					}
					ctx->cycle_count++;
					break;
				}
			}
			if ( i == MAX_PATH_LENGTH ) { printf ("Exceed MAX_PATH_LENGTH=%d starting from p=%ld\n", MAX_PATH_LENGTH, p); return 1; }
		}
	}
	return 1;
}


int main (int argc, char *argv[])
{
	time_t start, end;
	clock_t start_cpu;
	unsigned long startp, endp, flags;
	smalljac_curve_t curve;
	int child_rpipe[MAX_THREADS-1][2], child_wpipe[MAX_THREADS-1][2];			// we only need MAX_THREADS-1 children because the parent will also be working
	FILE *in[MAX_THREADS-1], *out[MAX_THREADS-1];
	pid_t child_pid[MAX_THREADS-1];
	struct callback_ctx context, child_context;
	long result;
	int i, k, err, threads, child;

printf ("max prime = %lu\n", MPZ_MAX_ENUM_PRIME);
	if ( argc < 4) { puts ("amicable startp endp curve [threads]");  return 0; }
	startp = atol_exp (argv[1]);
	endp = atol_exp (argv[2]);

	flags = SMALLJAC_PRIME_ORDER | SMALLJAC_LOW_ORDER | SMALLJAC_GOOD_ONLY;
	memset(&context,0,sizeof(context));

	// create the curve - give an explanation of the error if there is a problem
	curve = smalljac_curve_init (argv[3], &err);
	if ( ! curve ) { 
		switch (err) {
		case SMALLJAC_PARSE_ERROR: printf ("Unable to parse curve string: %s\n", argv[3]); break;
		case SMALLJAC_UNSUPPORTED_CURVE: puts ("Specified curve not supported\n");  break;
		case SMALLJAC_SINGULAR_CURVE: puts ("Specified curve is singular\n");  break;
		default: printf ("smalljac_curve_init returned error %d\n", err);
		}
		return 0;
	}

	if ( argc > 4 ) {
		threads = atoi(argv[4]);
	} else {
		threads = sysconf(_SC_NPROCESSORS_ONLN) ;
		if ( threads <= 0 ) { printf ("Couldn't determine onliine processor count, assuming just one.\n"); threads = 1; }
	}
	if ( threads > MAX_THREADS ) threads = MAX_THREADS;
	k = ui_len(threads)-1;
	threads = 1<<k;		// round threads down to nearest power of 2
	flags |= ((1UL<<k)-1) << (SMALLJAC_SPLIT_SHIFT+1);

	fflush(stdout);
	child = 0;
	for ( i = 0 ; i < threads-1 ; i++ ) {
		if  ( pipe (child_rpipe[i]) == -1 || pipe(child_wpipe[i]) == -1 ) { printf ("Error creating pipe: %d\n", errno);  return 0; }
		child_pid[i] = fork();
		if ( child_pid[i] ) {  // in parent
			close(child_rpipe[i][0]); out[i] = fdopen (child_rpipe[i][1], "w");  in[i] = fdopen (child_wpipe[i][0], "r"); close(child_wpipe[i][1]);
		} else {  // in child
			in[i] = fdopen (child_rpipe[i][0], "r"); close (child_rpipe[i][1]);  close (child_wpipe[i][0]);  out[i] = fdopen (child_wpipe[i][1], "w");
			flags |= (i+1) << (SMALLJAC_HIGH_SHIFT+1);		// parent thread handles 0 case
			child = 1; break;
		}
	}
	context.id = i;
	if ( ! child ) printf ("Beginning search of interval [%ld,%ld] for curve %s using %d threads...\n", startp, endp, argv[3], threads);
	start_cpu = clock(); start = time(0);
	result = smalljac_Lpolys (curve, startp, endp, flags, callback, (void*)&context);					// this is where it all happens
	if ( result < 0 ) {  printf ("smalljac_Lpolys returned error %ld\n", result);  return 0; }
	smalljac_curve_clear (curve);
	context.cpu_musecs = clock()-start_cpu;  end = time(0);
	
	if ( child ) { fwrite (&context, sizeof(context), 1, out[i]);  fclose (out[i]);  return 0; }
	
	for ( i = 0 ; i < threads-1 ; i++ ) {
		if ( ! fread(&child_context, sizeof(child_context), 1, in[i]) ) { printf ("error reading child context\n"); return 0; }
		fclose (in[i]);
		context.prime_count += child_context.prime_count;
		context.pair_count += child_context.pair_count;
		context.cycle_count += child_context.cycle_count;
		context.cpu_musecs += child_context.cpu_musecs;
		context.sum += child_context.sum;
	}

	printf ("...scan completed.\n");
	printf ("Examined %lu prime order curves with #E(Fp) < p, finding %lu amicable pairs and %lu aliquot cycles.\n", context.prime_count, context.pair_count, context.cycle_count);
	printf ("Cpu time: %.1f cpu seconds,  Elapsed time: %ld seconds\n", (double)context.cpu_musecs/1000000.0, end-start);
	printf ("Sum of prime N_p < p for p in [%lu,%lu] is %lu (mod 2^64)\n", startp, endp, context.sum);
}
