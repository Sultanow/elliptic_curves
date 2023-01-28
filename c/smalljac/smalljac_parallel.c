#include <stdio.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <string.h>
#include <errno.h>
#include "cstd.h"
#include "smalljac.h"
#include "ff_poly.h"

/*
    Copyright (c) 2012 Pavel Panchekha and Andrew V. Sutherland
    See LICENSE file for license details.
*/

#define MAX_THREADS	1	// Maximum number of threads.  Actual number will be set based on # cores available.
int smalljac_parallel_threads ()
{
	int k, threads;
	
	threads = sysconf(_SC_NPROCESSORS_ONLN);
	if ( threads <= 0 ) {
		printf ("Couldn't determine online processor count, assuming just one.\n");
		threads = 1;
	} else if ( threads > MAX_THREADS ) {
		threads = MAX_THREADS;
	}

	k = ui_len(threads) - 1;
	threads = 1 << k; // round down to nearest power of 2
	return threads;
}

int smalljac_parallel_callback (smalljac_curve_t curve, unsigned long p, int good, long a[], int n, void *arg)
{
	FILE *out = (FILE *) arg;
	fwrite(&p, sizeof(p), 1, out);
	fwrite(&good, sizeof(good), 1, out);
	fwrite(&n, sizeof(n), 1, out);
	if ( n ) fwrite(&(a[0]), sizeof(a[0]), n, out);
	return 1;
}

int smalljac_parallel_argmin(unsigned long *p, int n) {
	unsigned long min_val, min_arg;
	int i;
	min_val = (0 - 1);
	min_arg = -1;
	for ( i = 0 ; i < n ; i++ ) {
		if (p[i] < min_val) {
			min_val = p[i];
			min_arg = i;
		}
	}
	return min_arg;
}

long smalljac_parallel_Lpolys (smalljac_curve_t curve, unsigned long start, unsigned long end, unsigned long flags, int (*callback)(smalljac_curve_t, unsigned long, int, long[], int, void *), void *arg)
{
	int child_rpipe[MAX_THREADS][2], child_wpipe[MAX_THREADS][2];
	FILE *in[MAX_THREADS], *out[MAX_THREADS];
	pid_t child_pid[MAX_THREADS];
	int threads, child;
	long result;
	int i,status;

	unsigned long stream_heads[MAX_THREADS];
	// The p are sent in order, so we only need to remember the
	// filtration status of the last prime
	unsigned long p, last_filter; 							  
	int good, n;
	long a[2*SMALLJAC_MAX_GENUS];	// Modified 12/1/12 by AVS to use fixed array for a -- n should never exceed 2*MAX_GENUS

	
	last_filter = 0; // Start with everything unfiltered

	threads = smalljac_parallel_threads();

	if ( 1 == threads || end - start < 25 ) { // 25 was experimentally determined
		// Fast-path a one-thread case
		return smalljac_Lpolys(curve, start, end, flags, callback, arg);
	}

	flags |= (threads - 1) << (SMALLJAC_SPLIT_SHIFT + 1);

	fflush(0);	// clear i/o buffers before we fork anything
	
	// we only need threads children, because the parent will be aggregating
	child = 0;
	for ( i = 0 ; i < threads ; i++ ) {
		if  ( pipe (child_rpipe[i]) == -1 || pipe(child_wpipe[i]) == -1 ) {
			printf ("Error creating pipe: %d\n", errno);
			exit(1);
		}

		child_pid[i] = fork();
		if ( child_pid[i] ) {  // in parent
			close(child_rpipe[i][0]);
			out[i] = fdopen (child_rpipe[i][1], "w");
			in[i] = fdopen (child_wpipe[i][0], "r");
			close(child_wpipe[i][1]);
		} else {  // in child
			in[i] = fdopen (child_rpipe[i][0], "r");
			close (child_rpipe[i][1]);
			close (child_wpipe[i][0]);
			out[i] = fdopen (child_wpipe[i][1], "w");
			flags |= i << (SMALLJAC_HIGH_SHIFT + 1);
			child = 1;
			break; // Only parent creates child processes
		}
	}

	if (child) {
		result = smalljac_Lpolys (curve, start, end, flags, smalljac_parallel_callback, out[i]);
		p = -1;
		fwrite(&p, sizeof(p), 1, out[i]);
		fclose (out[i]);
		fclose (in[i]);
		if ( result < 0 ) {
			printf ("smalljac_Lpolys returned error %ld\n", result);
			exit(-result);
		}
		exit(0);
	} else {
		// Fill out heads
		for ( i = 0 ; i < threads ; i++ ) {
			result = fread(&(stream_heads[i]), sizeof(stream_heads[i]), 1, in[i]);
			if (!result) { printf("Unexpected EOF (init, %d)\n", i); goto die; }
		}

		for (;;) {
			// Determine least prime; could be faster, but not worth it
			i = smalljac_parallel_argmin(stream_heads, threads);

			// All heads are empty (value -1)
			if (i == -1) { p = end; goto early_exit; }

			// Read in data for prime
			p = stream_heads[i];
			result = fread(&good, sizeof(good), 1, in[i]);
			if (!result) { printf("Unexpected EOF ($good, %d)\n", i); goto die; }
			result = fread(&n, sizeof(n), 1, in[i]);
			if (!result) { printf("Unexpected EOF ($n, %d)\n", i); goto die; }
			if ( n < 0 || n > 2*SMALLJAC_MAX_GENUS ) { printf("Read unexpected value n=%d from stream %d\n", n, i); goto die; }

			if (n) {			
				result = fread(&(a[0]), sizeof(a[0]), n, in[i]);
				if (result < n) { printf("Unexpected EOF ($a[%ld], %d)\n", result, i); goto die; }
			}

			// Send to user
			if ( p != last_filter ) { // Don't call callback for filtered p
				result = callback(curve, p, good, a, n, arg);
			}

			if ( !result ) {
				if (good >= 0) {
					goto early_exit;
				} else { // Filtration
					last_filter = p;
				}
			} else {
				last_filter = 0;
			}

			// Update stream head
			result = fread(&(stream_heads[i]), sizeof(stream_heads[i]), 1, in[i]);
			if (!result) { printf("Unexpected EOF ($p, %d)\n", i); goto die; }
		}

		return 0;
	}

	printf("Should never happen, line %d, file %s\n", __LINE__, __FILE__);
	exit(1);

 early_exit:
 	result = (long)p;
	for ( i = 0 ; i < threads ; i++ ) {
		fclose(in[i]);
		fclose(out[i]);
		waitpid(child_pid[i],&status,0);
		if (!WIFEXITED(status)) {
			printf("Unexpected result from waitpid()\n");
			exit(1);
		}
		if (WEXITSTATUS(status))
			result = -(long)WEXITSTATUS(status);
	}
	return result;

 die:
#if 0
	for ( i = 0 ; i < threads ; i++ ) {
		fclose(in[i]);
		fclose(out[i]);
	}
	return SMALLJAC_NODATA;
#else
	exit(1);
#endif
}
