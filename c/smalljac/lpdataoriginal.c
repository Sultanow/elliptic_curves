#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ff_poly.h"
#include "mpzutil.h"
#include "smalljac.h"
#include "cstd.h"

/*
    Copyright (c) 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// Program to demonstrate basic use of the smalljac interface -- dumps Lpoly coefficient data to a file


#define delta_msecs(s,t)		(1000UL*(t-s)/CLOCKS_PER_SEC)

// Context structure for callback function - this is user defined, modify to suit
struct callback_ctx {
	unsigned long count;
	unsigned long missing_count;
	long trace_sum;
	FILE *fp;
};

// global variable

int sum = 0;
/*
	This callback function simply outputs L_p(T) coefficients a_1,...a_g to ctx->fp and computes a1_sum.
*/
int dump_lpoly (smalljac_curve_t curve, unsigned long p, int good, long a[], int n, void *arg)
{
	static int cnt;
	struct callback_ctx *ctx;

	ctx = (struct callback_ctx*) arg;
	ctx->count++;

	
	if ( ! n ) {
		printf ("Lpoly not computed at %ld\n", p); ctx->missing_count++; fprintf (ctx->fp, "%ld,?\n", p); 
		return 1;
	}
	ctx->trace_sum -= a[0];
	switch (n) {
	case 1: fprintf(ctx->fp,"%ld,%ld\n", p, a[0]); break;
	case 2: fprintf(ctx->fp,"%ld,%ld,%ld\n", p, a[0], a[1]); break;
	case 3: fprintf(ctx->fp,"%ld,%ld,%ld,%ld\n", p, a[0], a[1], a[2]); break;
	default: printf ("\rUnexpected number of Lpoly coefficients %d\n", n);  return 0;
	}
	if ( ! ((++cnt)&0xFFFF) ) printf ("\r%lu\r", p);  fflush(stdout);
	return 1;	
}


int main (int argc, char *argv[])
{
	time_t start_time, end_time;
	smalljac_curve_t curve;
	char filename[256];
	struct callback_ctx context;
	unsigned long flags;
	long result;
	int i, err, jobs, jobid;
	long minp, maxp;
	char *s,*t;
	
	if ( argc < 4 ) {
		printf ("Usage:    lpdata file-prefix curve-spec interval/bound [flags jobs job-id]\n\n");
		printf ("Examples:\n");
		printf ("          lpdata mycurve \"x^5 + 3145926x + 271828\" [1000..2000] 1\n");
		printf ("          lpdata 31a \"[1,-1-a,a,0,0] / (a^2-a-1)\" 10e6\n");
		printf ("          lpdata foo \"[x^3 + (z^3+z-1)*x + 3*z^2-4] / (z^4+z^3+z^2+z+1)\" 2e20 4\n");
		puts ("");
		printf ("flags&1 => SMALLJAC_GOOD_ONLY, flags&2 => SMALLJAC_A1_ONLY, flags&4 => SMALLJAC_DEGREE1.\n");
		printf ("smalljac version %s\n", SMALLJAC_VERSION_STRING);
		return 0;
	}
	
	flags = 0;
	jobs = jobid = 0;
	if ( argc > 4 ) {
		i = atol (argv[4]);
		if ( i < 0 || i > 7 ) { printf ("Unknown/unsupported flags specified\n"); return 0; }
		if ( i&1 ) flags |= SMALLJAC_GOOD_ONLY;
		if ( i&2 ) flags |= SMALLJAC_A1_ONLY;
		if ( i&4 ) flags |= SMALLJAC_DEGREE1_ONLY;
	}
	if ( argc > 5 ) {
		jobs = atoi (argv[5]);
		for ( i = 1 ; i <9 ;i++ ) if ( jobs == (1<<i) ) break;
		if ( i == 9 ) { printf ("When specified, jobs must be 2^k with 1 <= k <= 8\n"); return 0; }
		if ( argc < 7 ) { puts ("Please specify job-id also\n"); return 0; }
		flags |= ((1UL<<i)-1) << (SMALLJAC_SPLIT_SHIFT+1);
		jobid = (atoi(argv[6]) % jobs);
		flags |= jobid << (SMALLJAC_HIGH_SHIFT+1);
	}
	
	if ( (t = strstr(argv[3],"..")) ) {
		*t = '\0';
		s = argv[3]; if ( *s == '[' ) s++;
		minp = atol_exp(s);
		maxp = atol_exp(t+2);
	} else {
		minp = 1; maxp = atol_exp(argv[3]);
	}
	curve = smalljac_curve_init (argv[2], &err);
	if ( ! curve ) { 
		switch (err) {
		case SMALLJAC_PARSE_ERROR: printf ("Unable to parse curve string: %s\n", argv[2]); break;
		case SMALLJAC_UNSUPPORTED_CURVE: puts ("Specified curve not supported - should be degree 3, 5, or 7\n");  break;
		case SMALLJAC_SINGULAR_CURVE: puts ("Specified curve is singular\n");  break;
		default: printf ("smalljac_curve_init returned error %d\n", err);
		}
		return 0;
	}
	if ( smalljac_curve_nf_degree(curve) > 2 ) {
		printf ("Only computing L-poly coefficients at degree-1 primes over number field of degree %d\n", smalljac_curve_nf_degree(curve));
		flags |= SMALLJAC_DEGREE1_ONLY;
	}

	memset (&context,0,sizeof(context));
	
	if ( jobs ) sprintf (filename, "%s_lpdata_%d_%d.txt", argv[1], jobs, jobid); else sprintf (filename, "%s_lpdata.txt", argv[1]);
	context.fp = fopen (filename,"w");
	if ( ! context.fp ) { printf ("Error creating file %s\n", filename); return 0; }
	
	// write header line
	if ( argv[2][0] != '[' ) fprintf (context.fp, "[%s]", argv[2]); else fprintf (context.fp,"%s",argv[2]);
	fprintf (context.fp, " %ld %ld", minp, maxp);
	if ( jobs ) fprintf(context.fp, " %d %d", jobs, jobid);
	fprintf (context.fp, "\n");
	fflush (context.fp);		// important to flush before calling parallel Lpolys!!
	
	context.trace_sum = 0;

	// this is where everything happens...
	start_time = time(0);
	result = smalljac_parallel_Lpolys (curve, minp, maxp, flags, dump_lpoly, (void*)&context);
//	result = smalljac_Lpolys (curve, minp, maxp, flags, dump_lpoly, (void*)&context);
	end_time = time(0);
	
	fclose (context.fp);
	smalljac_curve_clear (curve);
	
	if ( result < 0 ) {  printf ("smalljac_Lpolys returned error %ld\n", result);  return 0; }
	printf ("trace sum is %ld\n", context.trace_sum);
	if ( context.missing_count ) printf ("%ld Lpolys not computed due to bad reduction\n", context.missing_count);
	printf ("Processed %ld primes in %ld seconds (%.3f ms/prime)\n", context.count,
		   end_time-start_time, (context.count? ((1000.0*(end_time-start_time))/context.count) : 0.0));
	printf ("Output written to file %s\n", filename);
}
