#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "/zfs/users/ahatzi/ahatzi/ff_poly_v1.2.7/ff_poly.h"
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

float sum3 = 0;
float sum0 = 0;
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
	case 1: sum3 += ((a[0]+2)*log(p)/(p+1+a[0]));
			sum0 += a[0]*log(p)/p;
			break;
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

	FILE *file = fopen("candidates.txt", "r");
	    if (file == NULL) {
	        printf("Error opening file!\n");
	        exit(1);
	    }

	    //open the new file to write the modified contents
	    FILE *new_file = fopen("candidates_modified.txt", "w");
	        if (new_file == NULL) {
	            printf("Error creating new file!\n");
	            exit(1);
	        }
	    
	
	    int data[5];
	    char curvename[1000];
	    while (fscanf(file, "[ %*d, %d, %d, %d, %d, %d, %*d ]\n", &data[0], &data[1], &data[2], &data[3], &data[4]) == 5) {
	    
	     sprintf(curvename, "[%d, %d, %d, %d, %d]", data[0], data[1], data[2], data[3], data[4]);
	    flags = 0;
	    jobs = jobid = 0;
	    minp=1;
	    maxp=131072;
	    curve = smalljac_curve_init (curvename, &err);

		memset (&context,0,sizeof(context));
		
	 //   if ( jobs ) sprintf (filename, "%testfile_lpdata_%d_%d.txt", jobs, jobid); else sprintf (filename, "%testfile_lpdata.txt");
	  // 	context.fp = fopen (filename,"w");
	//    if ( ! context.fp ) { printf ("Error creating file %s\n", filename); return 0; }
	    	
	    	// write header line
	    
	    	fflush (context.fp);		// important to flush before calling parallel Lpolys!!
			fflush(new_file);	    	
	    	context.trace_sum = 0;
	    
	    	// this is where everything happens...
	    	result = smalljac_parallel_Lpolys (curve, minp, maxp, flags, dump_lpoly, (void*)&context);
	    	sum3 = sum3/log(maxp);
	    	sum0 = sum0/log(maxp) + 0.5;

	    	 // concatenate the contents of the line with a certain string
      	  char modified_line[1024];
       	 sprintf(modified_line,"%s,%f,%f\n", curvename, sum0,sum3);
        // write the concatenated string to the new file
      	  printf( "%s", modified_line);
	      


	   	
		
	    	// write the concatenated string to the new file
	  	 //	fprintf(new_file, "%s\n", modified_line);
	    	


	    //	fclose (context.fp);
	    	smalljac_curve_clear (curve);

	       
	    }
	
	    fclose(file);
		fclose(new_file);
	
}

