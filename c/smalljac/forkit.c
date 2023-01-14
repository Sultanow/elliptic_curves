#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include "forkit.h"

/*
	Simple parallel processing for (not necessarily thread safe) applictions that process an input text file line by line and generate an output text file with one output line per input line.
	Ensures that line ordering in output file matches that of input file.
	
	This could all be done with pipes, but for the sake of simplicity (and debugging/monitoring) we use the filesystem.  We assume i/o is not the limiting factor.
*/

int forkit (char *infile, char *outfile, int threads, int bufsize, int (*child) (char *outbuf, char *inbuf, int id, void *ctx), void *ctx)
{
	FILE *in, *out, **infiles;
	char *inbuf, *outbuf;
	char *filename;
	int id, sts;
	char *s;
	long n;
	
	// do a quick sanity check on the input file before forking any children
	in = fopen (infile, "r");
	if ( ! in ) { fprintf (stderr, "Error opening file %s\n", infile); return -1; }
	fclose (in);

	if ( ! bufsize ) bufsize = FORKIT_DEFAULT_BUFSIZE;
	inbuf = malloc (bufsize);
	filename = malloc (strlen(outfile) + 32);

	if ( ! threads ) threads = sysconf(_SC_NPROCESSORS_ONLN);

	fflush(0);	// clear all buffers
	for ( id = 0 ; id < threads ; id++ ) if ( ! fork() ) break;
	
	// parent waits here for all childen to finish and then interleaves the results
	if ( id == threads ) {
		infiles = calloc (threads, sizeof(*infiles));

		n = 0;
		while ( wait (&sts) > 0 ) if ( sts ) n = -1;
		if ( n < 0 ) goto cleanup;

		for ( id = 0 ; id < threads ; id++ ) {
			sprintf (filename, "%s_%d_", outfile, id);
			infiles[id] = fopen (filename, "r");
			if ( ! infiles[id] ) { fprintf (stderr, "Error opening file %s\n", filename); n = -1; goto cleanup; }
		}
		out = fopen (outfile, "w");
		if ( ! out ) { fprintf (stderr, "Error creating output file %s\n", outfile); for ( id-- ; id >= 0 ; id-- ) fclose (infiles[id]);  free (filename);  n = -1; goto cleanup; }
		in = fopen (infile, "r");
		if ( ! in ) { fprintf (stderr, "Error opening file %s\n", infile); n = -1; goto cleanup; }
		for ( n = 0 ; fgets (inbuf, bufsize, in) ; n++ ) {
			if ( ! fgets (inbuf, bufsize, infiles[n%threads]) ) { fprintf (stderr, "Missing expected record %ld in forkit child output file\n", n); n = -1; goto cleanup; }
			fputs (inbuf, out);
		}
		fclose (out);
cleanup:
		for ( id = 0 ; id < threads ; id++ ) {
			if ( infiles[id] ) {
				if ( fgets (inbuf, bufsize, infiles[id]) ) { fprintf (stderr, "Unexpected record %s in forkit child %d output file\n", inbuf, id);  n = -1; goto cleanup; }
				fclose (infiles[id]);
			}
			sprintf (filename, "%s_%d_", outfile, id);
			remove (filename);
		}
		free (filename);
		free (infiles);
		free (inbuf);
		return n;
	}

	// only children get here
	outbuf = malloc (bufsize);
	in = fopen (infile, "r");
	if ( ! in ) { fprintf (stderr, "Error opening file %s\n", infile); exit(-1); }
	sprintf (filename, "%s_%d_", outfile, id);
	out = fopen (filename, "w");
	if ( ! out ) { fprintf (stderr, "Error creating file %s\n", filename); exit(-1); }
	free (filename);
	for ( n = 0 ; fgets (inbuf, bufsize, in) ; n++ ) {
		if ( (n%threads) != id ) continue;
		s = strchr(inbuf,'\n');
		if ( !s ) { fprintf (stderr, "forkit input buffer overflow at line %ld of file %s\n%s\n", n+1, infile, inbuf); exit(-1); }
		*s = '\0';
		if ( ! (*child) (outbuf, inbuf, id, ctx) ) { fprintf (stderr, "forkit child failed at line %ld of file %s\n%s\n", n+1, infile, inbuf); exit(-1); }
		// make sure we write exactly one line (we ignore anything after the first newline and add a newline if needed)
		for ( s = outbuf ; *s && *s != '\n' ; s++);
		if ( s-outbuf > bufsize-2 ) { fprintf (stderr, "forkit input buffer overflow at line %ld of file %s\n%s\n", n+1, infile, inbuf); exit(-1); }
		if ( ! *s ) *s ='\n';
		*++s = '\0';
		fputs (outbuf, out);
	}
	fclose (in); fclose (out);
	free (inbuf); free (outbuf);
	exit (0);
}
