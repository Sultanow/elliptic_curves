#ifndef _LPPLOT_INCLUDE_
#define _LPPLOT_INCLUDE_

/*
    Copyright (c) 2007-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

// YSCALE factors are relative to the height of the uniform distribution at 1.0
#define LPPLOT_GIF_G1_A1_YSCALE		1.5
#define LPPLOT_GIF_G1_S_YSCALE		1.5
#define LPPLOT_GIF_G2_A1_YSCALE		4.5
#define LPPLOT_GIF_G2_A2_YSCALE		5.5
#define LPPLOT_GIF_G2_S_YSCALE		4.5
#define LPPLOT_GIF_G3_A1_YSCALE		6.0
#define LPPLOT_GIF_G3_A2_YSCALE		15.0
#define LPPLOT_GIF_G3_A3_YSCALE		21.0
#define LPPLOT_GIF_G3_S_YSCALE		4.0

#define LPPLOT_GIF_Y				640						// fixed for all plots
#define LPPLOT_GIF_X				1000					// fixed for all plots
#define LPPLOT_GIF_XCROP			975						// crop value for making thumbnails
#define LPPLOT_GIF_YCROP1			50//40					// ditto
#define LPPLOT_GIF_YCROP2			620						// ditto
#define LPPLOT_GIF_TNAIL_DELAY		1.0						// delay between frames, in seconds
#define LPPLOT_GIF_DELAY			1.0						// ditto
#define LPPLOT_GIF_TNAIL_SCALE		0.05						// scaling factor for thumbnails

#define LPPLOT_HISTOGRAM_LINE_STYLE		3			// bluish
#define LPPLOT_UNIFORM_LINE_STYLE		0			// light grey
#define LPPLOT_DISTRIBUTION_LINE_STYLE	2			// light greenish

#define LPPLOT_MIN_LG_N	10
#define LPPLOT_MAX_LG_N	40
#define LPPLOT_MAX_LEVELS   (LPPLOT_MAX_LG_N-LPPLOT_MIN_LG_N+1)
#define LPPLOT_MOMENTS	11	// includes moment 0

// lpplot_nbuckets[i] = floor(sqrt(pi(2^i))) - OEIS sequence A133498
long lpplot_nbuckets[LPPLOT_MAX_LG_N+1] = { 0, 1, 1, 2, 2, 3, 4, 5, 7, 9, 13, 17, 23, 32, 43, 59, 80, 110, 151, 208, 286, 394, 544, 751, 1038, 1436, 1989, 2757,
									    3825, 5309, 7375, 10251, 14257, 19839, 27621, 38473, 53613, 74742, 104241, 145436, 202985 };

// lpplot_joint_nbuckets[i] = floor(cbrt(pi(2^i)))
long lpplot_joint_nbuckets[LPPLOT_MAX_LG_N+1] = { 0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 10, 12, 15, 18, 23, 28, 35, 43, 53, 66, 82, 102, 127, 158, 196, 244, 304, 378, 471, 587, 732, 913, 1139, 1421, 1774, 2214,  2765, 3453 };

double lpplot_a1_radius[5] = { 0.0, 2.0, 4.0, 6.0, 8.0 };	// binom (2g 1) - also the radius for s2, s3, ...
double lpplot_a2_radius[5] = { 0.0, 1.0, 6.0, 15.0, 28.0 };	// binom (2g 2)
double lpplot_a3_radius[5] = { 0.0, 0.0, 4.0, 30.0, 56.0 };	// binom (2g 3)

#ifdef __cplusplus
}
#endif

#endif
