#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpzpolyutil.h"
#include "igusa.h"
#include "cstd.h"

/*
    Copyright (c) 2011-2012 Andrew V. Sutherland
    See LICENSE file for license details.
*/

// g(x) = f(x+1) for f quintic
// f(x+1) = f5*x^5 + (f4 + 5*f5)*x^4 + (f3 + 4*f4 + 10*f5)*x^3 + (f2 + 3*f3 + 6*f4 + 10*f5)*x^2 + (f1 + 2*f2 + 3*f3 + 4*f4 + 5*f5)*x + f0 + f1 + f2 + f3 + f4 + f5
// g = f is ok
static void _mpz_poly_shift_d5 (mpz_t g[], mpz_t f[])
{
	static mpz_t w1, w2, w3, w4, w5;
	static int init;
	
	// use static allocation --- we may get called repeatedly
	if ( ! init ) { mpz_init (w1); mpz_init (w2); mpz_init (w3); mpz_init (w4); mpz_init (w5); init = 1; }
	mpz_mul_ui(w1,f[5],5);  mpz_mul_2exp(w2,f[4],2);  mpz_mul_ui(w3,f[3],3); mpz_mul_ui(w4,f[4],6); mpz_add (w5,w1,w1);
	mpz_add(g[0],f[0],f[1]); mpz_add(g[0],g[0],f[2]); mpz_add(g[0],g[0],f[3]); mpz_add(g[0],g[0],f[4]); mpz_add(g[0],g[0],f[5]); 
	mpz_add(g[1],f[1],f[2]);  mpz_add(g[1],g[1],f[2]); mpz_add(g[1],g[1],w3); mpz_add(g[1],g[1],w2); mpz_add(g[1],g[1],w1);
	mpz_add(g[2],f[2],w3); mpz_add(g[2],g[2],w4); mpz_add(g[2],g[2],w5);
	mpz_add(g[3],f[3],w2); mpz_add(g[3],g[3],w5);
	mpz_add(g[4],f[4],w1);
	mpz_set(g[5],f[5]);	
}

#define _mpz_mul2_ui(w1,w2,w3,c)	{ mpz_mul(w1,w2,w3); mpz_mul_ui(w1,w1,c); }

// computes the absolute igusa invariants for a genus 2 curve y^2=f(x) with f degree 5 or 6
int mpq_poly_igusa_inv (mpq_t I[3], mpz_t f[], int degree)
{
	static mpz_t f12, f22, f32, f42, f52, f62, f06, f13, f14, f15, f16, f23, f24, f25, f26, f34, f35, f36, f45, f46, I2, I4, I6, D, w1, w2, w3, w4, w5;
	static mpz_t g[7];
	static int init, ginit;

	// use static allocation --- we may get called repeatedly and we want this code to be as fast as possible
	if ( ! init ) { mpz_init (f12); mpz_init (f22); mpz_init (f32); mpz_init (f42); mpz_init (f52); mpz_init (f62); mpz_init (f06); mpz_init (f13); mpz_init (f14); mpz_init (f15); 
		           mpz_init (f23); mpz_init (f24); mpz_init (f25); mpz_init (f26); mpz_init (f34); mpz_init (f35); mpz_init (f45); mpz_init (f46); mpz_init(I2); mpz_init(I4); mpz_init(I6); mpz_init(D);
			   mpz_init (w1); mpz_init (w2); mpz_init (w3); mpz_init(w4); mpz_init(w5); init = 1; }
	// for degree 5 curves construct an isomorphic curve of degree 6
	if ( degree == 5 ) {
		int i;
		if ( ! ginit ) { for ( i = 0 ; i <= 6 ; i++ ) mpz_init(g[i]); ginit = 1; }
		// make constant term nonzero
		for ( i = 0 ; i <= 5 ; i++ ) mpz_set (g[i],f[i]);
		while ( ! mpz_sgn(g[0]) ) _mpz_poly_shift_d5 (g,g);
		// set g(x) = x^6*g(1/x), i.e. reverse coefficients and multiply by x
		mpz_set(g[6],g[0]);  mpz_swap(g[5],g[1]);  mpz_swap(g[4],g[2]); mpz_set_ui(g[0],0);
		return mpq_poly_igusa_inv(I,g,6);
	}
	if ( degree != 6 ) { err_printf ("Invalid degree %d in mpq_poly_igusa_inv, must be 5 or 6\n", degree);  abort(); }
//gmp_printf("%Zd*x^6+%Zd*x^5+%Zd*x^4+%Zd*x^3+%Zd*x^2+%Zd*x+%Zd\n", f[6], f[5], f[4], f[3], f[2], f[1], f[0]);
	/*
	    f12:=f1^2; f22:=f2^2; f32:=f3^2; f42:=f4^2; f52:=f5^2; f62:=f6^2;
	    f06:=f0*f6; f13:=f1*f3; f14:=f1*f4;  f15:=f1*f5; f16:=f1*f6;
	    f23:=f2*f3; f24:=f2*f4;  f25:=f2*f5; f26:=f2*f6; 
	    f34:=f3*f4; f35:=f3*f5; f36:=f3*f6; f45:=f4*f5; f46:=f4*f6;
	*/
	mpz_mul(f12,f[1],f[1]); mpz_mul(f22,f[2],f[2]); mpz_mul(f32,f[3],f[3]); mpz_mul(f42,f[4],f[4]); mpz_mul(f52,f[5],f[5]); mpz_mul(f62,f[6],f[6]);
	mpz_mul(f06,f[0],f[6]); mpz_mul(f13,f[1],f[3]); mpz_mul(f14,f[1],f[4]); mpz_mul(f15,f[1],f[5]); mpz_mul(f16,f[1],f[6]);
	mpz_mul(f23,f[2],f[3]); mpz_mul(f24,f[2],f[4]); mpz_mul(f25,f[2],f[5]); mpz_mul(f26,f[2],f[6]);
	mpz_mul(f34,f[3],f[4]); mpz_mul(f35,f[3],f[5]); mpz_mul(f36,f[3],f[6]); mpz_mul(f45,f[4],f[5]); mpz_mul(f46,f[4],f[6]);

	// I2:=-240*f06 + 40*f15 - 16*f24 + 6*f32;
	mpz_mul_ui(I2,f15,40); mpz_mul_ui(w1,f06,240); mpz_sub(I2,I2,w1); mpz_mul_2exp(w1,f24,4); mpz_sub(I2,I2,w1); mpz_mul_ui(w1,f32,6); mpz_add(I2,I2,w1);
//gmp_printf("I2 = %Zd\n", I2);

	// I4:=f0*(f6*(1620*f06 - 540*f15 - 504*f24 + 324*f32) + f5*(300*f25 - 180*f34) +  48*f4*f42) + f1*(f6*(300*f14 - 180*f23) - 12*f3*f42) + f15*(4*f24 - 80*f15 + 36*f32) + f22*(48*f26 - 12*f35 + 4*f42);
	if ( mpz_sgn(f[0]) ) {
		mpz_mul_ui(w1,f06,1620); mpz_mul_ui(w2,f15,540); mpz_sub(w1,w1,w2); mpz_mul_ui(w2,f24,504); mpz_sub(w1,w1,w2); mpz_mul_ui(w2,f32,324); mpz_add(w1,w1,w2); mpz_mul (w1,w1,f[6]);
		mpz_mul_ui(w2,f25,300); mpz_mul_ui(w3,f34,180); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[5]);  mpz_add(w1,w1,w2); mpz_mul(w2,f[4],f42); mpz_mul_ui(w2,w2,48); mpz_add(w1,w1,w2); mpz_mul(w1,w1,f[0]);
	} else {
		mpz_set_ui(w1,0);
	}
	if ( mpz_sgn(f[1]) ) { mpz_mul_ui(w2,f14,300); mpz_mul_ui(w3,f23,180); mpz_sub(w2,w2,w3);  mpz_mul(w2,w2,f[6]); _mpz_mul2_ui(w3,f[3],f42,12); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[1]); mpz_add(w1,w1,w2); }
	if ( mpz_sgn(f15) ) { mpz_mul_2exp(w2,f24,2); mpz_mul_ui(w3,f15,80); mpz_sub(w2,w2,w3); mpz_mul_ui(w3,f32,36); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f15); mpz_add(w1,w1,w2); }
	if ( mpz_sgn(f[2]) ) { mpz_mul_ui(w2,f26,48); mpz_mul_ui(w3,f35,12); mpz_sub(w2,w2,w3); mpz_mul_2exp(w3,f42,2); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f22); mpz_add(w1,w1,w2); }
	mpz_set (I4,w1);
//gmp_printf("I4 = %Zd\n", I4);
    
	/*
	I6:= f06*(f06*(59940*f15 - 119880*f06 + 20664*f24  - 10044*f32) +f0*(f5*(3060*f34 - 18600*f25) - 96*f4*f42) - 18600*f14*f16 - 2240*f15^2 + 3060*f13*f26 + 3472*f14*f25)
	       + f0*(f0*(2250*f35*f52 - 900*f42*f52) + f1*(1600*f25*f52 + 1818*f35*f36 - 876*f34*f46 - 1860*f34*f52 + 616*f42*f45)
		      + f2*(f26*(424*f42- 96*f26 - 876*f35) - 640*f24*f52 -  468*f32*f46 + 330*f32*f52 + f42*(492*f35 - 160*f42)) + f3*(f32*(162*f36 - 198*f45) + 60*f34*f42))
	       + f12*(f1*(f6*(2250*f36 +  1600*f45) - 320*f5*f52) - 900*f26^2 - 1860*f25*f36- 640*f26*f42 + 64*f24*f52 + 330*f32*f46 + 176*f35^2 + 26*f35*f42 - 36*f42^2)
	       + f1*(f22*(f6*(616*f25 + 492*f34) + f5*(26*f35 + 28*f42)) + f23*(f4*(76*f42-238*f35) - 198*f32*f6) + f3*f32*(72*f35 - 24*f42))
	       + f22*(f2*(f3*(60*f36 + 76*f45) - 24*f4*f42) - f22*(160*f46 + 36*f52) +  f32*(8*f42 - 24*f35));
	*/
	if ( mpz_sgn(f[0]) ) {
		mpz_mul_ui(w1,f15,59940); mpz_mul_ui(w2,f06,119880); mpz_sub(w1,w1,w2); mpz_mul_ui(w2,f24,20664); mpz_add(w1,w1,w2); mpz_mul_ui(w2,f32,10044); mpz_sub(w1,w1,w2); mpz_mul(w1,w1,f06);
		mpz_mul_ui(w2,f34,3060); mpz_mul_ui(w3,f25,18600); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[5]); _mpz_mul2_ui(w3,f[4],f42,96); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[0]); mpz_add(w1,w1,w2);
		_mpz_mul2_ui(w2,f14,f16,18600); mpz_sub(w1,w1,w2); _mpz_mul2_ui(w2,f15,f15,2240); mpz_sub(w1,w1,w2); _mpz_mul2_ui(w2,f13,f26,3060); mpz_add(w1,w1,w2);
		_mpz_mul2_ui(w2,f14,f25,3472); mpz_add(w1,w1,w2); mpz_mul(w1,w1,f06);
//gmp_printf ("w1=%Zd\n", w1);
		_mpz_mul2_ui(w2,f35,f52,2250); _mpz_mul2_ui(w3,f42,f52,900); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[0]); _mpz_mul2_ui(w3,f25,f52,1600); _mpz_mul2_ui(w4,f35,f36,1818); mpz_add(w3,w3,w4);
		_mpz_mul2_ui(w4,f34,f46,876); mpz_sub(w3,w3,w4); _mpz_mul2_ui(w4,f34,f52,1860); mpz_sub(w3,w3,w4); _mpz_mul2_ui(w4,f42,f45,616); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f[1]); mpz_add(w2,w2,w3);
		mpz_mul_ui(w3,f42,424); mpz_mul_ui(w4,f26,96); mpz_sub(w3,w3,w4); mpz_mul_ui(w4,f35,876); mpz_sub(w3,w3,w4); mpz_mul(w3,w3,f26); _mpz_mul2_ui(w4,f24,f52,640); mpz_sub(w3,w3,w4);
		_mpz_mul2_ui(w4,f32,f46,468); mpz_sub(w3,w3,w4); _mpz_mul2_ui(w4,f32,f52,330); mpz_add(w3,w3,w4); mpz_mul_ui(w4,f35,492); mpz_mul_ui(w5,f42,160); mpz_sub(w4,w4,w5); mpz_mul(w4,w4,f42); mpz_add(w3,w3,w4);
		mpz_mul(w3,w3,f[2]); mpz_add(w2,w2,w3); mpz_mul_ui(w3,f36,162); mpz_mul_ui(w4,f45,198); mpz_sub(w3,w3,w4); mpz_mul(w3,w3,f32); _mpz_mul2_ui(w4,f34,f42,60); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f[3]);
		mpz_add(w2,w2,w3); mpz_mul(w2,w2,f[0]); mpz_add(w1,w1,w2);
//gmp_printf("w2=%Zd\n", w2);
	} else {
		mpz_set_ui(w1,0);
	}
	mpz_mul_ui(w2,f36,2250); mpz_mul_ui(w3,f45,1600); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f[6]); _mpz_mul2_ui(w3,f[5],f52,320); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[1]); _mpz_mul2_ui(w3,f26,f26,900); mpz_sub(w2,w2,w3);
	_mpz_mul2_ui(w3,f25,f36,1860); mpz_sub(w2,w2,w3); _mpz_mul2_ui(w3,f26,f42,640); mpz_sub(w2,w2,w3); mpz_mul(w3,f24,f52); mpz_mul_2exp(w3,w3,6); mpz_add(w2,w2,w3); _mpz_mul2_ui(w3,f32,f46,330); mpz_add(w2,w2,w3);
	_mpz_mul2_ui(w3,f35,f35,176); mpz_add(w2,w2,w3); _mpz_mul2_ui(w3,f35,f42,26); mpz_add(w2,w2,w3); _mpz_mul2_ui(w3,f42,f42,36); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f12); mpz_add(w1,w1,w2);
//gmp_printf("w3=%Zd\n", w2);
	mpz_mul_ui(w2,f25,616); mpz_mul_ui(w3,f34,492); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f[6]); mpz_mul_ui(w3,f35,26); mpz_mul_ui(w4,f42,28); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f[5]); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f22);
	mpz_mul_ui(w3,f42,76); mpz_mul_ui(w4,f35,238); mpz_sub(w3,w3,w4); mpz_mul(w3,w3,f[4]); _mpz_mul2_ui(w4,f32,f[6],198); mpz_sub(w3,w3,w4); mpz_mul(w3,w3,f23); mpz_add(w2,w2,w3);
	mpz_mul_ui(w3,f35,72); mpz_mul_ui(w4,f42,24); mpz_sub(w3,w3,w4); mpz_mul(w3,w3,f32); mpz_mul(w3,w3,f[3]); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f[1]); mpz_add(w1,w1,w2);
//gmp_printf("w4=%Zd\n", w2);
	mpz_mul_ui(w2,f36,60); mpz_mul_ui(w3,f45,76); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f[3]); _mpz_mul2_ui(w3,f[4],f42,24); mpz_sub(w2,w2,w3); mpz_mul(w2,w2,f[2]);
	mpz_mul_ui(w3,f46,160); mpz_mul_ui(w4,f52,36); mpz_add(w3,w3,w4); mpz_mul(w3,w3,f22); mpz_sub(w2,w2,w3);
	mpz_mul_2exp(w3,f42,3); mpz_mul_ui(w4,f35,24); mpz_sub(w3,w3,w4); mpz_mul(w3,w3,f32); mpz_add(w2,w2,w3); mpz_mul(w2,w2,f22); mpz_add(w1,w1,w2);
//gmp_printf("w5=%Zd\n", w2);
	mpz_set (I6,w1);
//gmp_printf ("I6 = %Zd\n", I6);
	mpz_fast_poly_discriminant(D,f,degree);
	if ( ! mpz_sgn(D) ) return 0;
//gmp_printf("I10 = %Zd\n", D);
	mpz_mul(w1,I2,I2); mpz_mul(w2,w1,w1); mpz_mul(w2,w2,I2); mpq_set_num(I[0],w2); mpq_set_den(I[0],D); mpq_canonicalize(I[0]);
	mpz_mul(w2,w1,I2); mpz_mul(w2,w2,I4); mpq_set_num(I[1],w2); mpq_set_den(I[1],D); mpq_canonicalize(I[1]);
	mpz_mul(w2,w1,I6); mpq_set_num(I[2],w2); mpq_set_den(I[2],D); mpq_canonicalize(I[2]);
//gmp_printf("Absolute invariants %Qd, %Qd, %Qd\n", I[0], I[1], I[2]);
	return 1;
}
