August 29, 2012 smalljac version 4.0 README.txt

smalljac is a software library and a set of programs for computing zeta
functions and other invariants of low genus curves over finite fields.  More
specifically, given a curve C of genus 1 or 2 defined over Q or a quadratic
number field, smalljac can efficiently compute either:

   (a) the sequence of L-polynomials Lp(T) at each prime p of good reduction
        up to a specified norm bound.  Here Lp(T) is the numerator of the zeta
	function Z(C/Fp;T), equivalently, Lp(T)=T^{2g}chi(1/T), where chi(T)
	is the characteristic polynomial of the Frobenius endomorphism at p.

  (b) the sequence of isomorphism types of abelian groups corresponding
       to the Jacobian of the reduction of C at each prime of good reduction
       up to a specified norm bound.

For futher background on smalljac, see [1], which is the paper you should
cite if you use smalljac in your research.

New features in smalljac version 4:

1) Support for curves defined over quadratic number fields.

2) Support for genus 2 curves of the form y^2=f(x) where f(x) is an arbitrary
polynomial of degree 6.  Currently only L-polynomial computations are
supported for degree 6 curves; support for computing Jacobian group 
structures is planned for version 4.1.

3) Facilities that allow efficiently filtering for primes at which the reduction
of the Jacobian of the specified curve has prime order.  This is useful, for
example, when searching for amicable pairs of elliptic curves.

4) Sato-Tate group identification:  given a curve of genus 1 or 2, smalljac
can provisionally determine its Sato-Tate group out of the 3 possibilities in
genus 1 and the 52 possibilities in genus 2, according to the classification
given in [2].

5) Built-in support for parallel processing using up to 256 cores.

Version 4 also includes various bug fixes and many performance enhancements.
Support for genus 3 curves has been removed from version 4 (they were only
partially supported in version 3).  This was done to remove several external
dependencies.

IMPORTANT: the finite field arithmetic used by smalljac has now been moved
to a separate libary called "ff_poly".  You need to download the .tar files
for *both* ff_poly and smalljac, available at http://math.mit.edu/~drew

In order to build smalljac you need gcc installed on a 64-bit linux system
(life is too short to support 32-bit code).  You then need to install the
GMP library, which you can get at http://gmplib.org/.  Build ff_poly first,
and then build smalljac.  You will need to either install the header files
and library file for ff_poly in /usr/local/, or modify the smalljac makefile
to look elsewhere.

The interface to the smalljac library is specified in smalljac.h.  There are
also four programs included, that serve as examples of how to use
smalljac and are useful in their own right:

1) amicable: searches for amicable pairs and aliquot cycles related to an
elliptic curve over Q, as defined in [3].

2) lpdata: dumps L-polynomial data for a specified curve to a file.

3) lpoly: simply computes the L-polynomial of a specified curve at a
specified prime.

4) moments: computes moments of L-polynomial coefficients of a
specified curve and attempts to provisionally identify its Sato-Tate group.

The command line interface to each of the programs above can
be obtained by running the program with no arguments.

The smalljac software is licensed under GPL version 2 or later.


References:

[1] Kiran S. Kedlaya and Andrew V. Sutherland,
"Computing L-series of hyperelliptic curves", ANTS VIII, LNCS 5011,
Springer, 2008, pp. 312-326.
(preprint: http://arxiv.org/abs/arXiv:0801.2778)

[2] Francesc Fite, Kiran S. Kedlaya, Victor Rotger, and Andrew V.
Sutherland, "Sato-Tate distributions and Galois endomorphism
modules in genus 2", to appear in Compositio Mathematica,
http://arxiv.org/abs/1110.6638.

[3] Katherine E. Stange and Joseph H. Silverman,
"Amicable pairs and aliquot cycles for elliptic curves",
Experimental Mathematics 20 (2011), pp. 329--357.
(preprint: http://arxiv.org/abs/0912.1831)
