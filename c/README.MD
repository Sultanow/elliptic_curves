## A program that calculates the rank of an elliptic curve
Using [Sutherland's library smalljac](https://math.mit.edu/~drew/) we implement a `C` program that calculates the rank of elliptic curves in a parallel manner.

For example we call `result = smalljac_parallel_Lpolys(curve, minp, maxp, flags, dum)` where `curve=[0,1,1,-347,48061]` and `minp=1` and `maxp=100` and (by default) `flags=0` (possible values are SMALLJAC_GOOD_ONLY, SMALLJAC_A1_ONLY, SMALLJAC_DEGREE1_ONLY).

## Prerequisites
In order to get the program running, the library [GMP](https://gmplib.org/) (The GNU Multiple Precision Arithmetic Library) must be installed. Additionally [Sutherland's library ff_poly](https://math.mit.edu/~drew/) is required.

## Installation under Linux
Suppose `/data/ellipticcurverank` is the directory, into which we want to install our program.

Then download/checkout the folders *ff_poly* and *smalljac* into */data/ellipticcurverank/* and create the *bin* folders:

```bash
mkdir /data/ellipticcurverank/bin /data/ellipticcurverank/bin/ff_poly /data/ellipticcurverank/bin/ff_poly/lib /data/ellipticcurverank/bin/ff_poly/include
mkdir /data/ellipticcurverank/bin /data/ellipticcurverank/bin/smalljac /data/ellipticcurverank/bin/smalljac/lib /data/ellipticcurverank/bin/smalljac/include
```

To compile the program go into the directory */data/ellipticcurverank/ff_poly* and run `make install all`. After that go into the directory */data/ellipticcurverank/smalljac* and run again `make install all`.

## Running the program
In order to run the program call

```bash
./lpdata testf "[1,1,0,-15413,728554]" 100000
```

which yields writes output to testf_lpdata.txt and yields *"Output written to file testf_lpdata.txt"*.
