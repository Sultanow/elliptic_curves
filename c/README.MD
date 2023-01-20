# A program that calculates the rank of an elliptic curve
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

The folder structure looks as follows, where *ff_poly* and *smalljac* are the original folders containing the sources and where the folders of the same name inside the bin directory contain the compilations that are created (and copied thereto) by the corresponding build process.

```bash
├── ellipticcurverank
│   ├── bin
│   │   ├── ff_poly
│   │   │   ├── include
│   │   │   │   └── ff_poly
│   │   │   └── lib
│   │   └── smalljac
│   │       ├── include
│   │       └── lib
│   ├── ff_poly
│   └── smalljac
```

To compile the program go into the directory */data/ellipticcurverank/ff_poly* and run `make install all`. After that go into the directory */data/ellipticcurverank/smalljac* and run again `make install all`. After the *ff_poly* build process the header file *ff_poly.h* is located in *bin/ff_poly/include*:

```bash
├── ellipticcurverank
│   ├── bin
│   │   ├── ff_poly
│   │   │   ├── include
│   │   │   │   ├── ff_poly
│   │   │   │   │   ├── asm.h
│   │   │   │   │   ├── cstd.h
│   │   │   │   │   ├── ff.h
│   │   │   │   │   ├── ff2k.h
│   │   │   │   │   ├── ffext.h
│   │   │   │   │   ├── ffmontgomery64.h
│   │   │   │   │   ├── ffpoly.h
│   │   │   │   │   ├── ffpolybig.h
│   │   │   │   │   ├── ffpolyfromroots.h
│   │   │   │   │   ├── ffpolysmall.h
│   │   │   │   │   ├── ntutil.h
│   │   │   │   │   └── polyparse.h
│   │   │   │   └── ff_poly.h
│   │   │   └── lib
│   │   │       └── libff_poly.a
│   │   └── smalljac
│   │       ├── include
│   │       │   └── smalljac.h
│   │       └── lib
│   │           └── libsmalljac.a
```

## Running the program
In order to run the program call

```bash
./lpdata testf "[1,1,0,-15413,728554]" 100000
```

which writes an output file testf_lpdata.txt and yields *"Output written to file testf_lpdata.txt"*.
