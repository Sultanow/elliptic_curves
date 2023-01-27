# A program that calculates the rank of an elliptic curve
Using [Sutherland's library smalljac](https://math.mit.edu/~drew/) we implement a `C` program that calculates the rank of elliptic curves in a parallel manner.

For example we call `result = smalljac_parallel_Lpolys(curve, minp, maxp, flags, dum)` where `curve=[0,1,1,-347,48061]` and `minp=1` and `maxp=100` and (by default) `flags=0` (possible values are SMALLJAC_GOOD_ONLY, SMALLJAC_A1_ONLY, SMALLJAC_DEGREE1_ONLY).

## Prerequisites
In order to get the program running, the library [GMP](https://gmplib.org/) (The GNU Multiple Precision Arithmetic Library) must be installed. Additionally [Sutherland's library ff_poly](https://math.mit.edu/~drew/) is required.

In order to check the path of a required lib, for instance *gmp*, run the command `ld -lgmp --verbose`.

## Installation under Linux
Suppose `/data/ellipticcurverank` is the directory, into which we want to install our program.

Then download/checkout the folders *ff_poly* and *smalljac* into */data/ellipticcurverank/* and create the *bin* folders:

```bash
mkdir /data/ellipticcurverank/bin /data/ellipticcurverank/bin/ff_poly /data/ellipticcurverank/bin/ff_poly/lib /data/ellipticcurverank/bin/ff_poly/include
mkdir /data/ellipticcurverank/bin /data/ellipticcurverank/bin/smalljac /data/ellipticcurverank/bin/smalljac/lib /data/ellipticcurverank/bin/smalljac/include
```

The folder structure looks as follows, where *ff_poly* and *smalljac* are the original folders containing the sources and where the folders of the same name inside the bin directory contain the compilations that are created (and copied thereto) by the corresponding build process.

```bash
└── ellipticcurverank
    ├── bin
    │   ├── ff_poly
    │   │   ├── include
    │   │   │   └── ff_poly
    │   │   └── lib
    │   └── smalljac
    │       ├── include
    │       └── lib
    ├── ff_poly
    └── smalljac
```

To compile the program go into the directory */data/ellipticcurverank/ff_poly* and run `make install all`. After that go into the directory */data/ellipticcurverank/smalljac* and run again `make install all`. After the *ff_poly* build process the header file *ff_poly.h* is located in *bin/ff_poly/include* and all other header files are located in *bin/ff_poly/include/ff_poly* and the binaries are located in *bin/ff_poly/lib*:

```bash
└── ellipticcurverank
    └── bin
        └── ff_poly
            ├── include
            │   ├── ff_poly
            │   │   ├── asm.h
            │   │   ├── cstd.h
            │   │   ├── ff.h
            │   │   ├── ff2k.h
            │   │   ├── ffext.h
            │   │   ├── ffmontgomery64.h
            │   │   ├── ffpoly.h
            │   │   ├── ffpolybig.h
            │   │   ├── ffpolyfromroots.h
            │   │   ├── ffpolysmall.h
            │   │   ├── ntutil.h
            │   │   └── polyparse.h
            │   └── ff_poly.h
            └── lib
                └── libff_poly.a
```

The result of compiling *smalljac* are the files *smalljac.h* and *libsmalljac.a*:

```bash
└── ellipticcurverank
    └── bin
        └── smalljac
            ├── include
            │   └── smalljac.h
            └── lib
                └── libsmalljac.a
```

## Preparing the files
Suppose you have 10 files: prime_all_01990, prime_all_01991, prime_all_01992, prime_all_01993, prime_all_01994, prime_all_01995, prime_all_01996, prime_all_01997, prime_all_01998, prime_all_01999.

Concatenate them using the command `cat prime_all_* > prime_all_01990-01999` to one large file named prime_all_01990-01999. Then split this large file into 40 parts using the command `split -n l/40 -d -a 2 prime_all_01990-01999 subfile-`.

Run the `smalljac_runner.py` program as described in the section below using the command `python smalljac_runner.py > out.txt`.

In the next step concatenate all generated CSV files using the command `cat subfile-*.result.csv > prime_all_01990-01999.csv`.

Finally cleanup the subfiles using the command `rm -f subfile-*`.

## Running the program
In order to run the program call

```bash
./lpdata testf "[1,1,0,-15413,728554]" 100000
```

which writes an output file testf_lpdata.txt and yields *"Output written to file testf_lpdata.txt"*.

To let *lpdata* process a file line by line, for instance the file *candidates.txt* (see *data* folder in this repository), then call the Python script *smalljac_runner.py* (see *python* folder in this repository) as follows:

```bash
python smalljac_runner.py > out.txt
```
