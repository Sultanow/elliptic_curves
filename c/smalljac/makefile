CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -m64 -pedantic -std=gnu99
LDFLAGS = -static
INCLUDES = -I/usr/include -I/data/ellipticcurverank/bin/ff_poly/include
LIBS = -L/data/ellipticcurverank/bin/ff_poly/lib -lff_poly -lgmp -lm
INSTALL_ROOT = /data/ellipticcurverank/bin/smalljac

HEADERS = ecurve.h ecurve_ff2.h g2tor3poly.h hecurve.h hcpoly.h igusa.h jac.h jacorder.h lpplot.h nfpoly.h pointcount.h smalljac_g23.h smalljac_internal.h smalljactab.h bitmap.h cstd.h mpzpolyutil.h mpzutil.h ntutil.h polyparse.h prime.h
OBJECTS = ecurve.o ecurve_ladic.o ecurve_ff2.o hcpoly.o hecurve.o hecurve1.o hecurve2_ladic.o hecurve2.o igusa.o jac.o jacorder.o jacstructure.o nfpoly.o pointcount.o \
                  prime.o smalljac.o smalljac_moments.o smalljac_parallel.o smalljac_special.o smalljactab.o smalljac_g23.o smalljac_tiny.o STgroups.o  mpzpolyutil.o mpzutil.o polyparse.o
PROGRAMS = amicable lpdata lpoly moments

all: libsmalljac.a $(PROGRAMS)

clean:
	rm -f *.o
	rm -f libsmalljac.a $(PROGRAMS)

install: all
	cp -v smalljac.h $(INSTALL_ROOT)/include
	cp -v libsmalljac.a $(INSTALL_ROOT)/lib

##### smalljac library

libsmalljac.a: $(OBJECTS)
	ar -r libsmalljac.a $(OBJECTS)
	ranlib libsmalljac.a

##### executables

amicable: amicable.o libsmalljac.a smalljac.h
	$(CC) $(LDFLAGS) -o $@ $< libsmalljac.a $(LIBS)
	
lpoly: lpoly.o libsmalljac.a smalljac.h
	$(CC) $(LDFLAGS) -o $@ $< libsmalljac.a $(LIBS)
	
lpdata: lpdata.o libsmalljac.a smalljac.h
	$(CC) $(LDFLAGS) -o $@ $< libsmalljac.a $(LIBS)
	
moments: moments.o libsmalljac.a smalljac.h
	$(CC) $(LDFLAGS) -o $@ $< libsmalljac.a $(LIBS)

##### C modules
	
ecurve.o: ecurve.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ecurve_ladic.o: ecurve_ladic.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

ecurve_ff2.o : ecurve_ff2.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

amicable.o: amicable.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

lpoly.o : lpoly.c  smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

lpdata.o : lpdata.c  smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

lpplot.o : lpplot.c  smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

lprender.o : lprender.c  smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

moments.o : moments.c  smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hcpoly.o: hcpoly.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hecurve.o: hecurve.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hecurve1.o: hecurve1.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hecurve2.o: hecurve2.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

hecurve2_ladic.o: hecurve2_ladic.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

igusa.o: igusa.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

jac.o: jac.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

jacorder.o: jacorder.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

jacstructure.o: jacstructure.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	
nfpoly.o: nfpoly.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	
pointcount.o: pointcount.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

prime.o: prime.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac.o: smalljac.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljactab.o: smalljactab.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac_g23.o: smalljac_g23.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac_moments.o: smalljac_moments.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac_parallel.o: smalljac_parallel.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac_special.o: smalljac_special.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac_tiny.o: smalljac_tiny.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

STgroups.o: STgroups.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

mpzpolyutil.o: mpzpolyutil.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

mpzutil.o: mpzutil.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

polyparse.o: polyparse.c smalljac.h
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<

smalljac.h: $(HEADERS)
	touch smalljac.h
