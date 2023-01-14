CC = gcc
CFLAGS = -O3 -fomit-frame-pointer -funroll-loops -m64 -pedantic -std=gnu99
LDFLAGS = -static
INCLUDES = -I/usr/include
LIBS = -lgmp -lm 
INSTALL_ROOT = ~/bin/ff_poly

HEADERS = asm.h ff.h ffmontgomery64.h ff2k.h ffext.h ffpolyfromroots.h ffpolybig.h ffpolysmall.h ffpoly.h cstd.h ntutil.h polyparse.h
OBJECTS = ff.o ff2k.o  ffext.o  ffpolyfromroots.o ffpolysmall.o ffpoly_small.o polyparse.o

all: libff_poly.a

clean:
	rm -f *.o 
	rm -f libff_poly.a
	
install: all
	cp -v ff_poly.h $(INSTALL_ROOT)/include
	mkdir -p $(INSTALL_ROOT)/include/ff_poly
	cp -v $(HEADERS) $(INSTALL_ROOT)/include/ff_poly
	cp -v libff_poly.a $(INSTALL_ROOT)/lib
	
##### ff_poly library (standard version only supports small polys (degree < 256) and does not depend on zn_poly)

libff_poly.a: $(OBJECTS)
	ar -r libff_poly.a $(OBJECTS)
	ranlib libff_poly.a
	
##### executables

##### FF_POLY C modules

ff.o: ff.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o ff.o -c ff.c

ff2k.o : ff2k.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o ff2k.o -c ff2k.c

ffext.o: ffext.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o ffext.o -c ffext.c

ffpoly_small.o: ffpoly.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o ffpoly_small.o -c ffpoly.c

ffpolyfromroots.o: ffpolyfromroots.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o ffpolyfromroots.o -c ffpolyfromroots.c

ffpolysmall.o: ffpolysmall.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o ffpolysmall.o -c ffpolysmall.c

polyparse.o: polyparse.c ff_poly.h
	$(CC) $(CFLAGS) $(INCLUDES) -o polyparse.o -c polyparse.c

ff_poly.h: $(HEADERS)
	touch ff_poly.h
