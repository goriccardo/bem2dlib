LIBSRCS = pressure.f90 geom.f90 fieldgrid.f90 wake.f90 velocity.f90 geomwing.f90 integrals.f90 bem2d.f90 \
          geomcircle.f90 boundaryc.f90 ematrices.f90 clnaca00xx.f90 vgmeth.f90 theo.f90

SRCS = $(LIBSRCS) circle2d.f90 wing2d.f90 wing2dlap.f90 etest.f90 naca00xx.f90

OBJS = pressure.o geom.o fieldgrid.o wake.o velocity.o integrals.o bem2d.o boundaryc.o ematrices.o \
       geomwing.o vgmeth.o theo.o

COBJS = $(OBJS) circle2d.o geomcircle.o
WOBJS = $(OBJS) wing2d.o
WLOBJS = $(OBJS) wing2dlap.o
NOBJS = $(OBJS) naca00xx.o
EOBJS = $(OBJS) etest.o

LIBS = -llapack

F90 = gfortran
F90FLAGS = -O2 -g -Wall
LIBFLAGS = -shared -fPIC
LDFLAGS =

all: wing2d circle2d wing2dlap etest naca00xx libbem2d.so

wing2d: $(WOBJS)
	$(F90) $(LDFLAGS) -o $@ $(WOBJS) $(LIBS)

wing2dlap: $(WLOBJS)
	$(F90) $(LDFLAGS) -o $@ $(WLOBJS) $(LIBS)

circle2d: $(COBJS)
	$(F90) $(LDFLAGS) -o $@ $(COBJS) $(LIBS)

naca00xx: $(NOBJS)
	$(F90) $(LDFLAGS) -o $@ $(NOBJS) $(LIBS)

etest: $(EOBJS)
	$(F90) $(LDFLAGS) -o $@ $(EOBJS) $(LIBS)

libbem2d.so:  $(OBJS)
	$(F90) $(LDFLAGS) $(LIBFLAGS) -o $@ $(OBJS) $(LIBS)

dist:
	tar -cvzf bem2dlib.tar.gz Makefile LICENSE AUTHORS README $(SRCS)

clean:
	rm -f $(COBJS) $(WOBJS) $(WLOBJS) $(EOBJS)

python:
	f2py -c $(LIBS) -m bempy2d $(LIBSRCS)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<
