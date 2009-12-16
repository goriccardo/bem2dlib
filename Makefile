LIBSRCS = pressure.f90 geom.f90 fieldgrid.f90 wake.f90 velocity.f90 geomwing.f90 integrals.f90 bem2d.f90 \
          geomcircle.f90 crappyvel.f90 boundaryc.f90

SRCS = $(LIBSRCS) circle2d.f90 wing2d.f90

OBJS = pressure.o geom.o fieldgrid.o wake.o velocity.o integrals.o bem2d.o boundaryc.f90

COBJS = $(OBJS) circle2d.o geomcircle.o
WOBJS = $(OBJS) wing2d.o geomwing.o
WLOBJS = $(OBJS) wing2dlap.o geomwing.o

LIBS = -llapack

F90 = gfortran
F90FLAGS = -O2 -Wall
LDFLAGS =

all: wing2d circle2d wing2dlap

wing2d: $(WOBJS)
	$(F90) $(LDFLAGS) -o $@ $(WOBJS) $(LIBS)

wing2dlap: $(WLOBJS)
	$(F90) $(LDFLAGS) -o $@ $(WLOBJS) $(LIBS)

circle2d: $(COBJS)
	$(F90) $(LDFLAGS) -o $@ $(COBJS) $(LIBS)

clean:
	rm -f $(COBJS) $(WOBJS)

python:
	f2py -c $(LIBS) -m bempy2d $(LIBSRCS)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<
