SRCS = pressure.f90 geom.f90 fieldgrid.f90 wake.f90 velocity.f90 geomwing.f90 integrals.f90 bem2d.f90 \
       geomcircle.f90 circle2d.f90 crappyvel.f90

OBJS = pressure.o geom.o fieldgrid.o wake.o velocity.o integrals.o bem2d.o

COBJS = $(OBJS) circle2d.o geomcircle.o
WOBJS = $(OBJS) wing2d.o geomwing.o

LIBS = -llapack

F90 = gfortran
F90FLAGS = -O2 -Wall
LDFLAGS =

all: wing2d circle2d

wing2d: $(WOBJS)
	$(F90) $(LDFLAGS) -o $@ $(WOBJS) $(LIBS)
    
circle2s: $(COBJS)
	$(F90) $(LDFLAGS) -o $@ $(COBJS) $(LIBS)

clean:
	rm -f $(COBJS) $(WOBJS)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<
