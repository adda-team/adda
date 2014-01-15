# Makefile for calc_pip_arg. Uses the default gfortran, which should be present on any Unix. Other compilers and
# optimization flags may also be used, but this makes little sense, since the routine # is a trivial one and is executed
# only a few times.

CF      = gfortran
FFLAGS  = -O $(EXTRA_FLAGS)
PROG    = calc_pip_arg
FSOURCE = calc_pip_arg.f90

srcdir = .
vpath %.f90 $(srcdir)/
vpath Makefile $(srcdir)/

#=======================================================================================================================

.PHONY: all clean

all: $(PROG)

$(PROG): $(FSOURCE) Makefile
	$(CF) -o $@ $(FFLAGS) $<

clean:
	rm -f $(PROG) $(PROG).exe
