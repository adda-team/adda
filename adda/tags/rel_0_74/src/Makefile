# Main makefile for ADDA package
# Actual compiling goes in make_seq and make_mpi for
# sequential and parallel (mpi) version respectively
# AUTH: Maxim Yurkin

CSOURCE  = ADDAmain.c  calculator.c make_particle.c CalculateE.c \
           GenerateB.c debug.c timing.c linalg.c \
           iterative.c memory.c comm.c matvec.c param.c \
           fft.c crosssec.c Romberg.c io.c prec_time.c vars.c
FSOURCE  = cfft99D.f  

COBJECTS = $(CSOURCE:.c=.o)
FOBJECTS = $(FSOURCE:.f=.o)
LDLIBS   = -lfftw3 -lm
PROG     = adda
SHELL    = /bin/sh
FFLAGS   = -w

# information on the following optional flags 
# can be found in corresponding headers

# for faster recompiling turn on in debug.h
#CFLAGS += -DDEBUG             

# for faster recompiling turn on in fft.h
#CFLAGS += -DFFT_TEMPERTON

# for faster recompiling turn on in timing.h
#CFLAGS += -DPRECISE_TIMING

# for faster recompiling turn on in io.h
#CFLAGS += -DNOT_USE_LOCK

export CSOURCE
export COBJECTS
export FOBJECTS
export CFLAGS
export FFLAGS
export LDLIBS
export PROG
export SHELL

seq:
	$(MAKE) -f make_seq

mpi:
	$(MAKE) -f make_mpi

reset:
	rm -f ExpCount

clean:
	rm -f *.o
	rm -f *.d
	rm -f $(PROG)
