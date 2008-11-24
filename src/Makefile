# Actual compiling goes in make_seq and make_mpi for
# sequential and parallel (mpi) version respectively

CSOURCE = CD2main.c  calculator.c make_particle.c CalculateE.c \
          GenerateB.c cmplx.c debug.c\
          CGNR.c MatVec_nim.c memory.c comm.c matvec.c \
          fft.c crosssec.c Romberg.c read_cluster.c
FSOURCE = cfft99D.f

export COBJECTS= $(CSOURCE:.c=.o)
#export CFLAGS= -DDEBUG
export FOBJECTS= $(FSOURCE:.f=.o)


seq:
	$(MAKE) -f make_seq
mpi:
	$(MAKE) -f make_mpi
reset:
	rm -f ExpCount
clean:
	rm -f *.o
	rm -f scatter
