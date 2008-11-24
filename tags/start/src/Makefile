# Huh ? No compiling commands in Makefile ?
# Yep, that is right. At this moment the actual makefiles are
#	1) make_mpi_s (MPI in single precision)
#	2) make_mpi_d (MPI in double precision)
#	3) make_pvm_s (PVM in single precision)
#	4) make_seq_s (Sequential single precision)
#	5) make_seq_d (Sequential double precision)
# The motivation for this clumsily construct is one of a rather 
# inexperienced and hasty programmer. Merging these 5 + 1 
# different makefiles into is rather easy by giving each version
# a different name. 

reset:
	rm -f ExpCount
clean:
	rm -f *.o
	rm -f scatter
