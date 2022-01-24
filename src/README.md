Main source of ADDA, including Makefiles and other scripts. C++ and Fortran sources of auxiliary (optional) routines are placed in separate folders. The compilation itself happens in `seq`, `mpi`, or `ocl` folder depending on the mode - object and executable files are placed there.
* `cpp/` - C++ sources for Apple clFFT
* `fort/` - Fortran sources
  * `bessel.f90` - subroutines for calculating Bessel functions
  * `cfft99D.f` - source file for Temperton FFT
  * `propaesplibreintadda.f`, `d07hre.f`, `d09hre.f`, `d113re.f`, `d132re.f`, `dadhre.f`, `dchhre.f`, `dcuhre.f`, `dfshre.f`, `dinhre.f`, `drlhre.f`, `dtrhre.f` - subroutine to compute IGT and related numerical routines
* `mpi/Makefile` - makefile for MPI version (called from the main makefile)
* `ocl/Makefile` - makefile for OpenCL version (called from the main makefile)
* `seq/Makefile` - makefile for sequential version (called from the main makefile)
* `ADDAmain.c`, `CalculateE.c`, `calculator.c`, `chebyshev.c`, `cmplx.c/h`, `comm.c/h`, `const.h`, `crosssec.c/h`, `debug.c/h`, `fft.c/h`, `function.h`, `GenerateB.c`, `igt_so.c/h`, `interaction.c/h`, `io.c/h`, `iterative.c`, `linalg.c/h`, `make_particle.c`, `matvec.c`, `memory.c/h`, `oclcore.c/h`, `oclmatvec.c`, `os.h`, `param.c/h`, `parbas.h`, `prec_time.c/h`, `Romberg.c/h`, `sinint.c`, `sparse_ops.h`, `timing.c/h`, `types.h`, `vars.c/h` – C source and header files of ADDA (see [CodeDesign](https://github.com/adda-team/adda/wiki/CodeDesign))
* `Makefile` – main makefile
* `common.mk` - common part of child makefiles, including all compilation directives
* `iw_compile.bat` - batch script to compile ADDA with Intel compilers on Windows
* `mt19937ar.c/h` – source and header files for Mersenne twister random generator
* `oclkernels.cl` - OpenCL kernels
* `somnec.c` - routines to compute Sommerfeld integrals from NEC2C package
* `updgithash.sh` - script to obtain git hash of the source and store it into generated `githash.h`