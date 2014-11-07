# Main makefile for ADDA package
# Requires GNU make (at least 3.81) to execute. 
# Actual compiling goes in folders 'seq', 'mpi', and 'ocl' for sequential, parallel (MPI), and OpenCL (GPU accelerated)
# version respectively
# $Date::                            $
#
# Copyright (C) 2006-2014 ADDA contributors
# This file is part of ADDA.
#
# ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with ADDA. If not, see
# <http://www.gnu.org/licenses/>.

#=======================================================================================================================
# !!! This makefile accepts environmental variables CFLAGS, FFLAGS, CPPFLAGS, LDFLAGS, and LDLIBS. First three are
# passed to C, Fortran, and C++ compiler respectively, and last two are passed to the linker. Values specified below are
# added to the values given by the environment. This can be used in invocation of make as: 
# > export CFLAGS=...; make ...
# These variables can also be given in the command line to make, but it will override all values specified in this
# Makefile. This is especially critical for CFLAGS which should be overridden only with much care, including all
# relevant flags in the new value. Moreover, all other variables can be overridden by command line option to make. In
# particular, this can be relevant for specifying a different version of the compiler (see below).
#=======================================================================================================================

# Test version and features of GNU make
ifeq ($(.FEATURES),)
  $(error GNU make version 3.81 or higher is required)
endif
ifeq ($(filter else-if,$(.FEATURES)),)
  $(error make does not support else-if)
endif

#=======================================================================================================================
# Fixed definitions for non-compilation targets
#=======================================================================================================================

SHELL      := /bin/sh
# Separate directories to store Makefiles and .o and .d files for sequential and MPI executables
SEQ        := seq
MPI        := mpi
OCL        := ocl
# Names of executables
PROGSEQ    := adda
PROGMPI    := adda_mpi
PROGOCL    := adda_ocl
# Targets implying compilation (empty one is always included)
NONTRIVIAL := all seq mpi ocl
# Files to store last-used command lines for linker, C and Fortran compilers
LDOPTSFILE  := .ldopts
COPTSFILE   := .copts
CPPOPTSFILE := .cppopts
FOPTSFILE   := .fopts
# this is to simplify cleaning commands
OPTSFILES   := $(LDOPTSFILE) $(COPTSFILE) $(FOPTSFILE) $(CPPOPTSFILE)
# Common makefile to be included from the specific makefiles for seq and mpi version
COMMONMK    := ../common.mk

#=======================================================================================================================
# !!! Start of control section. Flags and options here can be modified by user. However, the default set of options may
# work out-of-box on some systems.

# First section contains general variables used by non-compilation targets (clean, install, etc.)
#=======================================================================================================================

# nothing here yet

#=======================================================================================================================
# Below are flags and options specific for compilation. If compiling with MPI or OpenCL support please also look at
# control sections in files 'mpi/Makefile' and 'ocl/Makefile' respectively.

# Everything from here down to the main action part of this makefile is skipped if targets specified in the command line
# to this makefile do not imply compilation.
ifneq ($(if $(MAKECMDGOALS),$(if $(filter $(NONTRIVIAL),$(MAKECMDGOALS)),1,),1),)
#=======================================================================================================================

# ---Optional flags---
# Uncomment those that you find appropriate; additional information is in specified headers. One can also export
# environmental variable OPTIONS before executing make, or define its value in the command line of make. For instance,
# 'make OPTIONS=DEBUG ...' or 'make OPTIONS+=DEBUG ...'. If several options need to be given, they should be given as
# one argument in quotes with its parts separated by spaces, e.g. 'make OPTIONS="DEBUG FFT_TEMPERTON" ...'. OPTIONS that
# are uncommented below are appended to the list specified elsewhere. Full list of possible options is the following:
VALID_OPTS := DEBUG DEBUGFULL FFT_TEMPERTON PRECISE_TIMING NOT_USE_LOCK ONLY_LOCKFILE NO_FORTRAN NO_CPP \
              OVERRIDE_STDC_TEST OCL_READ_SOURCE_RUNTIME CLFFT_APPLE SPARSE USE_SSE3 OCL_BLAS NO_SVNREV \
              ACCIMEXP
# Debug mode. By default, release configuration is used (no debug, no warnings, maximum optimization). DEBUG turns on
# producing debugging symbols (-g) and warnings and brings optimization down to O2 (this is required to produce all
# possible warnings by the compiler). DEBUGFULL turns off optimization completely (for more accurate debugging symbols)
# and turns on additional diagnostic messages in the source code (debug.c/h).
#override OPTIONS += DEBUG
#override OPTIONS += DEBUGFULL
#override OPTIONS += SPARSE ACCIMEXP
#USE_SSE3

# Temperton FFT (fft.h).
#override OPTIONS += FFT_TEMPERTON

# Precise timing (prec_timing.h).
#override OPTIONS += PRECISE_TIMING

# Controls the mode of file locking, if any (io.h). Use at maximum one of the following options.
#override OPTIONS += NOT_USE_LOCK
#override OPTIONS += ONLY_LOCKFILE

# This should be uncommented, if compiling of Fortran sources cause problems. However, currently a number of ADDA
# (optional) features rely on Fortran sources, in particular: IGT, Temperton FFT. This features will not work if line
# below is uncommented.
#override OPTIONS += NO_FORTRAN

# This should be uncommented, if compiling of C++ sources cause problems. However, currently a number of ADDA (optional)
# features rely on C++ sources, in particular: OpenCL version. This features will not work if line below is uncommented.
#override OPTIONS += NO_CPP

# ADDA code relies on certain parts of C99 standard, therefore support of this standard is checked at compile time (and
# error is produced if the test fails). However, if you believe that your compiler supports all the required features
# (listed in const.h), but do not define itself conforming to C99, you may uncomment the following option to override
# the test. Do it at your own risk!
#override OPTIONS += OVERRIDE_STDC_TEST

# By default OpenCL version of ADDA reads CL kernel source from file(s) *.cl during compilation and stores them as
# static string, which is compiled into the GPU code at runtime. Alternative option is to read these files at runtime.
# This creates hassle of keeping .cl files together with the executable, but may be a bit more robust/portable than
# the default approach. Uncomment the following option to enable reading of kernel sources at runtime.
#override OPTIONS += OCL_READ_SOURCE_RUNTIME

# APPLE clFFT (fft.h)
#override OPTIONS += CLFFT_APPLE

# used AMDs CLBlas library for iterative solver (currently just suported in bicg method, oo effect on the other methods)
#override OPTIONS += OCL_BLAS

# By default, this Makefile invokes a Bash script that will determine the current subversion revision and pass it to
# ADDA source (so it is produced in the output of 'adda -V'. If this causes any problems (including significant delays),
# the functionality can be disabled by uncommenting the following line.
#override OPTIONS += NO_SVNREV

# ---Compilers---
# Choose one of the following. Can also be specified from command line to make (see explanation above for OPTIONS),
# overriding definition below. To specify a different version of the compiler, e.g. 'gcc-4.7' instead of 'gcc', use
# command 'make CC=gcc-4.7' (analogously for Fortran and C++ compilers).
# gnu - tested for gcc 3.2.3 - 4.5.2 (most recent experience only with 4.3 and newer)
# intel - tested on icc 9.0 - 11.0
# compaq - tested on Compaq C V6.5-303 (dtk) - last tested in 2007
# ibm - tested on xlc 8.0 - last tested in 2008
# hpux - tested on ia64
# other -
# WARNING: Currently intel compiler version 11.0 has problems with CUDA version 3.2 (Nvidia) CL headers. In particular,
# cl_platform.h defines some types with __attribute__(vector...) not checking for support of these attributes (against
# gcc version). Compiling with option '-no-gcc' does not seem to work either. So the only feasible option to use icc
# with ocl seems to manually modify this (system) header or replace them with that from AMD SDK. Hopefully, this will
# be fixed in later versions of CUDA and/or icc.
COMPILER := gnu

# Additional options for compiler. Flags specified below are appended to the ones specified in the environment or
# command line of make (see explanation above for OPTIONS). The same value is used for all programming languages to be
# compiled and for linker. In particular, "-m32" may be used to force 32 bit compilation in 64 bit environment. However,
# it also requires proper libraries (especially, external ones, like FFTW3 or MPI) to be supplied.
override EXTRA_FLAGS +=

# --FFTW3 paths--
# Specify path to headers and libraries of FFTW3. Some systems do not need them at all, some specify special global
# variables (first or second 2 lines), on some - FFTW3 is installed under user account (next 2 lines). Under Windows it
# may be required to specify paths manually (last 2 lines). Relative (to location of this Makefile) paths should be
# immediately transformed into absolute ones using "$(abspath ...)". This two variables (FFTW3_INC_PATH and
# FFTW3_LIB_PATH) can also be defined in the enviroment or in the command line of make (see explanation above for
# OPTIONS).
#FFTW3_INC_PATH := $(FFTW_INC)
#FFTW3_LIB_PATH := $(FFTW_LIB)
#FFTW3_INC_PATH := $(FFTWINCLUDE)
#FFTW3_LIB_PATH := $(FFTWLIB)
#FFTW3_INC_PATH := $(HOME)/include
#FFTW3_LIB_PATH := $(HOME)/lib
#FFTW3_INC_PATH := "$(abspath ./../lib)"
#FFTW3_LIB_PATH := "$(abspath ./../lib)"

#=======================================================================================================================
# !!! End of control section. Everything below is not designed to be modified by user. However, advanced users may wish
# to modify some compilers flags below, especially when using 'other' compiler.

# Unconditional variables
#=======================================================================================================================

# C files are located in source folder (src/), other files may be added below
CSOURCE := ADDAmain.c CalculateE.c calculator.c chebyshev.c comm.c crosssec.c GenerateB.c interaction.c io.c \
           iterative.c linalg.c make_particle.c memory.c  mt19937ar.c param.c Romberg.c sinint.c somnec.c \
           timing.c vars.c
# Fortran files are located in src/fort folder, other files may be added below
FSOURCE := d07hre.f d09hre.f d113re.f d132re.f dadhre.f dchhre.f dcuhre.f dfshre.f dinhre.f drlhre.f dtrhre.f \
           propaesplibreintadda.f
F90SOURCE := 
# C++ files are located in src/cpp folder, other files may be added below
CPPSOURCE :=

# All object files and dependencies are defined in common.mk after the source lists are finalized

# Path to search for source files from the child Makefiles (in child folders) Those are used for vpath directives in
# child Makefiles
PARENT    := ..
FFOLDER   := fort
CPPFOLDER := cpp

LDLIBS    += -lm
DEPFLAG   := -MD
# Fortran and C++ sources generate a lot of warnings, we do not plan to investigate them
FWARN     := -w

#=======================================================================================================================
# Conditional variables that depend on the values of optional flags.
#=======================================================================================================================
#
# Process OPTIONS
# Produce warnings if any of OPTIONS are not present in VALID_OPTS list
$(foreach var,$(filter-out $(VALID_OPTS),$(OPTIONS)),\
              $(warning The option '$(var)' is not recognized. Please check the spelling.))

$(info --- Compilation options: ---)
ifneq ($(filter DEBUGFULL,$(OPTIONS)),)
  $(info Full debug mode)
  DBGLVL := 2
  CDEFS  += -DDEBUGFULL
  CSOURCE += debug.c
  ifneq ($(filter DEBUG,$(OPTIONS)),)
    $(warning DEBUG has no effect when DEBUGFULL is enabled)
  endif
else ifneq ($(filter DEBUG,$(OPTIONS)),)
  $(info Debug mode)
  DBGLVL := 1
  CDEFS  += -DDEBUG
else
  $(info Release mode)
  DBGLVL := 0
endif
ifneq ($(filter SPARSE,$(OPTIONS)),)
  $(info Sparse (non-FFT) mode)
  ifneq ($(filter FFT_TEMPERTON,$(OPTIONS)),)
    $(error SPARSE turns off all FFT-related code, so it is incompatible with FFT_TEMPERTON)
  endif
  ifneq ($(filter CLFFT_APPLE,$(OPTIONS)),)
    $(error SPARSE turns off all FFT-related code, so it is incompatible with CLFFT_APPLE)
  endif
  ifneq ($(filter PRECISE_TIMING,$(OPTIONS)),)
    $(error SPARSE is currently incompatible with PRECISE_TIMING)
  endif
  
  CDEFS += -DSPARSE
else
  CSOURCE += fft.c
  ifneq ($(filter FFT_TEMPERTON,$(OPTIONS)),)
    $(info Temperton FFT)
    CDEFS += -DFFT_TEMPERTON
    ifeq ($(filter NO_FORTRAN,$(OPTIONS)),)
      FSOURCE += cfft99D.f
    else
      $(error Temperton FFT (FFT_TEMPERTON) is implemented in Fortran, hence incompatible with NO_FORTRAN)
    endif
  else
    $(info FFTW3)
    LDLIBS += -lfftw3
    ifdef FFTW3_INC_PATH
      CFLAGS += -I$(FFTW3_INC_PATH)
    endif
    ifdef FFTW3_LIB_PATH
      LDFLAGS += -L$(FFTW3_LIB_PATH)
    endif
  endif
  ifneq ($(filter CLFFT_APPLE,$(OPTIONS)),)
    # Here only the info is printed, the main logic is in ocl/Makefile
    $(info Apple clFFT routines)
  endif
endif
ifneq ($(filter PRECISE_TIMING,$(OPTIONS)),)
  $(info Precise timing)
  CDEFS += -DPRECISE_TIMING
  CSOURCE += prec_time.c
endif
ifneq ($(filter NOT_USE_LOCK,$(OPTIONS)),)
  $(info No locks at all)
  CDEFS += -DNOT_USE_LOCK
  ifneq ($(filter ONLY_LOCKFILE,$(OPTIONS)),)
    $(warning ONLY_LOCKFILE has no effect when NOT_USE_LOCK is enabled)
  endif
else ifneq ($(filter ONLY_LOCKFILE,$(OPTIONS)),)
  $(info Only lock files (without system locks))
  CDEFS += -DONLY_LOCKFILE
endif
ifneq ($(filter NO_FORTRAN,$(OPTIONS)),)
  $(info Without Fortran sources)
  CDEFS += -DNO_FORTRAN
  FSOURCE =
  F90SOURCE =
endif
ifneq ($(filter NO_CPP,$(OPTIONS)),)
  $(info Without C++ sources)
  CDEFS += -DNO_CPP
  CPPSOURCE =
endif
ifneq ($(filter OVERRIDE_STDC_TEST,$(OPTIONS)),)
  $(info Overriding test for C99 conformance)
  CDEFS += -DOVERRIDE_STDC_TEST
endif
ifneq ($(filter USE_SSE3,$(OPTIONS)),)
  CDEFS += -DUSE_SSE3
  CFLAGS += -msse3
  $(info Using SSE3 optimizations)
endif
ifneq ($(filter OCL_READ_SOURCE_RUNTIME,$(OPTIONS)),)
  # Here only the info is printed, the main logic is in ocl/Makefile
  $(info Read CL sources at runtime)
endif
ifneq ($(filter OCL_BLAS,$(OPTIONS)),)
  # Here only the info is printed, the main logic is in ocl/Makefile
  $(info Using clBLAS)
endif
ifneq ($(filter NO_SVNREV,$(OPTIONS)),)
  CDEFS += -DNO_SVNREV
  $(info Skipping revision number)
else
  # Get revision number if possible
  REV=$(shell sh updsvnrev.sh)
  ifneq ($(REV),)
    $(info Revision $(REV))
  endif
endif
ifneq ($(filter ACCIMEXP,$(OPTIONS)),)
  CDEFS += -DACCIMEXP
  $(info Using accelerated imExp with precomputed tables)
endif
# Process EXTRA_FLAGS
ifneq ($(strip $(EXTRA_FLAGS)),)
  $(info Extra compiler options: '$(EXTRA_FLAGS)')
  CFLAGS   += $(EXTRA_FLAGS)
  FFLAGS   += $(EXTRA_FLAGS)
  CPPFLAGS += $(EXTRA_FLAGS)
  LDFLAGS  += $(EXTRA_FLAGS)
endif

# By default Fortran anc C++ optimization options are the same as C, but specific options can be used for some compilers
# below. This definition creates a permanent link, so further changes to COPT affect FOPT and CPPOPT as well.
FOPT   = $(COPT)
CPPOPT = $(COPT)
# This is for additional libraries that may be needed when using C linker on Fortran or C++ sources (using Fortran or
# C++ linker may also cause some problems, i.e. for MPI mode). Particular values are assigned for each compiler below.
FLIBS   +=
CPPLIBS +=
# Debugging flag is assumed the same for all compilers
DBGFLAG := -g
# Compiler-specific options
ifeq ($(COMPILER),gnu)
  CC      := gcc
  CSTD    := -std=c99
  COPT1   := -O2
  COPT2   := -O3 -ffast-math -funroll-loops
  CWARN   := -Wall -Wextra -Wcast-qual -Wpointer-arith -Wwrite-strings -Wstrict-prototypes \
             -Wstrict-aliasing=1 -Wshadow -Wcast-align -Wnested-externs -Wcomment -Wno-overlength-strings
  # gcc versions prior to 4.7.2 are affected by bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=7263 , which causes
  # -pedantic flag to generates warnings on every occurence of I (complex i)
  GCC_GTEQ_472 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' \
    -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40702)
  ifeq "$(GCC_GTEQ_472)" "1"
    CWARN += -pedantic
  endif
  
  # On Windows (with MinGW) there is certain motivation to link Fortran libraries statically, since otherwise three
  # DLLs (total size about 1.5MB) need to be distributed together with executables. Static linking increases the size
  # of ADDA by about 0.4MB, but we distribute at least 5 executables (not counting misc/). Moreover, complete static
  # linking is inconvenient since FFTW3 static libraries are not easily available for Windows (need to be recompiled
  # manually). Selective static linking can be achived by surrounding '-l...' flags in '-Wl,-Bstatic' ... 
  # '-Wl,-Bdynamic', but we do not do it here.
  #
  # Use gfortran if available (GCC 4 and later), otherwise try g77
  ifeq ($(shell which gfortran > /dev/null 2>&1 && echo 0),0)
    CF    := gfortran
    FLIBS += -lgfortran
  else
    # This is not expected to work for f90 sources but we keep it for now
    CF    := g77
    FLIBS += -lg2c
  endif

  CCPP    := g++
  CPPLIBS := -lstdc++
  # for now we do not want to investigate C++ warnings (since these sources are planned to be replaced by more advanced
  # routines), so we consider the following combination thorough enough
  CPPWARN := -Wall -Wextra
else ifeq ($(COMPILER),intel)
  CC    := icc
  CSTD  := -std=c99 -vec-report0 # the last flag is used to always remove vectorization remarks
  COPT1 := -O2
  COPT2 := -O3
  CWARN := -Wall -Wcheck -diag-disable 981,1418,1419,1498,1572,2259

  CF    := ifort
  FWARN += -vec-report0
  FLIBS += -lifcore

  CCPP  := icpc
  CPPWARN := -Wall -Wcheck -diag-disable 279,981,1418,1419
  # it seems that icpc relies on gcc stdc++ library anyway, but icc not always adds it during linking
  CPPLIBS += -lstdc++
  # if IPO is used, corresponding flags should be added to linker options: LDFLAGS += ...
else ifeq ($(COMPILER),compaq)
  # This compiler was not tested since 2007. In particular, warning options may not fit exactly the C99 standard, to
  # which the code was transferred. Its support for 64 bit compilations is also undefined. No C++ compiler is defined.
  # If you happen to use this compiler, please report results to the authors.
  CC    := cc
  CF    := f77
  #CSTD  = -std=c99
  COPT1 := -O2
  COPT2 := -fast
  CWARN := -w0 -msg_disable nestedcomment,unknownpragma,unreachcode
else ifeq ($(COMPILER),ibm)
  # This compiler was not tested since 2008. In particular, it is not clear, whether and what FLIBS should be used. No
  # C++ compiler is defined. If you happen to use this compiler, please report results to the authors.
  #
  # -O5 implies "-arch=auto", which tunes compilation exclusively for the host machine. However, it will not work in
  # some configurations. Then use '-O3 -qarch=... -qtune=...' instead
  CC    := xlc
  CF    := xlf
  CSTD  := -qlanglvl=extc99
  COPT1 := -O2
  COPT2 := -O3 -qcache=auto
  DEPFLAG := -qmakedep=gcc
  CWARN   := -qsuppress=1506-224:1506-342:1500-036
else ifeq ($(COMPILER),hpux)
  # This compiler was not tested since 2010. In particular, no C++ compiler is defined.  If you  happen to use this
  # compiler, please report results to the authors.
  CC    := cc
  CF    := f90
  CSTD  := -AC99
  COPT1 := +O2 +DD64
  COPT2 := +O3 +DD64
  DEPFLAG :+= +Md
  CWARN :=
  CFLAGS += -DNOT_USE_LOCK
else ifeq ($(COMPILER),other)
# add here definitions corresponding to 'other' compiler, if you want to use it.
else
  $(error Unknown compiler set '$(COMPILER)')
endif
$(info Compiler set '$(COMPILER)')

# if 'release' turn off warnings
ifeq ($(DBGLVL),0)
# the following two options (-w) are assumed universal across all compilers
  CWARN   := -w
  CPPWARN := -w
  LDFLAGS += -w
  COPT    := $(COPT2)
  DBGFLAG :=
else ifeq ($(DBGLVL),1)
  COPT := $(COPT1)
else ifeq ($(DBGLVL),2)
  COPT   :=
  FOPT   :=
  CPPOPT :=
endif
# Finalize option flags; these aggregates are used as a whole furher on
# LDFLAGS is finalized in common.mk
CFLAGS   += $(COPT) $(CWARN) $(DBGFLAG) $(CSTD) $(CDEFS)
FFLAGS   += $(FOPT) $(FWARN) $(DBGFLAG)
CPPFLAGS += $(CPPOPT) $(CPPWARN) $(DBGFLAG)

$(info ----------------------------)

#=======================================================================================================================
#end of check for targets implying compilation
endif
# Main action part
#=======================================================================================================================

.EXPORT_ALL_VARIABLES:
.PHONY: seq mpi ocl all cleanfull clean cleanruns cleanseq cleanmpi cleanocl cleanrunsseq cleanrunsmpi cleanrunsocl

all: seq mpi ocl

seq:
	@echo "Compiling sequential version of ADDA"
	$(MAKE) -C $(SEQ)

mpi:
	@echo "Compiling MPI version of ADDA"
	$(MAKE) -C $(MPI)

ocl:
	@echo "Compiling OpenCL version of ADDA"
	$(MAKE) -C $(OCL)

cleanfull: clean cleanruns

clean: cleanseq cleanmpi cleanocl

cleanruns: cleanrunsseq cleanrunsmpi cleanrunsocl

# 'clean' commands are executed here to keep them simple. Child makefiles are used only for compilation and thus contain
# quite heavy processing.
cleanseq:
	@echo "Removing sequential compiled files"
	cd $(SEQ) && rm -f *.o *.d $(OPTSFILES) $(PROGSEQ) $(PROGSEQ).exe

cleanmpi:
	@echo "Removing MPI compiled files"
	cd $(MPI) && rm -f *.o *.d $(OPTSFILES) $(PROGMPI) $(PROGMPI).exe

cleanocl:
	@echo "Removing OpenCL compiled files"
	cd $(OCL) && rm -f *.o *.d $(OPTSFILES) $(PROGOCL) $(PROGOCL).exe *.clstr
	
cleanrunsseq: 
	@echo "Removing output of sequential version of ADDA"
	cd $(SEQ) && rm -f -r ExpCount run* test*

cleanrunsmpi: 
	@echo "Removing output of MPI version of ADDA"
	cd $(MPI) && rm -f -r ExpCount run* test*
	
cleanrunsocl: 
	@echo "Removing output of OpenCL version of ADDA"
	cd $(OCL) && rm -f -r ExpCount run* test*
