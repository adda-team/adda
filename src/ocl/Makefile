# Makefile for ocl version of ADDA package
# All options are defined in Makefile
# $Date::                            $
#
# Copyright (C) 2010-2014 ADDA contributors
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
# !!! Start of control section. Flags and options here are designed to be modified by user to choose particular options
# for compilation. However, the default set of options may work out-of-box on some systems.
#=======================================================================================================================

# Depending on a particular OpenCL installation one may need to manually specify paths to corresponding headers and
# libraries. Path should be absolute or relative to the location of this Makefile.
# On Windows it is also possible to link to OpenCL.dll directly.
#CFLAGS += -I"C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v3.2/include"
#LDFLAGS += -L"C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v3.2/lib/Win32"

# This workaround is to disable linker warnings when linking against 32-bit OpenCL.dll on 64-bit Windows. The issue 
# appeared with Nvidia library supplied with driver 301.42. Uncomment the following 3 lines to enable it.
#ifneq ($(filter -m32,$(LDFLAGS)),)
#  LDFLAGS += -Wl,--enable-stdcall-fixup
#endif

ifneq ($(filter SPARSE,$(OPTIONS)),)
  $(error OPENCL version of sparse mode is not yet available)
else # the following is only for non-sparse mode
  ifneq ($(filter CLFFT_APPLE,$(OPTIONS)),)
    ifeq ($(filter NO_CPP,$(OPTIONS)),)
      CPPSOURCE += fft_execute.cpp fft_setup.cpp fft_kernelstring.cpp
      CFLAGS    += -DCLFFT_APPLE
    else
      $(error Apple clFFT (CLFFT_APPLE) is implemented in C++, hence incompatible with NO_CPP)
    endif
  else # AMD clFFT flags
    # Depending on a particular clFFT installation one may need to manually specify paths to its headers and 
    # libraries. Path should be absolute or relative to the location of this Makefile.
    # On Windows (as in example below) it is recommended to link to *.dll instead of *.lib
    #CFLAGS  += -I"C:\Program Files (x86)\AMD\clAmdFft\include"
    #LDFLAGS += -L"C:\Program Files (x86)\AMD\clAmdFft\bin64"
  
    LDLIBS += -lclFFT

    # On Unix clFFT may need pthread library, which is not linked by default; if yes, uncomment the following
    #LDLIBS += -lpthread
  
    # Uncomment the following line if you experience problems during linking, like
    # "...libclFFT.so: undefined reference to `clRetainContext@OPENCL_1.0'"
    # This is caused by Nvidia OpenCL libraries on Linux (hopefully it will disappear in next releases of this library)
    # The drawback of this workaround is that a warning would appear at each execution of adda_ocl. If you don't like 
    # it, provide another libOpenCL.so (e.g. that from AMD, or better - open-source one from 
    # https://forge.imag.fr/projects/ocl-icd/) for linking.
    #LDFLAGS +=-Wl,--unresolved-symbols=ignore-in-shared-libs
  endif
  ifneq ($(filter OCL_BLAS,$(OPTIONS)),)
    # Depending on a particular clBLAS installation one may need to manually specify paths to its headers and 
    # libraries. Path should be absolute or relative to the location of this Makefile.
    # On Windows (as in example below) it is recommended to link to *.dll instead of *.lib
    #CFLAGS  += -I"C:\Program Files (x86)\AMD\clAmdBlas\include"
    #LDFLAGS += -L"C:\Program Files (x86)\AMD\clAmdBlas\bin64"
    
    LDLIBS += -lclBLAS

    CFLAGS += -DOCL_BLAS
    # If you have any problems during linking, see above for issues related to clFFT. The same solutions may work
    # for clBLAS
  endif
endif

#=======================================================================================================================
# !!! End of control section. Everything below is not designed to be modified by user
#=======================================================================================================================

CFLAGS  += -DOPENCL 
LDLIBS  += -lOpenCL
CSOURCE += oclcore.c oclmatvec.c

ifneq ($(filter OCL_READ_SOURCE_RUNTIME,$(OPTIONS)),)
  CFLAGS += -DOCL_READ_SOURCE_RUNTIME
else
  # Each of the following .cl files is automatically transformed into .clstr files (located at ocl/) during compile
  # time. The latter contain all the content of the original file as a C string, and are included in the C source files
  # to define static strings. These strings define (macro) are further used as a source to compile OpenCL kernels.
  CLSOURCE := oclkernels.cl
endif

# This is a weird prefix, which is required since that is how .clstr files are quoted in .d files
CLSPREFIX  := ../$(OCL)/
CLSTRING   := $(CLSOURCE:.cl=.clstr)
# The list of non-existing .clstr files
MISSINGCLS := $(filter-out $(wildcard *.clstr),$(CLSTRING))
# This command performs transformation of source CL file into a C string. Call to sed quotes special characters '#' and
# '\', then encloses each line in "...\n".
CLSCOMMAND := sed -e 's/["\\]/\\&/g' -e 's/.*/"&\\n"/'
# Non-existing .clstr should be created right now, since otherwise automatic generation of dependencies will fail. The
# implementation is equivalent to the one used in the make pattern in common.mk (which is used for updating .clstr files
# when necessary), but more complicated since vpath (for .cl files) and automatic pattern substitution cannot be used
$(foreach file,$(MISSINGCLS),$(shell $(CLSCOMMAND) $(patsubst %.clstr,$(PARENT)/%.cl,$(file)) > $(file)))
vpath %.cl $(PARENT)

PROG   := $(PROGOCL)
MYCC   := $(CC)
MYCF   := $(CF)
MYCCPP := $(CCPP)

include $(COMMONMK)
