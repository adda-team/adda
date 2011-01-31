# Makefile for ocl version of ADDA package
# All options are defined in Makefile
# $Date::                            $
#
# Copyright (C) 2010-2011 ADDA contributors
# This file is part of ADDA.
#
# ADDA is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
#
# ADDA is distributed in the hope that it will be useful,but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with 
# ADDA. If not, see <http://www.gnu.org/licenses/>.

# !!! This file do not have any options designed to be changed by ADDA user



PROG := $(PROGOCL)
MYADDACOMP := $(CCPP)
MYCC := $(CC)
MYCF := $(CF)

CFLAGS += -DOPENCL

CPPSRCS = fft_execute.cpp fft_setup.cpp fft_kernelstring.cpp
CPPHEADERS = fft_internal.h fft_base_kernels.h clFFT.h
CPPFLAGS = -Wall -Werror -O3
CPPOBJECTS += fft_execute.o fft_setup.o fft_kernelstring.o

COBJECTS+= oclerrfunc.o oclgetdev.o oclvars.o
CDEPEND+= oclerrfunc.d oclgetdev.d oclvars.d

LDFLAGS += -lOpenCL


include $(COMMONMK)

#CPP source files of APPLE FFT
$(CPPOBJECTS): %.o: %.cpp
	$(CCPP) -c $(CPPFLAGS) $<

