# Makefile for MPI version of ADDA
# Most options are defined in Makefile
#
# Copyright (C) ADDA contributors
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

# Heurestics to determine if MPI wrapper is available. If $(MPICC) is defined it is used as is. Otherwise we try to
# locate 'mpicc' wrapper through system call, and reverse to $(CC) as the last resort. This should work fine in many
# cases, but can be replaced by a simple assignment of MPICC.
ifndef MPICC
  ifeq ($(shell which mpicc > /dev/null 2>&1 && echo 0),0)
    MPICC := mpicc
  else
    MPICC := $(CC)
  endif
endif

# Compiler dependent options
# These are options for a particular Alpha system
ifeq ($(COMPILER),compaq)
  MPICC := cc
  LDLIBS += -lmpi -lelan
endif

# If the compiler is used directly, a few additional options are needed
ifeq ($(MPICC),$(CC))
  # Microsoft MPI on Windows is detected automatically. If you have it installed, but want to use another 
  # implementation, delete the last symbol on the following line (then this case will not be satisfied)
  ifdef MSMPI_INC
    $(info Using Microsoft MPI)
    # Microsoft MPI does not properly define MSMPI_VER in mpi.h, so we have to get it from registry and then
    # transform it to hex format assumed in mpi.h (only two numbers are used, the rest have to be discarded.
    # If reg fails, we discard its stderr, but identify empty output.
    REG_KEY := $(shell reg.exe query "HKLM\\SOFTWARE\\Microsoft\\MPI" /v Version 2>/dev/null \
                       | awk '/Version/ {print $$NF}')
    ifneq ($(REG_KEY),)
      REG_KEY_SP := $(subst ., ,$(REG_KEY))
      MSMPI_MAJOR := $(firstword $(REG_KEY_SP))
      MSMPI_MINOR := $(word 2,$(REG_KEY_SP))
      MSMPI_VER := $(shell printf "0x%X" $$((256*$(MSMPI_MAJOR) + $(MSMPI_MINOR))))
      CFLAGS += -DMSMPI_VER='$(MSMPI_VER)'
    else
      $(warning Microsoft MPI version was not found in the registry. Hence the corresponding value \
                indicated by ADDA will probably be incorrect)
    endif
    
    CFLAGS += -I'$(MSMPI_INC)'
    LDFLAGS += -L'$(MSMPI_LIB64)'
    LDLIBS += -lmsmpi
  else
    # Depending on a particular MPI installation one may need to manually specify paths to MPI headers and libraries.
    # Path should be absolute or relative to the location of this Makefile.
    #CFLAGS += -I"C:/Program Files/MPICH2/include"
    #LDFLAGS += -L"C:/Program Files/MPICH2/lib"
    LDLIBS += -lmpi
  endif
# If MPI compiler wrapper is used, environmental variables are set to define C compiler and linker
else
  # for MPICH
  export MPICH_CC = $(CC)
  export MPICH_CLINKER = $(CC)
  # for Open MPI
  export OMPI_CC = $(CC)
endif

#=======================================================================================================================
# !!! End of control section. Everything below is not designed to be modified by user
#=======================================================================================================================

# Finalize flags
CFLAGS  += -DADDA_MPI
CSOURCE += matvec.c

PROG   := $(PROGMPI)
MYCC   := $(MPICC)
MYCF   := $(CF)
MYCCPP := $(CCPP)

include $(COMMONMK)
