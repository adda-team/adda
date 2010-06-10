# Makefile for MPI version of ADDA package
# Most options are defined in Makefile
# $Author$
# $Date::                            $
#
# Copyright (C) 2006-2008 University of Amsterdam
# Copyright (C) 2009,2010 Institute of Chemical Kinetics and Combustion & University of Amsterdam
# This file is part of ADDA.
#
# ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along with ADDA. If not, see
# <http://www.gnu.org/licenses/>.

#===================================================================================================
# !!! Start of control section. Flags and options here are designed to be modified by user to choose
# particular options for compilation. However, the default set of options may work out-of-box on
# some systems.
#===================================================================================================

# Heurestics to determine if MPI wrapper is available. If $(MPICC) is defined it is used as is. 
# Otherwise we try to locate 'mpicc' wrapper through system call, and reverse to $(CC) as the last 
# resort. This should work fine in many cases, but can be replaced by a simple assignment of MPICC.
ifndef MPICC
  ifeq ($(shell which mpicc > /dev/null 2>&1 && echo 0),0)
    MPICC = mpicc
  else
    MPICC = $(CC)
  endif
endif

# Compiler dependent options
# These are options for a particular Alpha system
ifeq ($(COMPILER),compaq)
  MPICC = cc
  LDLIBS += -lmpi -lelan
endif

# If the compiler is used directly, a few additional options are needed
ifeq ($(MPICC),$(CC))
  # Depending on a particular MPI installation one may need to manually specify paths to MPI headers
  # and libraries. Path should be absolute or relative to the location of this Makefile.
  #CFLAGS += -I"C:/Program Files/MPICH2/include"
  #LDFLAGS += -L"C:/Program Files/MPICH2/lib"

  LDLIBS += -lmpi
# If MPI compiler wrapper is used, a few environmental variables are set
else
  # These are used by MPICH, other implementations may also need something
  export MPICH_CC = $(CC)
  export MPICH_CLINKER = $(CC)
endif

#===================================================================================================
# !!! End of control section. Everything below is not designed to be modified by user
#===================================================================================================

# Finalize flags
CDEFS += -DADDA_MPI
vpath %.c $(CPATH)
vpath %.h $(HPATH)
vpath %.f $(FPATH)

#===================================================================================================
# Main action part
#===================================================================================================

.DELETE_ON_ERROR:

$(PROGMPI): $(COBJECTS) $(FOBJECTS)
	@echo "Building $@"
	$(MPICC) -o $@ $(LDFLAGS) $(COBJECTS) $(FOBJECTS) $(LDLIBS)

# Everything is recompiled when any of makefiles is changed)
$(COBJECTS): %.o: %.c %.d $(MFILES)
	$(MPICC) -c $(CFLAGS) $(CDEFS) $<

$(FOBJECTS): %.o: %.f $(MFILES)
	$(CF) -c $(FFLAGS) $<

# Dependencies are only generated for C sources; we assume that each Fortran file is completely
# independent or all of them are compiled at once

$(CDEPEND): %.d: %.c $(MFILES)
	$(MPICC) $(DEPFLAG) $(CFLAGS) $(CDEFS) $< $(DFFLAG) $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(CDEPEND)