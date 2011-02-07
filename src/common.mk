# Common part of makefiles for different versions of ADDA package
# All options are defined in Makefile and specific makefiles
# $Date::                            $
#
# Copyright (C) 2010-2011 ADDA contributors
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

# !!! This file do not have any options designed to be changed by ADDA user

# define objects and dependencies
COBJECTS   := $(CSOURCE:.c=.o)
CDEPEND    := $(CSOURCE:.c=.d)
FOBJECTS   := $(FSOURCE:.f=.o)
CPPOBJECTS := $(CPPSOURCE:.cpp=.o)

# if Fortran sources are used, Fortran libraries are added
ifneq ($(strip $(FSOURCE)),)
  LDLIBS += $(FLIBS)
endif
# if C++ sources are used, C++ libraries are added
ifneq ($(strip $(CPPSOURCE)),)
  LDLIBS += $(CPPLIBS)
endif
# finalize LDFLAGS, the rest was done in the main Makefile
LDFLAGS  += $(LDLIBS)

# The following are used to track whether recompilation of corresponding parts is required
LDCMD  := $(MYCC) $(LDFLAGS)
CCMD   := $(MYCC) $(CFLAGS)
FCMD   := $(MYCF) $(FFLAGS)
CPPCMD := $(MYCCPP) $(CPPFLAGS)

READ_FILE = $(shell if [ -f $(1) ]; then cat $(1); fi)

ifneq ($(call READ_FILE,$(LDOPTSFILE)),$(LDCMD))
  $(shell rm -f $(LDOPTSFILE))
endif
ifneq ($(call READ_FILE,$(COPTSFILE)),$(CCMD))
  $(shell rm -f $(COPTSFILE))
endif
ifneq ($(call READ_FILE,$(FOPTSFILE)),$(FCMD))
  $(shell rm -f $(FOPTSFILE))
endif

vpath %.c $(CPATH)
vpath %.cpp $(CPPPATH)
vpath %.h $(HPATH)
vpath %.f $(FPATH)

#===================================================================================================
# Main action part
#===================================================================================================

.DELETE_ON_ERROR:

# C linker is used, while Fortran and C++ is handled by addition of libraries for linking
$(PROG): $(COBJECTS) $(FOBJECTS) $(CPPOBJECTS) $(LDOPTSFILE)
	@echo "Building $@"
	$(MYCC) -o $@ $(COBJECTS) $(FOBJECTS) $(CPPOBJECTS) $(LDFLAGS)

$(COBJECTS): %.o: %.c %.d
	$(MYCC) -c $(CFLAGS) $<

$(FOBJECTS): %.o: %.f $(FOPTSFILE)
	$(MYCF) -c $(FFLAGS) $<

$(CPPOBJECTS): %.o: %.cpp $(CPPOPTSFILE)
	$(MYCCPP) -c $(CPPFLAGS) $<

# Dependencies are only generated for C sources; we assume that each Fortran or C++ file is 
# completely independent or all of the files from dependent set are compiled at once. 

$(CDEPEND): %.d: %.c $(COPTSFILE)
	$(MYCC) $(DEPFLAG) $(CFLAGS) $< $(DFFLAG) $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

$(LDOPTSFILE):
	@echo Linking needs to be redone
	echo -n "$(LDCMD)" > $@

$(COPTSFILE):
	@echo C sources need to be recompiled
	echo -n "$(CCMD)" > $@

$(FOPTSFILE):
	@echo Fortran sources need to be recompiled
	echo -n "$(FCMD)" > $@

$(CPPOPTSFILE):
	@echo C++ sources need to be recompiled
	echo -n "$(CPPCMD)" > $@

-include $(CDEPEND)
