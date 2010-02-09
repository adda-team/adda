# Makefile for sequential version of ADDA package
# All options are defined in Makefile
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

# !!! This file do not have any options designed to be changed by ADDA user

vpath %.c $(CPATH)
vpath %.h $(HPATH)
vpath %.f $(FPATH)

#===================================================================================================
# Main action part
#===================================================================================================

.DELETE_ON_ERROR:

$(PROGSEQ): $(COBJECTS) $(FOBJECTS) 
	@echo "Building $@"
	$(CC) -o $@ $(LDFLAGS) $(COBJECTS) $(FOBJECTS) $(LDLIBS)

# Everything is recompiled when any of makefiles is changed
$(COBJECTS): %.o: %.c %.d $(MFILES)
	$(CC) -c $(CFLAGS) $<

$(FOBJECTS): %.o: %.f $(MFILES)
	$(CF) -c $(FFLAGS) $<

# Dependencies are only generated for C sources; we assume that each Fortran file is completely 
# independent or all of them are compiled at once

$(CDEPEND): %.d: %.c $(MFILES)
	$(CC) $(DEPFLAG) $(CFLAGS) $< $(DFFLAG) $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

include $(CDEPEND)
