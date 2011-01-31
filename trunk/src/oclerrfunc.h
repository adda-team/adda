/* File: oclerrfunc.h
 * $Date::                            $
 * Descr: a Function to get useful information to error codes of OpenCL
 *
 * Copyright (C) 2010-2011 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */


#ifndef __oclerrfunc_h
#define __oclerrfunc_h
char *print_cl_errstring(cl_int err);
void checkErr(cl_int err, const char * name);
#endif
