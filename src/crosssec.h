/* File: crosssec.h
 * $Date::                            $
 * Descr: definitions of functions for calculation of different measured quantities
 *
 * Copyright (C) 2006,2008,2010,2012-2013 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __crosssec_h
#define __crosssec_h

// project headers
#include "const.h" // for enum definition
#include "types.h" // for doublecomplex

void CalcField(doublecomplex ebuff[static restrict 3],const double n[static restrict 3]);
void InitRotation(void);
double ExtCross(const double * restrict incPol);
double AbsCross(void);
double DecayCross(void);
double ScaCross(const char *f_suf);
void ReadAlldirParms(const char * restrict fname);
void ReadAvgParms(const char * restrict fname);
void ReadScatGridParms(const char * restrict fname);
void SetScatPlane(const double ct,const double st,const double phi,double robs[static restrict 3],
	double polPer[static restrict 3]);
void CalcAlldir(void);
void CalcScatGrid(enum incpol which);
void AsymParm_x(double *vec,const char *f_suf);
void AsymParm_y(double *vec,const char *f_suf);
void AsymParm_z(double *vec,const char *f_suf);
void Frp_mat(double Finc_tot[static restrict 3],double Fsca_tot[static restrict 3],double * restrict Frp);

#endif // __crosssec_h
