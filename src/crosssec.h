/* File: crosssec.h
 * $Date::                            $
 * Descr: definitions of functions for calculation of different measured quantities
 *
 * Copyright (C) 2006,2008,2010 ADDA contributors
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
#ifndef __crosssec_h
#define __crosssec_h

void CalcField(doublecomplex * restrict ebuff,const double * restrict n);
void InitRotation(void);
double ExtCross(const double * restrict incPol);
double AbsCross(void);
double ScaCross(char * restrict f_suf);
void ReadAlldirParms(const char * restrict fname);
void ReadAvgParms(const char * restrict fname);
void ReadScatGridParms(const char * restrict fname);
void CalcAlldir(void);
void CalcScatGrid(enum incpol which);
void AsymParm(double *vec,char * restrict f_suf);
void AsymParm_x(double *vec,char * restrict f_suf);
void AsymParm_y(double *vec,char * restrict f_suf);
void AsymParm_z(double *vec,char * restrict f_suf);
void Frp_mat(double Fsca_tot[3],double * restrict Fsca,double Finc_tot[3],double * restrict Finc,
	double Frp_tot[3],double * restrict Frp);

#endif // __crosssec_h
