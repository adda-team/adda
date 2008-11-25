/* File: crosssec.h
 * $Author$
 * $Date::                            $
 * Descr: definitions of functions for calculation of different measured quantities
 *
 * Copyright (C) 2006,2008 University of Amsterdam
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

#include "function.h" // for function attributes

void CalcField(doublecomplex *ebuff,const double *n);
void InitRotation(void);
double ExtCross(const double *incPol) ATT_PURE;
double AbsCross(void) ATT_PURE;
double ScaCross(char *f_suf);
void ReadAlldirParms(const char *fname);
void ReadAvgParms(const char *fname);
void ReadScatGridParms(const char *fname);
void CalcAlldir(void);
void CalcScatGrid(char which);
void AsymParm(double *vec,char *f_suf);
void AsymParm_x(double *vec,char *f_suf);
void AsymParm_y(double *vec,char *f_suf);
void AsymParm_z(double *vec,char *f_suf);
void Frp_mat(double Fsca_tot[3],double *Fsca,double Finc_tot[3],
             double *Finc,double Frp_tot[3],double *Frp);

#endif // __crosssec_h
