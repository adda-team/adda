/* FILE: crosssec.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        calculation of different measured quantities
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __crosssec_h
#define __crosssec_h

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
void Frp_mat(double Fsca_tot[3],double *Fsca,
	     double Finc_tot[3],double *Finc,
	     double Frp_tot[3],double *Frp);

#endif /*__crosssec_h*/
