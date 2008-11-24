/* FILE: crosssec.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        calculation of different measured quantities
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __crosssec_h
#define __crosssec_h

void CalcField(doublecomplex *ebuff,double *n);
void InitRotation(void);
double ExtCross(double *incPol);
double AbsCross(void);
double ScaCross(void);
void ReadAlldirParms(char *fname);
void ReadAvgParms(char *fname);
void ReadScatGridParms(char *fname);
void CalcAlldir(void);
void CalcScatGrid(char which);
void AsymParm(double *vec);
void AsymParm_x(double *vec);
void AsymParm_y(double *vec);
void AsymParm_z(double *vec);
void Frp_mat(double Fsca_tot[3],double *Fsca,
	     double Finc_tot[3],double *Finc,
	     double Frp_tot[3],double *Frp);

#endif /*__crosssec_h*/
