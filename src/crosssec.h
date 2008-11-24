/* FILE: crosssec.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        calculation of different measured quantities
 */
#ifndef __crosssec_h
#define __crosssec_h

void calc_field (doublecomplex *ebuff,double *n);

void init_rotation (void);

double Ext_cross(double *incPol);
double Abs_cross(void);
double Sca_cross(void);

void set_Parms(void);
void free_parms(void); 
void fill_tab(void);
void finish_int(void);

void calc_alldir(char which);
void Asym_parm(double *vec);
void Asym_parm_x(double *vec);
void Asym_parm_y(double *vec);
void Asym_parm_z(double *vec);

void Frp_mat(double Fsca_tot[3],double *Fsca,
	     double Finc_tot[3],double *Finc,
	     double Frp_tot[3],double *Frp);
void read_avg_parms(char *fname);
	     
#endif /*__crosssec_h*/
