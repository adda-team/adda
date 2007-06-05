#ifndef __crosssec_h
#define __crosssec_h

void calc_field (doublecomplex *x, doublecomplex *ebuff, double **rdip, 
                  double *n, double kk, int nlocalDip);

int n_non_void_sites(int mat_count[],int Nmat);
void init_rotation (void);

double Ext_cross(doublecomplex *x,double *incPol,double **rdip,double k);
double Abs_cross(doublecomplex *x,double k);

void set_Parms(void);
void fill_tab(void);
void finish_int(void);

void calc_alldir(doublecomplex *x,double **rdip,double k,char which);
void Sca_cross(double k,double *res);
void Asym_parm(double k,double vec[]);
void Asym_parm_x(double k,double *vec);
void Asym_parm_y(double k,double *vec);
void Asym_parm_z(double k,double *vec);

void Frp_mat(double Fsca_tot[3],double *Fsca,
	     double Finc_tot[3],double *Finc,
	     double Frp_tot[3],double *Frp,
	     doublecomplex *x,double **rdip,double k);

void read_avg_parms(char *fname);
	     
#endif /*__crosssec_h*/