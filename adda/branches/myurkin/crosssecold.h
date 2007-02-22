int n_non_void_sites(int mat_count[],int Nmat);

double Ext_cross(dcomplex *x,dcomplex *Einc,double k);
double Abs_cross(dcomplex *x,double k);
double Sca_diff(dcomplex *x,dcomplex *r,double k);

void set_Parms(void);
void fill_tab(void);
void finish_int(void);

void calc_alldir(dcomplex *x,REAL **rdip,double k);
void Sca_cross(double k,double *res);
void Asym_parm(double k,double vec[]);
void Asym_parm_x(double k,double *vec);
void Asym_parm_y(double k,double *vec);
void Asym_parm_z(double k,double *vec);

void Frp_mat(double Fsca_tot[3],REAL *Fsca,
	     double Finc_tot[3],REAL *Finc,
	     double Frp_tot[3],REAL *Frp,
	     dcomplex *x,REAL **rdip,double k);
