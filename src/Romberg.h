#define THETA 0
#define PHI   1
#define Narg  2


void nrerror(char error_text[]);
REAL *alloc_Rvector(int nl,int nh);
void freeRv(REAL *v,int nl,int nh);

void polint(REAL xa[],REAL ya[],int n,REAL x,REAL *y,REAL *dy);


void inner_trapzd(Parms_1D input[Narg],Rvector (*FUNC)(int theta,int phi),int fixed,Rvector *s,int n);
Rvector inner_qromb(Parms_1D input[Narg],Rvector (*FUNC)(int theta,int phi),int fixed);
void outer_trapzd(Parms_1D input[Narg],Rvector (*FUNC)(int theta,int phi),Rvector *s,int n);
Rvector outer_qromb(Parms_1D input[Narg],Rvector (*FUNC)(int theta,int phi),int dim,char logfile[]);

