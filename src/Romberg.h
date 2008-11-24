#ifndef __Romberg_h
#define __Romberg_h

#define THETA 0
#define PHI   1
#define Narg  2
#define KMIN  3

double Romberg1D(Parms_1D param,int size,double *data,double *ss); 

void Romberg2D(Parms_1D parms_input[Narg],void (*func_input)(int theta,int phi,double *res),
		    int dim_input, double *ss, char *fname);

#endif /*__Romberg_h*/