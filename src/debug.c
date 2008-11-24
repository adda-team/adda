/* FILE: debug.c
 * AUTH: vesseur
 * DESCR: Functions for printing debugging information when compiling 
 *        with option -DDEBUG
 *
 *        Currently is developed by Maxim Yurkin
 */
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include "cmplx.h"
#include "types.h"
#include "comm.h"

/*============================================================*/

void debug_printf (char *fname, int line, char *fmt, ... )
{ 
  va_list args;
  char pr_line[255];
  extern int ringid;
  
  va_start(args, fmt);
  sprintf(pr_line, "%s:%d: ", fname, line);
  vsprintf(pr_line+strlen(pr_line), fmt, args);
  strcat(pr_line, "\n");
/*  assert(isasync(stderr)); */
  printf("(ringID=%i) %s",ringid,pr_line);
  fflush(stdout);
  va_end(args);
}

/*=======================================================*/

void FieldPrint (doublecomplex *x) /* print current field at certain dipole
				not used; left for deep debug */
{
  extern FILE *logfile;
  extern int *material;
  extern int Nmat;
  extern double **DipoleCoord;
  int i,j=0;
  for (i=0;i<local_Ndip && j<9810;++i)
     if (material[i]<Nmat-1) j++;
  fprintf(logfile,"Dipole ccords = %.7e, %.7e, %.7e\n",DipoleCoord[i][0],DipoleCoord[i][1],DipoleCoord[i][2]);
  i*=3;
  fprintf(logfile,"E = %.7e+%.7ei,  %.7e+%.7ei, %.7e+%.7ei\n",
	x[i].r, x[i].i,x[i+1].r, x[i+1].i,x[i+2].r, x[i+2].i);
}
