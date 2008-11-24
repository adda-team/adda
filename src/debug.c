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
  extern double *DipoleCoord;
  int i=9810;

  i*=3;
  fprintf(logfile,"Dipole ccords = %.7e, %.7e, %.7e\n",DipoleCoord[i],DipoleCoord[i+1],DipoleCoord[i+2]);
  fprintf(logfile,"E = %.7e+%.7ei,  %.7e+%.7ei, %.7e+%.7ei\n",
	x[i][re], x[i][im],x[i+1][re], x[i+1][im],x[i+2][re], x[i+2][im]);
}
