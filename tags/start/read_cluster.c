#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cmplx.h"
#include "read_cluster.h"

REAL scale_fac=0.25;

void format(/* Change format of Carstens input-file
	       to \t`a`\t"`b`\t`c`\t`d`\t`e`\t'f'.
	       a : particle-label
	       b : material-type-label
	       c : particle-radius
	       d : x-coordinate
	       e : y-coordinate
	       f : z-coordinate
	       Store the result in output-file */
	    char input[],char output[])
{
  int c[2];
  FILE *in,*out;

  in = fopen(input,"r");
  out = fopen(output,"w");

  while ((c[1]=getc(in)) != EOF)
    {
      if (c[1] == ' ')
	;
      else
	{
	  if (c[0] == ' ')
	    putc('\t',out);

	  putc(c[1],out);
	}

      c[0] = c[1];
    }

  fclose(in);
  fclose(out);
}

void read_cluster(/* Read radius, material-type and 
		     grain-coordinates from file 'in' */
		  char input[],
		  int *ngrains,
		  REAL *coords[3],
		  REAL **radius,
		  int **grain_material)
{
  int 
    i,
    label;
  FILE *in;

#if defined(SINGLE)
  char format[] = "\t%d\t%d\t%g\t%g\t%g\t%g\n";
#elif defined(DOUBLE)
  char format[] = "\t%d\t%d\t%lg\t%lg\t%lg\t%lg\n";
#endif

  in = fopen(input,"r");

  fscanf(in,"%d\n",ngrains);

  for (i=0;i<3;++i)
    coords[i] = (REAL *) calloc((*ngrains),sizeof(REAL));
  *grain_material = (int *) calloc((*ngrains),sizeof(int));
  *radius = (REAL *) calloc((*ngrains),sizeof(REAL));

  for (i=0;i<*ngrains;++i)
    {
      fscanf(in,format,
	     &label,&(*grain_material)[i],&(*radius)[i],
	     &coords[0][i],&coords[1][i],&coords[2][i]);
      (*radius)[i] *= scale_fac;
      coords[0][i] *= scale_fac;
      coords[1][i] *= scale_fac;
      coords[2][i] *= scale_fac;
    }

  fclose(in);
}

void setup_box(/* Determine the size of the box
		  circumscribing the aggregate */
	       int ngrains,
	       REAL *coords[3],
	       REAL *max) 
{
  int
    j,comp;
  REAL
    dummy=0;

  /* Search diagonal corners of the rectangular box 
     enclosing the aggregate */
  *max = 0;
  for (j=0;j<ngrains;++j)
    {
      for (comp=0;comp<3;++comp)
	*max = MAX(*max,fabs(coords[comp][j]));
    }
}

void free_aggregate(REAL *coords[3],REAL *radius,int *grain_material)
{
  free(grain_material);
  free(radius);
  free(coords[0]);
  free(coords[1]);
  free(coords[2]);
}
