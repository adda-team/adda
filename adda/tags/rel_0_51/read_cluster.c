/* FILE: read_cluster.c
 * AUTH: Maxim Yurkin
 * DESCR: Reading of clusters (dipole coordinates) from file
 *
 *        reading of old format was implemented by Michel Grimminck
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cmplx.h"
#include "read_cluster.h"
#include "types.h"
#include "comm.h"
#include "const.h"

double scale_fac=0.25;
extern int jagged;

/*============================================================*/

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

  if ((in=fopen(input,"r"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Could not open file %s",input);
  if ((out=fopen(output,"w"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Could not open file %s",output);

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

/*=====================================================*/

void read_cluster(/* Read radius, material-type and 
		     grain-coordinates from file 'in' */
		  char input[],
		  int *ngrains,
		  double *coords[3],
		  double **radius,
		  int **grain_material)
{
  int 
    i,
    label;
  FILE *in;

  char format[] = "\t%d\t%d\t%lg\t%lg\t%lg\t%lg\n";

  if ((in=fopen(input,"r"))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not open file %s",input);

  fscanf(in,"%d\n",ngrains);

  for (i=0;i<3;++i)
    coords[i] = (double *) calloc((*ngrains),sizeof(double));
  *grain_material = (int *) calloc((*ngrains),sizeof(int));
  *radius = (double *) calloc((*ngrains),sizeof(double));

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

/*=======================================================*/

void setup_box(/* Determine the size of the box
		  circumscribing the aggregate */
	       int ngrains,
	       double *coords[3],
	       double *max) 
{
  int
    j,comp;
  double
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

/*==========================================================*/

void free_aggregate(double *coords[3],double *radius,int *grain_material)
{
  free(grain_material);
  free(radius);
  free(coords[0]);
  free(coords[1]);
  free(coords[2]);
}

/*===========================================================*/

void InitDipFile(char *filename,int *maxX, int *maxY, int *maxZ)
   /* read dipole file first to determine box sizes */
{
  FILE *input;
  char format[] = "%d   %d   %d\n";
  int x, y, z;

  if ((input=fopen(filename,"r"))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not open dipole file - %s",filename);
  
  *maxX=*maxY=*maxZ=0;
  while(!feof(input)) {
    if (fscanf(input,format,&x,&y,&z)!=3)
      LogError(EC_ERROR,ONE,POSIT,"Could not scan from dipole file - %s - offset=%d",filename,ftell(input));
    
    if (x>*maxX) *maxX=x;
    if (y>*maxY) *maxY=y;
    if (z>*maxZ) *maxZ=z;
  }
  *maxX=jagged*(*maxX+1);
  *maxY=jagged*(*maxY+1);
  *maxZ=jagged*(*maxZ+1);
  fclose(input);
}

/*===========================================================*/

void ReadDipFile(char *filename)
   /* read dipole file first to determine box sizes */
{
  FILE *input;
  char format[] = "%d   %d   %d\n";
  int x, y, z, x0, y0, z0;
  int index;
  
  extern int *material;
  extern int boxX, boxY;
 
  if ((input=fopen(filename,"r"))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not open dipole file - %s",filename);
  
  while(!feof(input)) {
    if (fscanf(input,format,&x0,&y0,&z0)!=3)
      LogError(EC_ERROR,ONE,POSIT,"Could not scan from dipole file - %s - offset=%d",filename,ftell(input));
    
    for (z=jagged*z0;z<jagged*(z0+1);z++) if (z>=local_z0 && z<local_z1) 
      for (x=jagged*x0;x<jagged*(x0+1);x++) for (y=jagged*y0;y<jagged*(y0+1);y++) {
        index=(z-local_z0)*boxX*boxY+y*boxX+x;
        material[index]=0;
    }
  }
  fclose(input);
}
