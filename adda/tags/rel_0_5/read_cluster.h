#ifndef __read_cluster_h
#define __read_cluster_h

void format(char input[],char output[]);
void read_cluster(char input[],int *ngrains,
		  double *coords[3],double **radius,int **grain_material);
void setup_box(int ngrains,double *coords[3],double *max);
void free_aggregate(double *coords[3],double *radius,int *grain_material);

void InitDipFile(char *filename,int *maxX, int *maxY, int *maxZ);
void ReadDipFile(char *filename);

#endif /*__read_cluster_h*/