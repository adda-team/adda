#define MIN(A,B) (((A) > (B)) ? (B) : (A))
#define MAX(A,B) (((A) < (B)) ? (B) : (A))

void format(char input[],char output[]);
void read_cluster(char input[],int *ngrains,
		  REAL *coords[3],REAL **radius,int **grain_material);
void setup_box(int ngrains,REAL *coords[3],REAL *max);
void free_aggregate(REAL *coords[3],REAL *radius,int *grain_material);
