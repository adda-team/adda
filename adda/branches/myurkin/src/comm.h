#define LOCAL 0
#define GLOBAL 1
#define GLOB_LOC 2
 
#if (defined(PVM) || defined(MPI))
#define PARALLEL
#endif

extern int local_d0,local_d1,local_Ndip,Ndip;
extern int local_z0,local_z1,local_Nz;
extern int local_xs_unit;

extern int nprocs;

extern int RingId;
extern int MyProcId;
extern int me,mytid;
extern int tids[256];            /* array of task id */
extern int ringid;

unsigned long extime(void);
void LogError (int ErrCode, char *FName, char *Format,...);

void init_comm(int *argc,char ***argv);
void stop(int);
void par_setup(void);
void synchronize(void);

void Bcast_parms(Parms_1D *parms);

void block_transpose(REAL *X);
void accumulate(REAL *,int);
void my_inner_product(double *a);

void all_gather(void *x_from,void *x_to,char data_type[],int n_elem);
void all_gather_dcomplex(void *x_from,void *x_to,int n_elem);
void all_gather_int(void *x_from,void *x_to,int n_elem);
void all_gather_REAL(void *x_from,void *x_to,int n_elem);

#ifdef PARALLEL
#define printz if (ringid==0) printf
#define fprintz if (ringid==0) fprintf
#define sprintz if (ringid==0) sprintf
#define systemz if (ringid==0) system
#define fclosez if (ringid==0) fclose
#define fopenz(a,b) (ringid==0)?fopen(a,b):NULL;
#else
#define fclosez fclose
#define fopenz fopen
#define systemz system
#define printz printf
#define fprintz fprintf
#define sprintz sprintf
#endif

struct RingData {
  int size;
  int id;
};
typedef struct RingData RingData_t;
 
/* readibility defines...
 */
#define PUBLIC
#define PRIVATE  static

#define ALLNODES 0
#define NULLPTR (void *) 0
 
struct nodenv {
  int procnum;
  int nprocs;
  int host;
  int taskid;
};

/**GD************   Global Defines and Data structures   *****************/
 
#define EC_MASK  0xF0000000
#define EC_FATAL 0xE0000000
#define EC_CRIT  0xD0000000
#define EC_ERROR 0xC0000000
#define EC_WARN  0xB0000000
#define EC_DEBUG 0xA0000000
#define EC_INFO  0x90000000
#define EC_MESS  0x80000000


