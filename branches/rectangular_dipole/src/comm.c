/* FILE : comm.c
 * $Date::                            $
 * Descr: incorporates all parallelization related code, so most of it is directly involved in or closely related to
 *        interprocess communication
 *
 * Copyright (C) 2006-2013 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include "const.h" // keep this first
#include "comm.h" // corresponding header
// project headers
#include "cmplx.h"
#include "debug.h"
#include "fft.h"
#include "function.h"
#include "io.h"
#include "memory.h"
#include "parbas.h"
#include "timing.h"
#include "vars.h"
// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef ADDA_MPI
MPI_Datatype mpi_dcomplex,mpi_int3,mpi_double3,mpi_dcomplex3; // combined datatypes
int *recvcounts,*displs; // arrays of size ringid required for AllGather operations
bool displs_init=false;  // whether arrays above are initialized
#endif

/* whether a synchronize call should be performed before parallel timing. It makes communication timing more accurate,
 * but may deteriorate overall performance by introducing unnecessary delays (test showed only slight difference for
 * granule generator)
 */
#define SYNCHRONIZE_TIMING

#ifdef PARALLEL
#ifndef SPARSE

// SEMI-GLOBAL VARIABLES

// defined and allocated in fft.c
extern double * restrict BT_buffer, * restrict BT_rbuffer;
// defined and initialized in timing.c
extern TIME_TYPE Timing_InitDmComm;

// LOCAL VARIABLES

static int Ntrans;         // number of transmissions; used in CalcPartner
static int * restrict gr_comm_size;  // sizes of transmissions for granule generator communications
static int * restrict gr_comm_overl; // shows whether two sequential transmissions overlap
static unsigned char * restrict gr_comm_ob; // buffer for overlaps
static bool * restrict gr_comm_buf;         // buffer for MPI transfers
#endif // !SPARSE


// First several functions are defined only in parallel mode
//======================================================================================================================

static void RecoverCommandLine(int *argc_p,char ***argv_p)
// eliminate all NULL pointers from argv, shift the rest, and adjust argc accordingly. Used in InitComm
{
	int i,j;

	for (i=0,j=0;i<(*argc_p);i++) {
		if ((*argv_p)[i]==NULL) j++;
		else if (j!=0) (*argv_p)[i-j]=(*argv_p)[i];
	}
	(*argc_p)-=j;
}
//======================================================================================================================

#ifndef SPARSE

static inline size_t IndexBlock(const size_t x,const size_t y,const size_t z,const size_t lengthY)
// index block; used in BlockTranspose
{
	return((z*lengthY+y)*gridX+x);
}

//======================================================================================================================

static inline int CalcPartner(const int tran)
/* calculate ringid of partner processor for current transmission; used in BlockTranspose. Many different
 * implementations are possible; the only requirements are
 * 1) f(tran,f(tran,ringid))=ringid
 * 2) f({1,2,Ntrans},ringid)={0,1,Ntrans}\ringid
 * where f=nprocs is equivalent to skipping this transmission (relevant for odd nprocs)
 */
{
	int part;

	if (ringid==0) part=tran;
	else if (ringid==tran) part=0;
	else {
		part=2*tran-ringid;
		if (part<=0) part+=Ntrans;
		else if (part>Ntrans) part-=Ntrans;
	}
	return part;
}
#endif // !SPARSE

//======================================================================================================================

void CatNFiles(const char * restrict dir,const char * restrict tmpl,const char * restrict dest)
/* cat several temporary files (one for each processor, names defined by the template 'tmpl' that should contain %d to
 * be replaced by ringid). Files are located in directory 'dir'. Combined into 'dest' in the same directory. Afterwards
 * temporary files are removed.
 */
{
	int i,c;
	FILE * restrict in,* restrict out;
	char fname_out[MAX_TMP_FNAME],fname_in[MAX_TMP_FNAME];
	size_t shift;

	// produce full path of destination file and open it
	SnprintfErr(ONE_POS,fname_out,MAX_TMP_FNAME,"%s/%s",dir,dest);
	out=FOpenErr(fname_out,"w",ONE_POS);
	for (i=0;i<nprocs;i++) {
		// produce full path of tmp file and open it
		shift=SnprintfErr(ONE_POS,fname_in,MAX_TMP_FNAME,"%s/",dir);
		/* the following will cause warning by GCC if -Wformat-nonliteral (or -Wformat=2) is used, but we do not know
		 * any other convenient way to make this function work for different file name templates.
		 */
		SnprintfShiftErr(ONE_POS,shift,fname_in,MAX_TMP_FNAME,tmpl,i);
		in=FOpenErr(fname_in,"r",ONE_POS);
		// copy file in to out
		while((c=getc(in))!=EOF) putc(c,out);
		// close and remove tmp file
		FCloseErr(in,fname_in,ONE_POS);
		RemoveErr(fname_in,ONE_POS);
	}
	// close destination file
	FCloseErr(out,fname_out,ONE_POS);
}

//======================================================================================================================

#ifdef ADDA_MPI
static MPI_Datatype MPIVarType(var_type type,bool reduce,int *mult)
/* Chooses MPI datatype corresponding to 'type'. When only copying operations are required (like cast, gather, etc.)
 * exact correspondence is possible. But when reduce operations are used ('reduce'=true), only built-in datatypes can be
 * directly used. In this case we emulate more complex datatypes through multiplication of double, and additional
 * variable 'mult' is returned to account for this factor.
 *
 * Reduction of complex numbers is emulated if not supported; C99 implies that this emulation is portable. Anyway, it
 * is only used for old (less modern) MPI implementations.
 */
{
	if (reduce) *mult=1; // default value when direct correspondence is possible
	switch (type) {
		case uchar_type: return MPI_UNSIGNED_CHAR;
		case int_type: return MPI_INT;
		case sizet_type: return MPI_SIZE_T;
		case double_type: return MPI_DOUBLE;
		case cmplx_type:
#ifndef SUPPORT_MPI_COMPLEX_REDUCE
			if (reduce) {
				*mult=2;
				return MPI_DOUBLE;
			}
			else
#endif
			return mpi_dcomplex;
		case int3_type:
			if (reduce) {
				*mult=3;
				return MPI_INT;
			}
			else return mpi_int3;
		case double3_type:
			if (reduce) {
				*mult=3;
				return MPI_DOUBLE;
			}
			else return mpi_double3;
		case cmplx3_type:
			if (reduce) {
#ifdef SUPPORT_MPI_COMPLEX_REDUCE
				*mult=3;
				return mpi_dcomplex;
#else
				*mult=6;
				return MPI_DOUBLE;
#endif
			}
			else return mpi_dcomplex3;
	}
	LogError(ONE_POS,"Variable type %d is not supported",(int)type);
}

//======================================================================================================================

void InitDispls(void)
// initialize arrays recvcounts and displs once, further calls have no effect
{
	if (!displs_init) {
		if (nvoid_Ndip>INT_MAX)
			LogError(ONE_POS,"int overflow in MPI function for number of non-void dipoles (%zu)",nvoid_Ndip);
		MALLOC_VECTOR(recvcounts,int,nprocs,ALL);
		MALLOC_VECTOR(displs,int,nprocs,ALL);
		// !!! TODO: check for overflow of int
		recvcounts[ringid]=local_nvoid_Ndip;
		displs[ringid]=local_nvoid_d0;
		MPI_Allgather(MPI_IN_PLACE,0,MPI_INT,recvcounts,1,MPI_INT,MPI_COMM_WORLD);
		MPI_Allgather(MPI_IN_PLACE,0,MPI_INT,displs,1,MPI_INT,MPI_COMM_WORLD);
		displs_init=true;
	}
}

#endif // ADDA_MPI

//======================================================================================================================

void AllGather(void * restrict x_from,void * restrict x_to,const var_type type,TIME_TYPE *timing)
/* Gather distributed arrays; works for all types; x_from can be NULL then the data from x_to is used (in_place)
 * increments 'timing' (if not NULL) by the time used
 */
{
#ifdef ADDA_MPI
	MPI_Datatype mes_type;
	TIME_TYPE tstart;

	// redundant initialization to remove warnings
	tstart=0;
	if (timing!=NULL) {
#ifdef SYNCHRONIZE_TIMING
		MPI_Barrier(MPI_COMM_WORLD);  // synchronize to get correct timing
#endif
		tstart=GET_TIME();
	}
	InitDispls(); // actually initialization is done only once
	mes_type=MPIVarType(type,false,NULL);
	if (x_from==NULL) MPI_Allgatherv(MPI_IN_PLACE,0,mes_type,x_to,recvcounts,displs,mes_type,MPI_COMM_WORLD);
	else MPI_Allgatherv(x_from,local_nvoid_Ndip,mes_type,x_to,recvcounts,displs,mes_type,MPI_COMM_WORLD);
	if (timing!=NULL) (*timing)+=GET_TIME()-tstart;
#endif
}

#endif // PARALLEL

//======================================================================================================================

void InitComm(int *argc_p UOIP,char ***argv_p UOIP)
// initialize communications in the beginning of the program
{
#ifdef ADDA_MPI
	int ver,subver;

	/* MPI_Init may alter argc and argv and interfere with normal parsing of command line parameters. The way of
	 * altering is implementation depending. MPI searches for MPI parameters in the command line and removes them (we
	 * assume some kind of removing does take place - otherwise ADDA will produce error 'unknown parameter'). The best
	 * would be to change argc and argv so that they look like no special command line arguments are present. However,
	 * MPICH 1.2.5, for example, just replaces corresponding parameters by NULLs. To incorporate it we introduce special
	 * function to restore the command line
	 */
	MPI_Init(argc_p,argv_p);
	tstart_main = GET_TIME(); // initialize program time
	RecoverCommandLine(argc_p,argv_p);
	// initialize ringid and nprocs
	MPI_Comm_rank(MPI_COMM_WORLD,&ringid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
#ifndef SPARSE
	// initialize Ntrans
	if (IS_EVEN(nprocs)) Ntrans=nprocs-1;
	else Ntrans=nprocs;
#endif
	// define a few derived datatypes
#ifdef SUPPORT_MPI_COMPLEX
	mpi_dcomplex = MPI_C_DOUBLE_COMPLEX; // use built-in datatype if supported
#else
	MPI_Type_contiguous(2,MPI_DOUBLE,&mpi_dcomplex);
	MPI_Type_commit(&mpi_dcomplex);
#endif
	MPI_Type_contiguous(3,MPI_INT,&mpi_int3);
	MPI_Type_commit(&mpi_int3);
	MPI_Type_contiguous(3,MPI_DOUBLE,&mpi_double3);
	MPI_Type_commit(&mpi_double3);
	MPI_Type_contiguous(3,mpi_dcomplex,&mpi_dcomplex3);
	MPI_Type_commit(&mpi_dcomplex3);
	// check MPI version at runtime
	MPI_Get_version(&ver,&subver);
	if (!GREATER_EQ2(ver,subver,RUN_MPI_VER_REQ,RUN_MPI_SUBVER_REQ)) {
		if (GREATER_EQ2(ver,subver,MPI_VER_REQ,MPI_SUBVER_REQ)) // library fits into minimum requirements
			LogError(ONE_POS,"MPI library version (%d.%d) is too old for current ADDA executable. Version %d.%d or "
				"newer is required. Alternatively, you may recompile ADDA using this version of the library.",
				ver,subver,RUN_MPI_VER_REQ,RUN_MPI_SUBVER_REQ);
		else if (RUN_MPI_VER_REQ==MPI_VER_REQ && RUN_MPI_SUBVER_REQ==MPI_SUBVER_REQ) // executable requires minimum
			LogError(ONE_POS,"MPI library version (%d.%d) is too old. Version %d.%d or newer is required",
				ver,subver,RUN_MPI_VER_REQ,RUN_MPI_SUBVER_REQ);
		else // two upgrade options
			LogError(ONE_POS,"MPI library version (%d.%d) is too old. Version %d.%d or newer is required for the "
				"current ADDA executable or at least %d.%d, if ADDA is recompiled using such library.",
				ver,subver,RUN_MPI_VER_REQ,RUN_MPI_SUBVER_REQ,MPI_VER_REQ,MPI_SUBVER_REQ);
	}
	D("MPI library version: %d.%d",ver,subver);
	// if MPI crashes, it happens here
	Synchronize();
#elif !defined(PARALLEL)
	nprocs=1;
	ringid=ADDA_ROOT;
#endif

#ifndef SPARSE // CheckNprocs does not exist in sparse mode
	// check if weird number of processors is specified; called even in sequential mode to initialize weird_nprocs
	CheckNprocs();
#endif // !SPARSE
}

//======================================================================================================================

void Stop(const int code)
// stops the program with exit 'code'
{
#ifdef ADDA_MPI
	if (code!=EXIT_SUCCESS) { // error occurred
		fflush(stdout);
		fprintf(stderr,"Aborting process %d\n",ringid);
		fflush(stderr);
		MPI_Abort(MPI_COMM_WORLD,code);
	}
	else { // regular termination
		// clean MPI constructs and some memory
#ifndef SUPPORT_MPI_COMPLEX
		MPI_Type_free(&mpi_dcomplex);
#endif
		MPI_Type_free(&mpi_int3);
		MPI_Type_free(&mpi_double3);
		MPI_Type_free(&mpi_dcomplex3);
		if (displs_init) {
			Free_general(recvcounts);
			Free_general(displs);
		}
		// wait for all processors
		fflush(stdout);
		Synchronize();
		// finalize MPI communications
		MPI_Finalize();
	}
#endif
	exit(code);
}

//======================================================================================================================

void Synchronize(void)
// synchronizes all processes
{
#ifdef ADDA_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//======================================================================================================================

void MyBcast(void * restrict data UOIP,const var_type type UOIP,const size_t n_elem UOIP,TIME_TYPE *timing UOIP)
/* casts values stored in '*data' from root processor to all other; works for all types;
 * increments 'timing' (if not NULL) by the time used
 */
{
#ifdef ADDA_MPI
	TIME_TYPE tstart=0; // redundant initialization to remove warnings

	if (n_elem>INT_MAX) LogError(ONE_POS,"int overflow in MPI function (%zu)",n_elem);
	if (timing!=NULL) {
#ifdef SYNCHRONIZE_TIMING
		MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
		tstart=GET_TIME();
	}
	MPI_Bcast(data,n_elem,MPIVarType(type,false,NULL),ADDA_ROOT,MPI_COMM_WORLD);
	if (timing!=NULL) (*timing)+=GET_TIME()-tstart;
#endif
}

//======================================================================================================================

void BcastOrient(int *i UOIP, int *j UOIP, int *k UOIP)
/* cast current orientation angle (in orientation averaging) to all processes from root. This can be done as atomic
 * operation using some complex MPI datatype, but it is not worth it.
 */
{
#ifdef ADDA_MPI
	int buf[3];

	if (IFROOT) {
		buf[0]=*i;
		buf[1]=*j;
		buf[2]=*k;
	}
	MPI_Bcast(buf,3,MPI_INT,ADDA_ROOT,MPI_COMM_WORLD);
	if (!IFROOT) {
		*i=buf[0];
		*j=buf[1];
		*k=buf[2];
	}
#endif
}

//======================================================================================================================

double AccumulateMax(double data UOIP,double *max UOIP)
// given a single double on each processor, accumulates their sum (returns) and maximum on root processor
{
#ifdef ADDA_MPI
	double buf;
	// potentially can be optimized by combining into one operation
	MPI_Reduce(&data,&buf,1,MPI_DOUBLE,MPI_SUM,ADDA_ROOT,MPI_COMM_WORLD);
	MPI_Reduce(&data,max,1,MPI_DOUBLE,MPI_MAX,ADDA_ROOT,MPI_COMM_WORLD);
	return buf;
#else
	return data;
#endif
}

//======================================================================================================================

void Accumulate(void * restrict data UOIP,const var_type type UOIP,size_t n UOIP,TIME_TYPE *timing UOIP)
// Gather and add complex vector on processor root; total time is saved in timing (NOT incremented).
// Can be easily made to accept any variable type
{
#ifdef ADDA_MPI
	MPI_Datatype mes_type;
	int mult;
	TIME_TYPE tstart;

	if (n>INT_MAX) LogError(ONE_POS,"int overflow in MPI function (%zu)",n);
#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	mes_type=MPIVarType(type,true,&mult);
	n*=mult;
	// Strange, but MPI 2.2 doesn't seem to support calling the following the same way on all processes
	if (IFROOT) MPI_Reduce(MPI_IN_PLACE,data,n,mes_type,MPI_SUM,ADDA_ROOT,MPI_COMM_WORLD);
	else MPI_Reduce(data,NULL,n,mes_type,MPI_SUM,ADDA_ROOT,MPI_COMM_WORLD);
	(*timing)=GET_TIME()-tstart;
#endif
}

//======================================================================================================================

void MyInnerProduct(void * restrict data UOIP,const var_type type UOIP,size_t n UOIP,TIME_TYPE *timing UOIP)
/* gather values stored in *data, add them and return them in *data; works for all types; increments 'timing' (if not
 * NULL) by the time used;
 * Similar to Accumulate(), but returns the result to _all_ processors, and timing is _incremented_
 */
{
#ifdef ADDA_MPI
	MPI_Datatype mes_type;
	int mult;
	TIME_TYPE tstart=0; // redundant initialization to remove warnings

	if (n>INT_MAX) LogError(ONE_POS,"int overflow in MPI function (%zu)",n);
	if (timing!=NULL) {
#ifdef SYNCHRONIZE_TIMING
		MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
		tstart=GET_TIME();
	}
	mes_type=MPIVarType(type,true,&mult);
	n*=mult;
	MPI_Allreduce(MPI_IN_PLACE,data,n,mes_type,MPI_SUM,MPI_COMM_WORLD);
	if (timing!=NULL) (*timing)+=GET_TIME()-tstart;
#endif
}

//======================================================================================================================

void ParSetup(void)
// initialize common parameters; need to do in the beginning to enable call to MakeParticle
{
#ifndef SPARSE // FFT mode initialization
#	ifdef PARALLEL
	int unitZ,unitX;
#	endif
	// calculate size of 3D grid
	gridX=fftFit(2*boxX,nprocs);
	gridY=fftFit(2*boxY,1);
	gridZ=fftFit(2*boxZ,2*nprocs);
	// initialize some variables
	smallY=gridY/2;
	smallZ=gridZ/2;
	/* if this check is passed then all other multiplications of 2 grids are OK, except for XY values, used in granule
	 * generator
	 */
	gridYZ=MultOverflow(gridY,gridZ,ALL_POS,"gridYZ");
#	ifdef PARALLEL
	unitZ=smallZ/nprocs; // this should always be an exact division
	local_z0=ringid*unitZ;
	local_z1=(ringid+1)*unitZ;
	if (local_z1 > boxZ) local_z1_coer=boxZ;
	else local_z1_coer=local_z1;
	unitX=gridX/nprocs;
	local_x0=ringid*unitX;
	local_x1=(ringid+1)*unitX;
#	else
	local_z0=0;
	local_z1=smallZ;
	local_z1_coer=boxZ;
	local_x0=0;
	local_x1=gridX;
#	endif
	if (local_z1_coer<=local_z0) {
		LogWarning(EC_INFO,ALL_POS,"No real dipoles are assigned");
		local_z1_coer=local_z0;
	}
	local_Nz=local_z1-local_z0;
	local_Nx=local_x1-local_x0;
	boxXY=boxX*(size_t)boxY; // overflow check is covered by gridYZ above
	local_Ndip=MultOverflow(boxXY,local_z1_coer-local_z0,ALL_POS,"local_Ndip");
	D("%i :  %i %i %i %zu %zu \n",ringid,local_z0,local_z1_coer,local_z1,local_Ndip,local_Nx);
#else // SPARSE
	/* For sparse mode, nvoid_Ndip is defined in InitDipFile(), and here we define local_nvoid_d0 and local_nvoid_d1,
	 * since they are required already in ReadDipFile()
	 */
	D("%zu",nvoid_Ndip);
#	ifdef PARALLEL
	double unit_f=1.0/nprocs;
	//making sure that the first local_f0 and last local_f1 are 0.0 and 1.0
	double local_f0 = (ringid == 0) ? 0.0 : ringid*unit_f;
	double local_f1 = (ringid == nprocs-1) ? 1.0 : (ringid+1)*unit_f;
	local_nvoid_d0 = local_f0*nvoid_Ndip;
	local_nvoid_d1 = local_f1*nvoid_Ndip;
#	else
	local_nvoid_d0=0;
	local_nvoid_d1=nvoid_Ndip;
#	endif // PARALLEL
#endif
}

//======================================================================================================================

void SetupLocalD(void)
// initialize local_nvoid_d0 and local_nvoid_d1 in FFT mode. In sparse mode those are set earlier in ParSetup()
{
#ifdef ADDA_MPI
	/* use of exclusive scan (MPI_Exscan) is logically more suitable, but it has special behavior for the ringid=0. The
	 * latter would require special additional arrangements.
	 */
	MPI_Scan(&local_nvoid_Ndip,&local_nvoid_d1,1,MPI_SIZE_T,MPI_SUM,MPI_COMM_WORLD);
	local_nvoid_d0=local_nvoid_d1-local_nvoid_Ndip;
#else
	local_nvoid_d0=0;
	local_nvoid_d1=local_nvoid_Ndip;
#endif
}

//======================================================================================================================

void ReadField(const char * restrict fname,doublecomplex *restrict field)
/* Reads a complex field from file 'fname' and stores into 'field'.
 *
 * In MPI mode the algorithm is very far from optimal, since the whole file is read in total nprocs/2 times.
 * However, this is mainly the limitation of the text file (need to test for blank lines and format consistency).
 * The only feasible way to improve it is to use binary format like NetCDF (which may work on top of MPI_IO)
 */
{
	char linebuf[BUF_LINE];
	TIME_TYPE tstart=GET_TIME();
	FILE *file=FOpenErr(fname,"r",ALL_POS);
	// the same format as used for saving the beam by StoreFields(...) in make_particle.c
	const char format[]="%*f %*f %*f %*f %lf %lf %lf %lf %lf %lf";
	const char test_form[]="%*f"; // a quick test for a blank line, without actual assignment
	const int mustbe=6;
	int scanned;
	size_t i,j;
	double buf[6];

#if defined(ADDA_MPI) && defined(SYNCHRONIZE_TIMING)
	MPI_Barrier(MPI_COMM_WORLD);  // synchronize to get correct timing
#endif
	// skips first line with headers and any comments, if present
	size_t line=SkipNLines(file,1);
	line+=SkipComments(file);
	i=j=0;
	while(FGetsError(file,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL) {
		// scan numbers in a line
		if (i<local_nvoid_d0) { // just count non-blank lines
			if (sscanf(linebuf,test_form)!=EOF) i++;
		}
		else if (i==nvoid_Ndip) { // tests that file doesn't contains extra data rows
			if (sscanf(linebuf,test_form)!=EOF) LogError(ALL_POS,"Field file %s contains more data rows than number of "
				"dipoles (%zu) in the particle",fname,nvoid_Ndip);
		}
		else { // here local_nvoid_d0 <= i < local_nvoid_d1
			scanned=sscanf(linebuf,format,buf,buf+1,buf+2,buf+3,buf+4,buf+5);
			field[j] = buf[0] + I*buf[1];
			field[j+1] = buf[2] + I*buf[3];
			field[j+2] = buf[4] + I*buf[5];
			if (scanned==EOF) {
				/* in most cases EOF indicates blank line (which we just skip), but it can also be a short (ill-format)
				 * line, for which we test below. Otherwise, such lines may be interpreted differently b different
				 * processors
				 */
				if (sscanf(linebuf,test_form)!=EOF)
					LogError(ALL_POS,"Error occurred during scanning of line %zu from field file %s",line,fname);
			}
			else {
				if (scanned!=mustbe) // this in most cases indicates wrong format
					LogError(ALL_POS,"Error occurred during scanning of line %zu from field file %s",line,fname);
				j+=3;
				i++;
				/* all processors stop reading file as soon as possible, but the last processor reads one more line to
				 * test (above) for extra strings in the file
				 */
				if (i==local_nvoid_d1 && i!=nvoid_Ndip) break;
			}
		}
	}
#if defined(ADDA_MPI) && defined(SYNCHRONIZE_TIMING)
	MPI_Barrier(MPI_COMM_WORLD);  // synchronize to get correct timing
#endif
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

#ifndef SPARSE

void BlockTranspose(doublecomplex * restrict X UOIP,TIME_TYPE *timing UOIP)
/* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ; specializes at Xmatrix; do 3 components in one
 * message; increments 'timing' (if not NULL) by the time used
 *
 *  !!! TODO: Although size_t is used for bufsize,etc., MPI functions take int as arguments. This limits the largest
 *  possible size to some extent. Moreover, the size of int is not really well predicted. The exact implications of this
 *  are still unclear.
 */
{
#ifdef ADDA_MPI
	TIME_TYPE tstart;
	size_t bufsize,msize,posit,step,y,z;
	int transmission,part,Xpos,Xcomp;
	MPI_Status status;

	// redundant initialization to remove warnings
	tstart=0;

	if (timing!=NULL) {
#ifdef SYNCHRONIZE_TIMING
		MPI_Barrier(MPI_COMM_WORLD);  // synchronize to get correct timing
#endif
		tstart=GET_TIME();
	}
	step=2*local_Nx;
	msize=local_Nx*sizeof(doublecomplex);
	bufsize=6*local_Nz*smallY*local_Nx;
	if (bufsize>INT_MAX)
		LogError(ALL_POS,"int overflow in MPI function for BT buffer (%zu)",bufsize);

	for(transmission=1;transmission<=Ntrans;transmission++) {
		// if part==nprocs then skip this transmission
		if ((part=CalcPartner(transmission))!=nprocs) {
			posit=0;
			Xpos=local_Nx*part;
			for(Xcomp=0;Xcomp<3;Xcomp++) for(z=0;z<local_Nz;z++) for(y=0;y<smallY;y++) {
				memcpy(BT_buffer+posit,X+Xcomp*local_Nsmall+IndexBlock(Xpos,y,z,smallY),msize);
				posit+=step;
			}

			MPI_Sendrecv(BT_buffer, bufsize, MPI_DOUBLE, part, 0,
				BT_rbuffer, bufsize, MPI_DOUBLE, part, 0,
				MPI_COMM_WORLD,&status);

			posit=0;
			Xpos=local_Nx*part;
			for(Xcomp=0;Xcomp<3;Xcomp++) for(z=0;z<local_Nz;z++) for(y=0;y<smallY;y++) {
				memcpy(X+Xcomp*local_Nsmall+IndexBlock(Xpos,y,z,smallY),BT_rbuffer+posit,msize);
				posit+=step;
			}
		}
	}
	if (timing!=NULL) (*timing)+=GET_TIME()-tstart;
#endif
}

//======================================================================================================================

void BlockTranspose_DRm(doublecomplex * restrict X UOIP,const size_t lengthY UOIP,const size_t lengthZ UOIP)
/* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ; specialized for D or R matrix. It can be
 * updated to accept timing argument for generality. But, since this is a specialized function, we keep the timing
 * variable hard-wired in the code.
 */
{
#ifdef ADDA_MPI
	TIME_TYPE tstart;
	size_t bufsize,msize,posit,step,y,z;
	int transmission,part,Xpos;
	MPI_Status status;

#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	step=2*local_Nx;
	msize=local_Nx*sizeof(doublecomplex);
	bufsize = 2*lengthZ*lengthY*local_Nx;
	if (bufsize>INT_MAX)
		LogError(ALL_POS,"int overflow in MPI function for BT buffer (%zu)",bufsize);

	for(transmission=1;transmission<=Ntrans;transmission++) {
		if ((part=CalcPartner(transmission))!=nprocs) {
			posit=0;
			Xpos=local_Nx*part;
			for(z=0;z<lengthZ;z++) for(y=0;y<lengthY;y++) {
				memcpy(BT_buffer+posit,X+IndexBlock(Xpos,y,z,lengthY),msize);
				posit+=step;
			}

			MPI_Sendrecv(BT_buffer,bufsize,MPI_DOUBLE,part,0,
				BT_rbuffer,bufsize,MPI_DOUBLE,part,0,
				MPI_COMM_WORLD,&status);

			posit=0;
			Xpos=local_Nx*part;
			for(z=0;z<lengthZ;z++) for(y=0;y<lengthY;y++) {
				memcpy(X+IndexBlock(Xpos,y,z,lengthY),BT_rbuffer+posit,msize);
				posit+=step;
			}
		}
	}
	Timing_InitDmComm += GET_TIME() - tstart;
#endif
}

//======================================================================================================================

#ifdef PARALLEL
void CalcLocalGranulGrid(const double z0,const double z1,const double gdZ,const int gZ,const int id,int *lz0,int *lz1)
/* calculates starting and ending (+1) cell of granule grid (lz0 & lz1) on a processor with ringid=id
 */
{
	int dzl,dzh; // similar to local_z0 and local_z1

	dzl=local_Nz*id;
	// should not be coerced because the result differs only for dzh>boxZ, then dzh-1>z0
	dzh=dzl+local_Nz;
	if (dzl>z1) *lz0=*lz1=gZ;
	else {
		if (dzl>z0) *lz0=(int)floor((dzl-z0)/gdZ);
		else *lz0=0;
		if (dzh>z1) *lz1=gZ;
		else if (dzh-1>z0) *lz1=(int)floor((dzh-z0-1)/gdZ)+1;
		else *lz1=0;
	}
}
#endif

//======================================================================================================================

void SetGranulComm(const double z0 UOIP,const double z1 UOIP,const double gdZ UOIP,const int gZ,const size_t gXY UOIP,
	size_t max_gran UOIP,int *lz0,int *lz1,const int sm_gr UOIP)
/* sets communication for granule generator; max_gran - maximum number of granules in one set (used to allocate buffer);
 * sm_gr - whether granules are small (simpler)
 */
{
#ifdef PARALLEL
	int i,loc0,loc1,loc1_prev=0;

	MALLOC_VECTOR(gr_comm_buf,bool,max_gran,ALL);
	if (!sm_gr) {
		if (IFROOT) {
			MALLOC_VECTOR(gr_comm_size,int,nprocs,ONE);
			MALLOC_VECTOR(gr_comm_overl,int,nprocs-1,ONE);
			// always allocated, not to mess with its freeing
			MALLOC_VECTOR(gr_comm_ob,uchar,gXY,ONE);
			/* The following is very inefficient (may be significantly optimized), but using one
			 * common function is more robust.
			 */
			for (i=0;i<nprocs;i++) {
				CalcLocalGranulGrid(z0,z1,gdZ,gZ,i,&loc0,&loc1);
				if (i==ADDA_ROOT) {
					*lz0=loc0;
					*lz1=loc1;
				}
				gr_comm_size[i]=loc1-loc0;
				if (i!=0) gr_comm_overl[i-1]=(loc0<loc1_prev);
				loc1_prev=loc1;;
			}
		}
		else CalcLocalGranulGrid(z0,z1,gdZ,gZ,ringid,lz0,lz1);
	}
#else
	*lz0=0;
	*lz1=gZ;
#endif
}

//======================================================================================================================

void CollectDomainGranul(unsigned char * restrict dom UOIP,const size_t gXY UOIP,const int lz0 UOIP,
	const int locgZ UOIP,TIME_TYPE *timing UOIP)
// collects the map of domain for granule generator on the root processor; timing is incremented by the total time used
{
#ifdef ADDA_MPI
	int i,unit,index;
	size_t j;
	MPI_Status status;
	TIME_TYPE tstart;

#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	unit=gXY*sizeof(char);
	if (IFROOT) {
		index=(lz0+gr_comm_size[ADDA_ROOT])*gXY;
		for (i=ADDA_ROOT+1;i<nprocs;i++) {
			if (gr_comm_size[i]!=0) {
				if (gr_comm_overl[i-1]) {
					index-=gXY;
					memcpy(gr_comm_ob,dom+index,unit);
				}
				MPI_Recv(dom+index,unit*gr_comm_size[i],MPI_UNSIGNED_CHAR,i,0,MPI_COMM_WORLD,&status);
				if (gr_comm_overl[i-1]) for (j=0;j<gXY;j++) dom[index+j]|=gr_comm_ob[j];
				index+=gXY*gr_comm_size[i];
			}
		}
		// that is only needed when ADDA_ROOT!=0; kind of weird but should work
		index=lz0*gXY;
		for (i=ADDA_ROOT-1;i>=0;i--) {
			if (gr_comm_size[i]!=0) {
				if (gr_comm_overl[i]) {
					memcpy(gr_comm_ob,dom+index,unit);
					index+=gXY;
				}
				MPI_Recv(dom+index-gXY*gr_comm_size[i],unit*gr_comm_size[i],MPI_UNSIGNED_CHAR,i,0,MPI_COMM_WORLD,
					&status);
				if (gr_comm_overl[i]) for (j=0;j<gXY;j++) dom[index-gXY+j]|=gr_comm_ob[j];
				index-=gXY*gr_comm_size[i];
			}
		}
	}
	else if (locgZ!=0) {
		// the test here implies the test for above MPI_Recv as well
		size_t size=(size_t)unit*(size_t)locgZ;
		if (size>INT_MAX) LogError(ALL_POS,"int overflow in MPI function (%zu)",size);
		MPI_Send(dom,size,MPI_UNSIGNED_CHAR,ADDA_ROOT,0,MPI_COMM_WORLD);
	}
	(*timing)+=GET_TIME()-tstart;
#endif
}

//======================================================================================================================

void FreeGranulComm(const int sm_gr UOIP)
// frees all additional memory used for communications of granule generator; simpler if small granules
{
#ifdef PARALLEL
	Free_general(gr_comm_buf);
	if (!sm_gr && IFROOT) {
		Free_general(gr_comm_size);
		Free_general(gr_comm_overl);
		Free_general(gr_comm_ob);
	}
#endif
}

//======================================================================================================================

void ExchangeFits(bool * restrict data UOIP,const size_t n UOIP,TIME_TYPE *timing UOIP)
/* performs a collective AND operation on the (vector) data; timing is incremented by the total time used.
 * Starting from MPI 2.2 (with errata) MPI_C_BOOL should occupy 1 byte, the same when substitute is used
 */
{
#ifdef ADDA_MPI
	TIME_TYPE tstart;

	if (n>INT_MAX) LogError(ONE_POS,"int overflow in MPI function (%zu)",n);
#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	MPI_Allreduce(data,gr_comm_buf,n,mpi_bool,MPI_LAND,MPI_COMM_WORLD);
	memcpy(data,gr_comm_buf,n*sizeof(bool));
	(*timing)+=GET_TIME()-tstart;
#endif
}

//======================================================================================================================

#ifdef PARALLEL

bool ExchangePhaseShifts(doublecomplex * restrict bottom, doublecomplex * restrict top,TIME_TYPE *timing)
/* propagates slice of complex values from bottom to top. In the beginning 'top' contains phase shift over the current
 * processor, at the end 'bottom' and 'top' contain phase shifts from the bottom of the first processor to bottom and
 * top of the current processor. Potentially, can be optimized by some tree algorithm. However, there seem to be no
 * ready MPI function available.
 */
{
#ifdef ADDA_MPI
	MPI_Status status;
	TIME_TYPE tstart;
	size_t i;

	if (2*boxXY>INT_MAX) LogError(ONE_POS,"int overflow in MPI function (%zu)",2*boxXY);
#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	// receive slice from previous processor and increment own slice by these values
	if (ringid>0) { // It is important to use 0 instead of ROOT
		MPI_Recv(bottom,2*boxXY,MPI_DOUBLE,ringid-1,0,MPI_COMM_WORLD,&status);
		for (i=0;i<boxXY;i++) top[i]+=bottom[i];
	}
	// send updated slice to previous processor
	if (ringid<(nprocs-1)) MPI_Send(top,2*boxXY,MPI_DOUBLE,ringid+1,0,MPI_COMM_WORLD);
#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	(*timing)+=GET_TIME()-tstart;
	return (ringid!=0);
#endif
}

#endif // PARALLEL

#endif // !SPARSE
