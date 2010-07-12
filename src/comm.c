/* FILE : comm.c
 * $Author$
 * $Date::                            $
 * Descr: incorporates all parallelization related code, so most of it is directly involved in or
 *        closely related to interprocess communication
 *
 *        Previous versions were by Martijn Frijlink
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * Copyright (C) 2009 Institute of Chemical Kinetics and Combustion & University of Amsterdam
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vars.h"
#include "comm.h"
#include "const.h"
#include "io.h"
#include "fft.h"
#include "memory.h"
#include "timing.h"
#include "function.h"
#include "parbas.h"

#ifdef ADDA_MPI
MPI_Datatype mpi_dcomplex;
#endif

/* whether a synchronize call should be performed before parallel timing. It makes communication
 * timing more accurate, but may deteriorate overall performance by introducing unnecessary
 * delays (test showed only slight difference for granule generator) */
#define SYNCHRONIZE_TIMING

// SEMI-GLOBAL VARIABLES

// defined and allocated in fft.c
extern double *BT_buffer, *BT_rbuffer;
// defined and initialized in timing.c
extern TIME_TYPE Timing_Dm_Init_comm;

// LOCAL VARIABLES

#ifdef PARALLEL
static int Ntrans;         // number of transmissions; used in CalcPartner
static int *gr_comm_size;  // sizes of transmissions for granule generator communications
static int *gr_comm_overl; // shows whether two sequential transmissions overlap
static unsigned char *gr_comm_ob; // buffer for overlaps
static void *gr_comm_buf;         // buffer for MPI transfers

// First several functions are defined only in parallel mode
//===========================================

static void RecoverCommandLine(int *argc_p,char ***argv_p)
/* eliminate all NULL pointers from argv, shift the rest, and adjust argc accordingly.
 * Used in InitComm
 * */
{
	int i,j;

	for (i=0,j=0;i<(*argc_p);i++) {
		if ((*argv_p)[i]==NULL) j++;
		else if (j!=0) (*argv_p)[i-j]=(*argv_p)[i];
	}
	(*argc_p)-=j;
}
//============================================================

INLINE size_t IndexBlock(const size_t x,const size_t y,const size_t z,const size_t lengthY)
// index block; used in BlockTranspose
{
	return((z*lengthY+y)*gridX+x);
}

//============================================================

INLINE int CalcPartner(const int tran)
/* calculate ringid of partner processor for current transmission; used in BlockTranspose. Many
 * different implementations are possible; the only requirements are
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

//============================================================

void CatNFiles(const char *dir,const char *tmpl,const char *dest)
/* cat several temporary files (one for each processor, names defined by the template 'tmpl' that
 * should contain %d to be replaced by ringid). Files are located in directory 'dir'. Combined into
 * 'dest' in the same directory. Afterwards temporary files are removed.
 */
{
	int i,c;
	FILE *in,*out;
	char fname_out[MAX_TMP_FNAME],fname_in[MAX_TMP_FNAME];

	// produce full path of destination file and open it
	sprintf(fname_out,"%s/%s",directory,dest);
	out=FOpenErr(fname_out,"w",ONE_POS);
	for (i=0;i<nprocs;i++) {
		// produce full path of tmp file and open it
		sprintf(fname_in,"%s/",directory);
		sprintf(fname_in+strlen(fname_in),tmpl,i);
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
#endif

//============================================================

void InitComm(int *argc_p,char ***argv_p)
// initialize communications in the beginning of the program
{
#ifdef ADDA_MPI
	int dcmplx_blocklength[2]={1,1},ver,subver;
	MPI_Aint dcmplx_displs[2] = {0,1};
	MPI_Datatype dcmplx_type[2];

	/* MPI_Init may alter argc and argv and interfere with normal parsing of command line
	 * parameters. The way of altering is implementation depending. MPI searches for MPI parameters
	 * in the command line and removes them (we assume some kind of removing does take place -
	 * otherwise ADDA will produce error 'unknown parameter'). The best would be to change argc and
	 * argv so that they look like no special command line arguments are present. However,
	 * MPICH 1.2.5, for example, just replaces corresponding parameters by NULLs. To incorporate it
	 * we introduce special function to restore the command line
	 */
	MPI_Init(argc_p,argv_p);
	tstart_main = GET_TIME(); // initialize program time
	RecoverCommandLine(argc_p,argv_p);
	// initialize ringid and nprocs
	MPI_Comm_rank(MPI_COMM_WORLD,&ringid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	// initialize Ntrans
	if (IS_EVEN(nprocs)) Ntrans=nprocs-1;
	else Ntrans=nprocs;
	// Create MPI-type for sending dcomplex-variables
	dcmplx_type[0] = MPI_DOUBLE; dcmplx_type[1] = MPI_DOUBLE;
	MPI_Type_struct(2,dcmplx_blocklength,dcmplx_displs,dcmplx_type,&mpi_dcomplex);
	MPI_Type_commit(&mpi_dcomplex);
	// check MPI version at runtime
	MPI_Get_version(&ver,&subver);
	if ((ver<MPI_VER_REQ) || ((ver==MPI_VER_REQ) && (subver<MPI_SUBVER_REQ))) LogError(EC_ERROR,
		ONE_POS,"MPI version (%d.%d) is too old. Version %d.%d or newer is required",ver,subver,
		MPI_VER_REQ,MPI_SUBVER_REQ);
	// if MPI crashes, it happens here
	Synchronize();
#elif !defined(PARALLEL)
	nprocs=1;
	ringid=ADDA_ROOT;
#endif
	/* check if weird number of processors is specified; called even in sequential mode to
	 * initialize weird_nprocs
	 */
	CheckNprocs();
}

//============================================================

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
		// wait for all processors
		fflush(stdout);
		Synchronize();
		// finalize MPI communications
		MPI_Finalize();
	}
#endif
	exit(code);
}

//============================================================

void Synchronize(void)
// synchronizes all processes
{
#ifdef ADDA_MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//============================================================

void MyBcast(void *data,const var_type type,size_t n_elem,TIME_TYPE *timing)
/* casts values stored in '*data' from root processor to all other; works for all types; increments
 * 'timing' (if not NULL) by the time used
 */
{
#ifdef ADDA_MPI
	MPI_Datatype mes_type;
	TIME_TYPE tstart;

	if (timing!=NULL) {
#ifdef SYNCHRONIZE_TIMING
		MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
		tstart=GET_TIME();
	}
	if (type==char_type) mes_type=MPI_CHAR;
	if (type==int_type) mes_type=MPI_INT;
	else if (type==double_type) mes_type=MPI_DOUBLE;
	else if (type==cmplx_type) {
		mes_type=MPI_DOUBLE;
		n_elem*=2;
	}
	else LogError(EC_ERROR,ONE_POS,"MyBcast: variable type %u is not supported",type);

	MPI_Bcast(data,n_elem,mes_type,ADDA_ROOT,MPI_COMM_WORLD);
	if (timing!=NULL) (*timing)+=GET_TIME()-tstart;
#endif
}

//============================================================

void BcastOrient(int *i, int *j, int *k)
// cast current orientation angle (in orientation averaging) to all processes from root
{
#ifdef ADDA_MPI
	int buf[3];

	if (ringid==ADDA_ROOT) {
		buf[0]=*i;
		buf[1]=*j;
		buf[2]=*k;
	}
	MPI_Bcast(buf,3,MPI_INT,ADDA_ROOT,MPI_COMM_WORLD);
	if (ringid!=ADDA_ROOT) {
		*i=buf[0];
		*j=buf[1];
		*k=buf[2];
	}
#endif
}

//============================================================

void AccumulateMax(double *data,double *max)
// given a single double on each processor, accumulates their sum and maximum on root processor
{
#ifdef ADDA_MPI
	double buf;
	// potentially can be optimized by combining into one operation
	MPI_Reduce(data,&buf,1,MPI_DOUBLE,MPI_SUM,ADDA_ROOT,MPI_COMM_WORLD);
	MPI_Reduce(data,max,1,MPI_DOUBLE,MPI_MAX,ADDA_ROOT,MPI_COMM_WORLD);
	if (ringid==ADDA_ROOT) *data=buf;
#endif
}

//============================================================

void Accumulate(double *vector,const size_t n,double *buf,TIME_TYPE *timing)
// gather and add double vector on processor root; total time is saved in timing (NOT incremented)
{
#ifdef ADDA_MPI
	TIME_TYPE tstart;

#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	MPI_Reduce(vector,buf,n,MPI_DOUBLE,MPI_SUM,ADDA_ROOT,MPI_COMM_WORLD);
	if (ringid==ADDA_ROOT) memcpy(vector,buf,n*sizeof(double));
	(*timing)=GET_TIME()-tstart;
#endif
}

//============================================================

void MyInnerProduct(void *data,const var_type type,size_t n_elem,TIME_TYPE *timing)
/* gather values stored in *data, add them and return them in *data; works for all types; increments
 * 'timing' (if not NULL) by the time used; not optimized for long data (allocates memory at every
 * call)
 */
{
#ifdef ADDA_MPI
	size_t size;
	MPI_Datatype mes_type;
	void *temp;
	TIME_TYPE tstart;

	if (timing!=NULL) {
#ifdef SYNCHRONIZE_TIMING
		MPI_Barrier(MPI_COMM_WORLD);  // synchronize to get correct timing
#endif
		tstart=GET_TIME();
	}
	if (type==int_type) {
		mes_type=MPI_INT;
		size=n_elem*sizeof(int);
	}
	else if (type==double_type) {
		mes_type=MPI_DOUBLE;
		size=n_elem*sizeof(double);
	}
	else if (type==cmplx_type) {
		mes_type=MPI_DOUBLE;
		n_elem*=2;
		size=n_elem*sizeof(double);
	}
	else LogError(EC_ERROR,ONE_POS,"MyInnerProduct: variable type %u is not supported",type);

	MALLOC_VECTOR(temp,void,size,ALL);
	MPI_Allreduce(data,temp,n_elem,mes_type,MPI_SUM,MPI_COMM_WORLD);
	memcpy(data,temp,size);
	Free_general(temp);
	if (timing!=NULL) (*timing)+=GET_TIME()-tstart;
#endif
}

//============================================================

void BlockTranspose(doublecomplex *X)
/* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ; specializes at Xmatrix;
 *  do 3 components in one message;
 */
{
#ifdef ADDA_MPI
	TIME_TYPE tstart;
	size_t bufsize,msize,posit,step,y,z;
	int transmission,part,Xpos,Xcomp;
	MPI_Status status;

#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	step=2*local_Nx;
	msize=local_Nx*sizeof(doublecomplex);
	bufsize=6*local_Nz*smallY*local_Nx;

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
	Timing_OneIterComm += GET_TIME() - tstart;
#endif
}

//============================================================

void BlockTranspose_Dm(doublecomplex *X,const size_t lengthY,const size_t lengthZ)
// do the data-transposition, i.e. exchange, between fftX and fftY&fftZ; specialized for D matrix
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
	Timing_Dm_Init_comm += GET_TIME() - tstart;
#endif
}

//============================================================

void ParSetup(void)
// initialize common parameters; need to do in the beginning to enable call to MakeParticle
{
#ifdef PARALLEL
	int unitZ,unitX;
#endif
	// calculate size of 3D grid
	gridX=fftFit(2*boxX,nprocs);
	gridY=fftFit(2*boxY,1);
	gridZ=fftFit(2*boxZ,2*nprocs);
	// initialize some variables
	smallY=gridY/2;
	smallZ=gridZ/2;
	/* if this check is passed then all other multiplications of 2 grids are OK,
	 * except for XY values, used in granule generator
	 */
	gridYZ=MultOverflow(gridY,gridZ,ALL_POS,"gridYZ");
#ifdef PARALLEL
	unitZ=smallZ/nprocs; // this should always be an exact division
	local_z0=ringid*unitZ;
	local_z1=(ringid+1)*unitZ;
	if (local_z1 > boxZ) local_z1_coer=boxZ;
	else local_z1_coer=local_z1;
	unitX=gridX/nprocs;
	local_x0=ringid*unitX;
	local_x1=(ringid+1)*unitX;
#else
	local_z0=0;
	local_z1=smallZ;
	local_z1_coer=boxZ;
	local_x0=0;
	local_x1=gridX;
#endif
	if (local_z1_coer<=local_z0) {
		LogError(EC_INFO,ALL_POS,"No real dipoles are assigned");
		local_z1_coer=local_z0;
	}
	local_Nz=local_z1-local_z0;
	local_Nx=local_x1-local_x0;
	local_Ndip=MultOverflow(boxX*(size_t)boxY,local_z1_coer-local_z0,ALL_POS,"local_Ndip");
	/* conversions to (unsigned long) are needed (to remove warnings) because %z printf argument is
	 * not yet supported by all target compiler environments
	 */
	printf("%i :  %i %i %i %lu %lu \n",
		ringid,local_z0,local_z1_coer,local_z1,(unsigned long)local_Ndip,(unsigned long)local_Nx);
}

//============================================================

void AllGather(void *x_from,void *x_to,const var_type type,size_t n_elem)
// Gather distributed arrays; works for all types
{
#ifdef ADDA_MPI
	// TODO: need to be rewritten when n_elem are unequal on each processor
	MPI_Datatype mes_type;

	if (type==char_type) mes_type = MPI_CHAR;
	else if (type==int_type) mes_type = MPI_INT;
	else if (type==double_type) mes_type = MPI_DOUBLE;
	else if (type==cmplx_type) {
		mes_type = MPI_DOUBLE;
		n_elem *= 2;
	}
	else LogError(EC_ERROR,ONE_POS,"AllGather: variable type %u is not supported",type);

	MPI_Allgather(x_from,n_elem,mes_type,x_to,n_elem,mes_type,MPI_COMM_WORLD);
#endif
}


//============================================================

#ifdef PARALLEL
void CalcLocalGranulGrid(const double z0,const double z1,const double gdZ,const int gZ,
	const int id,int *lz0,int *lz1)
/* calculates starting and ending (+1) cell of granule grid (lz0 & lz1) on a processor with
 * ringid=id
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

//============================================================

void SetGranulComm(const double z0,const double z1,const double gdZ,const int gZ,const size_t gXY,
	size_t max_gran,int *lz0,int *lz1,const int sm_gr)
/* sets communication for granule generator; max_gran - maximum number of granules in one set
 * (used to allocate buffer); sm_gr - whether granules are small (simpler)
 */
{
#ifdef PARALLEL
	int i,loc0,loc1,loc1_prev=0;

	MALLOC_VECTOR(gr_comm_buf,void,max_gran*sizeof(char),ALL);
	if (!sm_gr) {
		if (ringid==ADDA_ROOT) {
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

//============================================================

void CollectDomainGranul(unsigned char *dom,const size_t gXY,const int lz0,
	const int locgZ,TIME_TYPE *timing)
/* collects the map of domain for granule generator on the root processor;
 * timing is incremented by the total time used
 */
{
#ifdef PARALLEL
	int i,unit,index;
	size_t j;
	MPI_Status status;
	TIME_TYPE tstart;

#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	unit=gXY*sizeof(char);
	if (ringid==ADDA_ROOT) {
		index=(lz0+gr_comm_size[ADDA_ROOT])*gXY;
		for (i=ADDA_ROOT+1;i<nprocs;i++) {
			if (gr_comm_size[i]!=0) {
				if (gr_comm_overl[i-1]) {
					index-=gXY;
					memcpy(gr_comm_ob,dom+index,unit);
				}
				MPI_Recv(dom+index,unit*gr_comm_size[i],MPI_UNSIGNED_CHAR,i,0,MPI_COMM_WORLD,
					&status);
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
				MPI_Recv(dom+index-gXY*gr_comm_size[i],unit*gr_comm_size[i],MPI_UNSIGNED_CHAR,i,0,
					MPI_COMM_WORLD,&status);
				if (gr_comm_overl[i]) for (j=0;j<gXY;j++) dom[index-gXY+j]|=gr_comm_ob[j];
				index-=gXY*gr_comm_size[i];
			}
		}
	}
	else if (locgZ!=0) {
		MPI_Send(dom,unit*locgZ,MPI_UNSIGNED_CHAR,ADDA_ROOT,0,MPI_COMM_WORLD);
	}
	(*timing)+=GET_TIME()-tstart;
#endif
}

//============================================================

void FreeGranulComm(const int sm_gr)
/* frees all additional memory used for communications of granule generator;
 * simpler if small granules
 */
{
#ifdef PARALLEL
	Free_general(gr_comm_buf);
	if (!sm_gr && ringid==ADDA_ROOT) {
		Free_general(gr_comm_size);
		Free_general(gr_comm_overl);
		Free_general(gr_comm_ob);
	}
#endif
}

//============================================================

void ExchangeFits(char *data,const size_t n,TIME_TYPE *timing)
/* performs a collective AND operation on the (vector) data; timing is incremented by the total
 * time used.
 * TODO: When MPI_C_BOOL data type will become widely supported, this function should be rewritten
 * using bool input data. However, this may be memory-inefficient if sizeof(MPI_BOOL)>1. This is
 * currently discussed as an errata to MPI 2.2.
 */
{
#ifdef ADDA_MPI
	TIME_TYPE tstart;

#ifdef SYNCHRONIZE_TIMING
	MPI_Barrier(MPI_COMM_WORLD); // synchronize to get correct timing
#endif
	tstart=GET_TIME();
	MPI_Allreduce(data,gr_comm_buf,n,MPI_SIGNED_CHAR,MPI_LAND,MPI_COMM_WORLD);
	memcpy(data,gr_comm_buf,n*sizeof(char));
	(*timing)+=GET_TIME()-tstart;
#endif
}
