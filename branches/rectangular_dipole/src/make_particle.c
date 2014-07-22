/* File: make_particle.c
 * $Date::                            $
 * Descr: this module initializes the dipole set, either using predefined shapes or reading from a file;
 *        includes granule generator
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
// project headers
#include "cmplx.h"
#include "comm.h"
#include "debug.h"
#include "io.h"
#include "memory.h"
#include "param.h"
#include "timing.h"
#include "types.h"
#include "vars.h"
// 3rd party headers
#include "mt19937ar.h"
// system headers
#include <float.h> // for DBL_MAX
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // for time and clock (used for random seed)

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const enum sh shape;
extern const double lambda;
extern double sizeX,dpl,a_eq;
extern const int jagged;
extern const char *shape_fname;
extern const char *shapename;
extern const bool volcor,save_geom;
extern const opt_index opt_sh;
#ifndef SPARSE
extern const int sh_Npars;
extern const double sh_pars[];
extern const char *save_geom_fname;
extern const double gr_vf;
extern double gr_d;
extern const int gr_mat;
extern enum shform sg_format;
extern const bool store_grans;
#endif
// defined and initialized in timing.c
extern TIME_TYPE Timing_Particle;
#ifndef SPARSE
extern TIME_TYPE Timing_Granul,Timing_GranulComm;
#endif

// used in interaction.c
double ZsumShift; // distance between the lowest (in Z) dipoles and its image (in units of d)
// used in param.c
bool volcor_used;                // volume correction was actually employed
const char *sh_form_str1,*sh_form_str2; // strings for log file with shape parameters
size_t gr_N;                     // number of granules
double gr_vf_real;               // actual granules volume fraction
size_t mat_count[MAX_NMAT+1];    // number of dipoles in each domain

// LOCAL VARIABLES

static const char geom_format[]="%d %d %d\n";              // format of the geom file
static const char geom_format_ext[]="%d %d %d %d\n";       // extended format of the geom file
/* DDSCAT shape formats; several format are used, since first variable is unpredictable and last two are not actually
 * used (only to produce warnings)
 */
static const char ddscat_format_read1[]="%*s %d %d %d %d %d %d\n";
static const char ddscat_format_read2[]="%*s %d %d %d %d";
#ifndef SPARSE
static const char ddscat_format_write[]="%zu %d %d %d %d %d %d\n";
#endif
// ratio of scatterer volume to enclosing cube; used for dpl correction and initialization by a_eq
static double volume_ratio;
static double Ndip;             // total number of dipoles (in a circumscribing cube)
static double dpl_def;          // default value of dpl
static int minX,minY,minZ;      // minimum values of dipole positions in dipole file
static FILE * restrict dipfile; // handle of dipole file
static enum shform read_format; // format of dipole file, which is read
static double cX,cY,cZ;         // center for DipoleCoord, it is sometimes used in PlaceGranules

#ifndef SPARSE

// shape parameters
static double coat_x,coat_y,coat_z,coat_r2;
static double ad2,egnu,egeps; // for egg
static double chebeps,r0_2; // for Chebyshev
static int chebn; // for Chebyshev
static double hdratio,invsqY,invsqY2,invsqZ,invsqZ2,haspY,haspZ;
static double xcenter,zcenter; // coordinates of natural particle center (in units of Dx)
static double rc_2,ri_2; // squares of circumscribed and inscribed spheres (circles) in units of Dx
static double boundZ,zcenter1,zcenter2,ell_rsq1,ell_rsq2,ell_x1,ell_x2;
static double rbcP,rbcQ,rbcR,rbcS; // for RBC
static double prang; // for prism
// for axisymmetric; all coordinates defined here are relative
static double * restrict contSegRoMin,* restrict contSegRoMax,* restrict contRo,* restrict contZ;
static double contCurRo,contCurZ;
static int contNseg;
struct segment {
	bool single;           // whether segment consists of a single joint
	int first;             // index of the first point in the segment
	int last;              // index of the last point in the segment
	double zmin;           // minimum z-coordinate of the segment points
	double zmax;           // maximum z-coordinate of the segment points
	double romid;          // ro-coordinate of the point in the middle
	struct segment *left;  // pointer to left subsegment
	struct segment *right; // pointer to right subsegment
	double slope;          // only for single; (z[i+1]-z[i])/(ro[i+1]-ro[i])
	double add;            // only for single; ro[i](1-slope);
};
struct segment * restrict contSeg;

/* TO ADD NEW SHAPE
 * Add here all internal variables (aspect ratios, etc.), which you initialize in InitShape() and use in MakeParticle()
 * afterwards. If you need local, intermediate variables, put them into the beginning of the corresponding function. Add
 * descriptive comments, use 'static'.
 */

// temporary arrays before their real counterparts are allocated
static unsigned char * restrict material_tmp;
static unsigned short * restrict position_tmp;

#endif // !SPARSE

#ifndef SPARSE //much of the functionality here is disabled in sparse mode

// EXTERNAL FUNCTIONS

// chebyshev.c
void ChebyshevParams(double eps_in,int n_in,double *dx,double *dz,double *sz,double *vr);

//======================================================================================================================

static void SaveGeometry(void)
// saves dipole configuration to a file
{
	char fname[MAX_FNAME];
	FILE * restrict geom;
	size_t i,j;
	int mat;
	/* TO ADD NEW FORMAT OF SHAPE FILE
	 * Add code to this function to save geometry in new format. It should consist of:
	 * 1) definition of default filename (by supplying an appropriate extension);
	 * 2) writing header to the beginning of the file;
	 * 3) writing a single line for each dipole (this is done in parallel).
	 * Each part is done in corresponding switch-case sequence
	 */

	TIME_TYPE tstart=GET_TIME();
	// create save_geom_fname if not specified, by adding extension to the shapename
	if (save_geom_fname[0]==0) {
		const char *ext;
		// choose extension
		switch (sg_format) {
			case SF_TEXT:
			case SF_TEXT_EXT: ext="geom"; break;
			case SF_DDSCAT6:
			case SF_DDSCAT7: ext="dat"; break;
			default: LogError(ONE_POS,"Unknown format for saved geometry file (%d)",(int)sg_format);
				// no break
		}
		save_geom_fname=dyn_sprintf("%s.%s",shapename,ext);
	}
	// automatically change format if needed
	if (sg_format==SF_TEXT && Nmat>1) sg_format=SF_TEXT_EXT;
	// choose filename
#ifdef PARALLEL
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/"F_GEOM_TMP,directory,ringid);
#else
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/%s",directory,save_geom_fname);
#endif
	geom=FOpenErr(fname,"w",ALL_POS);
	// print head of file
#ifdef PARALLEL
	if (ringid==0) { // this condition can be different from being root
#endif
		switch (sg_format) {
			case SF_TEXT:
			case SF_TEXT_EXT:
				fprintf(geom,"#generated by ADDA v."ADDA_VERSION"\n"
				             "#shape: '%s'\n"
				             "#box size: %dx%dx%d\n",shapename,boxX,boxY,boxZ);
				if (sg_format==SF_TEXT_EXT) fprintf(geom,"Nmat=%d\n",Nmat);
				break;
			case SF_DDSCAT6:
			case SF_DDSCAT7:
				fprintf(geom,"shape: '%s'; box size: %dx%dx%d; generated by ADDA v."ADDA_VERSION"\n"
				             "%zu = NAT\n"
				             "1 0 0 = A_1 vector\n"
				             "0 1 0 = A_2 vector\n"
				             "1 1 1 = lattice spacings (d_x,d_y,d_z)/d\n",shapename,boxX,boxY,boxZ,nvoid_Ndip);
				if (sg_format==SF_DDSCAT7) fprintf(geom,"%g %g %g = coordinates (x0,y0,z0)/d of the zero dipole "
					"(IX=IY=IZ=0)\n",(1-boxX)/2.0,(1-boxY)/2.0,(1-boxZ)/2.0);
				fprintf(geom,"JA  IX  IY  IZ ICOMP(x,y,z)\n");
				break;
		}
#ifdef PARALLEL
	} // end of if
#endif
	// save geometry
	for(i=0;i<local_nvoid_Ndip;i++) {
		j=3*i;
		switch (sg_format) {
			case SF_TEXT:
				fprintf(geom,geom_format,position[j],position[j+1],position[j+2]);
				break;
			case SF_TEXT_EXT:
				fprintf(geom,geom_format_ext,position[j],position[j+1],position[j+2],material[i]+1);
				break;
			case SF_DDSCAT6:
			case SF_DDSCAT7:
				mat=material[i]+1;
				fprintf(geom,ddscat_format_write,i+local_nvoid_d0+1,position[j],position[j+1],position[j+2],mat,mat,
					mat);
				break;
		}
	}
	FCloseErr(geom,fname,ALL_POS);
#ifdef PARALLEL
	// wait for all processes to save their part of geometry
	Synchronize();
	// combine all files into one and clean
	if (IFROOT) CatNFiles(directory,F_GEOM_TMP,save_geom_fname);
#endif
	if (IFROOT) printf("Geometry saved to file\n");
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================
#define ALLOCATE_SEGMENTS(N) (struct segment *)voidVector((N)*sizeof(struct segment),ALL_POS,"contour segment");

void InitContourSegment(struct segment * restrict seg,const bool increasing)
/* recursively initialize a segment of a contour: allocates memory and calculates all elements some elements are
 * calculated during the forward sweep (from long segments to short), others - during the backward sweep.
 * Recursive function calls incurs certain overhead, however here it is not critical.
 */
{
	int i;
	struct segment * restrict s1,* restrict s2;

	/* Remove constant parts in the beginning and end of segment, if present. After this procedure 'first' is guaranteed
	 * to be less than 'last' by definition of the segment
	 */
	while (contRo[seg->first]==contRo[seg->first+1]) (seg->first)++;
	while (contRo[seg->last-1]==contRo[seg->last]) (seg->last)--;
	if (seg->first+1 == seg->last) { // segment with a single fragment
		seg->single=true;
		seg->zmin=MIN(contZ[seg->first],contZ[seg->last]);
		seg->zmax=MAX(contZ[seg->first],contZ[seg->last]);
		seg->slope=(contZ[seg->last]-contZ[seg->first])/(contRo[seg->last]-contRo[seg->first]);
		seg->add=contZ[seg->first]-contRo[seg->first]*seg->slope;
	}
	else { // divide segment into two, and initialize each of them
		seg->single=false;
		i=(seg->first+seg->last)/2;
		seg->romid=contRo[i];
		// construct subsegments
		s1=ALLOCATE_SEGMENTS(1);
		s2=ALLOCATE_SEGMENTS(1);
		s1->first=seg->first;
		s1->last=s2->first=i;
		s2->last=seg->last;
		// initialize subsegments
		InitContourSegment(s1,increasing);
		InitContourSegment(s2,increasing);
		// calculate zmax and zmin
		seg->zmax=MAX(s1->zmax,s2->zmax);
		seg->zmin=MIN(s1->zmin,s2->zmin);
		// assign new segments to left and right based on 'increasing'
		if (increasing) {
			seg->left=s1;
			seg->right=s2;
		}
		else {
			seg->left=s2;
			seg->right=s1;
		}
	}
}

//======================================================================================================================
#define CHUNK_SIZE 128 // how many numbers are allocated at once for adjustable arrays

static void InitContour(const char *fname,double *ratio,double *shSize)
/* Reads a contour from the file, rotates it so that it starts from a local minimum in ro, then divides it into
 * monotonic (over ro) segments. It produces data, which are later used to test each dipole for being inside the
 * contour. Segments are either increasing or non-decreasing.
 */
{
	size_t line; // current line number
	int nr;      // number of contour points read from the file
	int size;    // current size of the allocated memory for contour
	int i,j,scanned;
	double *bufRo,*bufZ; // temporary buffers
	int *index;
	double ro,z,romin,romax,zmin,zmax,mult,zmid;
	FILE* file;
	bool increasing;
	char linebuf[BUF_LINE];

	D("InitContour has started");
	// Read contour from file
	TIME_TYPE tstart=GET_TIME();
	file=FOpenErr(fname,"r",ALL_POS);
	line=SkipComments(file);
	size=CHUNK_SIZE;
	MALLOC_VECTOR(bufRo,double,size,ALL);
	MALLOC_VECTOR(bufZ,double,size,ALL);
	nr=0;
	// this depends on variables been declared double
	romin=zmin=DBL_MAX;
	romax=zmax=-DBL_MAX;
	// reading is performed in lines
	while(FGetsError(file,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL) {
		// scan numbers in a line
		scanned=sscanf(linebuf,"%lf %lf",&ro,&z);
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			if (scanned!=2) // this in most cases indicates wrong format
				LogError(ONE_POS,"Error occurred during scanning of line %zu in contour file %s",line,fname);
			// check for consistency of input
			if (ro<0) LogError(ONE_POS,"Negative ro-coordinate is found on line %zu in contour file %s",line,fname);
			// update extreme values
			if (z>zmax) zmax=z;
			if (z<zmin) zmin=z;
			if (ro>romax) romax=ro;
			if (ro<romin) romin=ro;
			// add allocated memory to buf, if needed
			if (nr >= size) {
				size+=CHUNK_SIZE;
				REALLOC_VECTOR(bufRo,double,size,ALL);
				REALLOC_VECTOR(bufZ,double,size,ALL);
			}
			bufRo[nr]=ro;
			bufZ[nr]=z;
			nr++;
		}
	}
	FCloseErr(file,fname,ALL_POS);
	Timing_FileIO+=GET_TIME()-tstart;
	// Check number of points read
	if (nr<3) LogError(ONE_POS,"Contour from file %s contains less than three points",fname);

	// Determine initial point with local minimum ro[i-1]>=ro[i]<ro[i+1]
	i=0;
	while (i<nr-1 && bufRo[i]>=bufRo[i+1]) i++;
	if (i==0) { // first point is a minimum candidate
		if (bufRo[0]>bufRo[nr-1]) { // if required, search backwards; guaranteed to converge
			i=nr-1;
			while (bufRo[i]>bufRo[i-1]) i--;
		}
	}
	// if the whole contour is non-decreasing, check for constancy
	else if (i==nr-1 && bufRo[nr-1]==bufRo[0])
		LogError(ONE_POS,"Contour from file %s has zero area. Hence the scatterer is void",fname);
	/* Construct working contour so that its first point = last and is a local minimum. It is done by rotating buf and
	 * adding one extra point. Then free the buffer.
	 */
	MALLOC_VECTOR(contRo,double,nr+1,ALL);
	memcpy(contRo,bufRo+i,(nr-i)*sizeof(double));
	memcpy(contRo+nr-i,bufRo,i*sizeof(double));
	contRo[nr]=contRo[0];
	Free_general(bufRo);
	// same for Z vectors
	MALLOC_VECTOR(contZ,double,nr+1,ALL);
	memcpy(contZ,bufZ+i,(nr-i)*sizeof(double));
	memcpy(contZ+nr-i,bufZ,i*sizeof(double));
	contZ[nr]=contZ[0];
	Free_general(bufZ);
	// scale coordinates to be relative to total diameter, and centered (by z) around 0
	mult=1/(2*romax);
	zmid=(zmax+zmin)/2;
	*ratio=(zmax-zmin)*mult;
	*shSize=2*romax;
	ri_2=romin*romin*mult*mult;
	for (i=0;i<=nr;i++) {
		contRo[i]*=mult;
		contZ[i]=(contZ[i]-zmid)*mult;
	}

	// divide the contour into the segments; actually only the index is constructed marking end points of the segments
	MALLOC_VECTOR(index,int,nr+1,ALL); // this is enough, even if all segments are of one joint
	index[0]=0;
	i=j=1;
	increasing=true;
	while (i<nr) {
		while (i<nr && (increasing == (contRo[i]<contRo[i+1]))) i++;
		index[j]=i;
		j++;
		increasing=!increasing;
	}
	contNseg=j-1;
	/* Calculate maximum and minimum ro for segments; We implicitly use that first segment is increasing, second -
	 * decreasing, and so on.
	 */
	MALLOC_VECTOR(contSegRoMin,double,contNseg,ALL);
	MALLOC_VECTOR(contSegRoMax,double,contNseg,ALL);
	for (j=0;j<contNseg;j++) {
		if (IS_EVEN(j)) {
			contSegRoMin[j]=contRo[index[j]];
			contSegRoMax[j]=contRo[index[j+1]];
		}
		else {
			contSegRoMin[j]=contRo[index[j+1]];
			contSegRoMax[j]=contRo[index[j]];
		}
	}

	// Construct a tree of segments
	contSeg=ALLOCATE_SEGMENTS(contNseg);
	for (j=0;j<contNseg;j++) {
		contSeg[j].first=index[j];
		contSeg[j].last=index[j+1];
		InitContourSegment(contSeg+j,IS_EVEN(j));
	}

	Free_general(index);
	Free_general(contRo);
	Free_general(contZ);
	D("InitContour has finished");
	D("Nseg=%d",contNseg);
	D("minroSq="GFORM_DEBUG,ri_2);
}
#undef CHUNK_SIZE
#undef ALLOCATE_SEGMENTS

//======================================================================================================================

bool CheckContourSegment(struct segment * restrict seg)
/* Checks, whether point is under or above the segment, by traversing the tree of segments. It returns true, if
 * intersecting z value is larger than given z, and false otherwise. Point is defined by local variables contCurRo and
 * contCurZ.
 */
{
	while (true) {
		if (contCurZ < seg->zmin) return true;
		else if (contCurZ > seg->zmax) return false;
		else if (seg->single) return (contCurZ < seg->add + contCurRo*seg->slope);
		else seg=(contCurRo<seg->romid ? seg->left : seg->right);
	}
}

//======================================================================================================================

void FreeContourSegment(struct segment * restrict seg)
/* recursively frees memory allocated for contour segments;
 * Recursive function calls incurs certain overhead, however here it is not critical.
 */
{
	if (!(seg->single)) {
		FreeContourSegment(seg->left);
		FreeContourSegment(seg->right);
	}
}

//======================================================================================================================

#define KEY_LENGTH 2            // length of key for initialization of random generator
#define MAX_ZERO_FITS 1E4       // maximum number of zero fits in a row (each - many granules)
#define MAX_FALSE_SKIP 10       // number of false skips in granule placement to complete the set
#define MAX_FALSE_SKIP_SMALL 10 // the same for small granules
#define MAX_GR_SET USHRT_MAX    // maximum size of granule set
#define MIN_CELL_SIZE 4.0       // minimum cell size for small granules
#define CHECK_CELL(a) CheckCell(gr,vgran,tree_index,Di2,occup[a],&fits) // macro for simplicity
#define CHECK_CELL_TEST(a) (CHECK_CELL(a),fits)                         // ... combined with test for 'fits'

static inline int CheckCell(const double * restrict gr,const double * restrict vgran,
	const unsigned short * restrict tree_index,const double Di2,const int start,bool * restrict fits)
// function that checks whether granule intersects anything in the cell
{
	int index,last,index1;
	double t1,t2,t3;

	last=index=start;
	while (index!=MAX_GR_SET && (*fits)) {
		last=index;
		index1=3*index;
		t1=gr[0]-vgran[index1];
		t2=gr[1]-vgran[index1+1];
		t3=gr[2]-vgran[index1+2];
		if ((t1*t1+t2*t2+t3*t3)<Di2) *fits=false;
		index=tree_index[index];
	}
	return last;
}
//======================================================================================================================

static size_t PlaceGranules(void)
/* Randomly places granules inside the specified domain; Mersenne Twister is used for generating random numbers
 *
 * A simplest algorithm is used: to place randomly a sphere, and see if it overlaps with any dipoles (more exactly:
 * centers of dipoles) of not correct domain; if not, accept it and fill all this dipoles with granules' domain.
 * Optimized to perform in two steps: First it places of set of not-intersecting granules and do only a quick check
 * against the "domain pattern" - coarse representation of the domain. On the second step granules of the whole set are
 * thoroughly checked against the whole domain on each processor. When small granules are used, no domain pattern is
 * used - makes it simpler. Intersection of two granules between the sets is checked only through dipoles, which is not
 * exact, however it allows considering arbitrary complex domains, which is described only by a set of occupied dipoles.
 *
 * This algorithm is unsuitable for high volume fractions, it becomes very slow and for some volume fractions may fail
 * at all (Metropolis algorithm should be more suitable, however it is hard to code for arbitrary domains). Moreover,
 * statistical properties of the obtained granules distribution may be not perfect, however it seems good enough for our
 * applications.
 *
 * Currently it is not working with jagged. That should be improved, by rewriting the jagged calculation throughout the
 * program
 */
{
	int i,j,k,zerofit,last;
	size_t n,count,count_gr,false_count,ui;
	size_t nd;                           // number of dipoles occupied by granules
	int index,index1,index2;             // indices for dipole grid
	int dom_index,dom_index1,dom_index2; // indices for auxiliary grid
	int gX,gY,gZ;                        // auxiliary grid dimensions
	size_t gXY,gr_gN;                    // ... and their products
	size_t avail;                        // number of available (free) domain cells
	int gX2,gY2,gZ2,locgZ2;
	int i0,i1,j0,j1,k0,k1;
	bool fits;
	int cur_Ngr,ig,max_Ngr; // number of granules in a current set, index, and maximum set size
	double gdX,gdY,gdZ,gdXh,gdYh,gdZh; // auxiliary grid cell sizes and their halfs (h)
	int locz0,locz1,locgZ,gr_locgN;
	double R,R2,Di,Di2;          // radius and diameter of granule, and their squares
	double x0,x1,y0,y1,z0,z1;    // where to put random number (inner box)
	int id0,id1,jd0,jd1,kd0,kd1; // dipoles limit that fall inside inner box
	int Nfit;        // number of successfully placed granules in a current set
	double overhead; // estimate of the overhead needed to have exactly needed N of granules
	double tmp1,tmp2,t1,t2,t3;
		// maximum shifts for checks of neighboring cells in auxiliary grid; for 'small' it is the shift in index
	int sx,sy,sz;
	unsigned long key[KEY_LENGTH];   // key to initialize random number generator
	unsigned char * restrict dom;    // information about the domain on a granule grid
	unsigned short * restrict occup; // information about the occupied cells
	int sm_gr;                       // whether granules are small (then simpler algorithm is used)
	unsigned short * restrict tree_index; // index for traversing granules inside one cell (small)
	double * restrict vgran;              // coordinates of a set of granules
	bool * restrict vfit;                 // results of granule fitting on the grid (boolean)
	int * restrict ginX,* restrict ginY,* restrict ginZ; // indices to find dipoles inside auxiliary grid
	int indX,indY,indZ;    // indices for doubled auxiliary grid
	int bit;               // bit position in char of 'dom'
	double gr[3];          // coordinates of a single granule
	FILE * restrict file;  // file for saving granule positions
	char fname[MAX_FNAME]; // filename of file
	double minval;         // minimum size of auxiliary grid

	/* redundant initialization to remove warnings; most of this is due to the fact that code for small and large
	 * granules is largely independent (although there are some common parts, which motivates against complete
	 * separation of them into two functions).
	 */
	zerofit=gX2=gY2=gZ2=locgZ2=id0=id1=jd0=jd1=kd0=kd1=indZ=locgZ=gr_locgN=sx=sy=sz=0;
	gdXh=gdYh=gdZh=0;
	ginX=ginY=ginZ=NULL;
	dom=NULL;
	tree_index=occup=NULL;
	file=NULL;

	// prepare granule file for saving if needed
	if (store_grans && IFROOT) {
		SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_GRANS,directory);
		file=FOpenErr(fname,"w",ONE_POS);
		fprintf(file,"#generated by ADDA v."ADDA_VERSION"\n"
		             "#granule diameter = "GFORM"\n",gr_d);
	}
	// set variables; consider jagged
	Di=gr_d/(gridspace*jagged);
	if (Di<1) LogWarning(EC_WARN,ONE_POS,
		"Granule diameter is smaller than dipole size. It is recommended to increase resolution");
	R=Di/2;
	R2=R*R;
	Di2=4*R2;
	// inner box
	D("gr_N=%zu, Di="GFORMDEF,gr_N,Di);
	x0=R-0.5;
	x1=boxX-R-0.5;
	y0=R-0.5;
	y1=boxY-R-0.5;
	z0=R-0.5;
	z1=boxZ-R-0.5;
	minval=MIN(x1-x0,MIN(y1-y0,z1-z0));
	if (minval<=0) LogError(ONE_POS,"Granule size must be smaller than minimum particle dimension");
	/* initialize auxiliary grid; grid size is chosen to be always <= D/sqrt(3). Thus we can be sure that each cell
	 * contain no more than 1 granule.
	 */
	CheckOverflow(MAX(boxX,MAX(boxY,boxZ))*10/Di,ONE_POS_FUNC);
	tmp1=sqrt(3)/Di;
	gX=(int)ceil((x1-x0)*tmp1);
	gY=(int)ceil((y1-y0)*tmp1);
	gZ=(int)ceil((z1-z0)*tmp1);
	// this should occur only as a consequence of float inaccuracy in error check above
	if (gX==0 || gY==0 || gZ==0) LogError(ONE_POS,"Granule size is too close to minimum particle dimension");
	gdX=(x1-x0)/gX;
	gdY=(y1-y0)/gY;
	gdZ=(z1-z0)/gZ;
	tmp1=MAX(2*Di,MIN_CELL_SIZE);
	/* sets the discrimination for small or large granules; it should guarantee that no divisions by zero occur
	 * afterwards in the code for small granules
	 */
	sm_gr=minval>tmp1 && (gdX<2 || gdY<2 || gdZ<2);
	if (sm_gr) {
		if (IFROOT) printf("Using algorithm for small granules\n");
		/* redefine auxiliary grid; now the grid size is chosen to be always >= 2*D. Thus we can be sure that granule
		 * can intersect with granule in either left or right (x=+-1) cell but not both (and analogously with other
		 * coordinates).
		 */
		tmp1=1/tmp1;
		gX=(int)floor((x1-x0)*tmp1);
		gdX=(x1-x0)/gX;
		gY=(int)floor((y1-y0)*tmp1);
		gdY=(y1-y0)/gY;
		gZ=(int)floor((z1-z0)*tmp1);
		gdZ=(z1-z0)/gZ;
	}
	else { // large granules
		if (IFROOT) printf("Using algorithm for large granules\n");
		gX2=2*gX;
		gdXh=gdX/2;
		gY2=2*gY;
		gdYh=gdY/2;
		gZ2=2*gZ;
		gdZh=gdZ/2;
		/* this sets maximum distance of neighboring cells to check; condition gdX<R can only occur if gX<=7, which is
		 * quite rare, so no optimization is performed. sx>3 can only occur if gX<=2 and then it doesn't make sense to
		 * take bigger sx. Absolutely analogous for y, z.
		 */
		if (gdX<R) sx=3;
		else sx=2;
		if (gdY<R) sy=3;
		else sy=2;
		if (gdZ<R) sz=3;
		else sz=2;
	}
	gXY=MultOverflow(gX,gY,ONE_POS_FUNC);
	gr_gN=MultOverflow(gXY,gZ,ONE_POS_FUNC);
	// calculate maximum number of granules in a grid; crude estimate
	tmp2=(ceil((x1-x0)/Di)+1)*(ceil((y1-y0)/Di)+1)*(ceil((z1-z0)/Di)+1);
	max_Ngr=MIN(MAX_GR_SET,tmp2);
	// local z grid + initialize communications
	D("gr_gN=%zu, max_Ngr=%d",gr_gN,max_Ngr);
	SetGranulComm(z0,z1,gdZ,gZ,gXY,max_Ngr,&locz0,&locz1,sm_gr);
	if (!sm_gr) {
		locgZ=locz1-locz0;
		locgZ2=2*locgZ;
		gr_locgN=gXY*locgZ;
	}
	if (IFROOT) {
		// initialize random generator
		key[0]=(unsigned long)time(NULL);
		key[1]=(unsigned long)(clock());
		init_by_array(key,KEY_LENGTH);
		// allocate memory
		MALLOC_VECTOR(occup,ushort,gr_gN,ONE);
		if (sm_gr) MALLOC_VECTOR(tree_index,ushort,max_Ngr,ONE);
		else MALLOC_VECTOR(dom,uchar,gr_gN,ALL);
	}
	else if (!sm_gr && locgZ!=0) MALLOC_VECTOR(dom,uchar,gr_locgN,ALL);
	MALLOC_VECTOR(vgran,double,3*max_Ngr,ALL);
	MALLOC_VECTOR(vfit,bool,max_Ngr,ALL);
	if (!sm_gr && locgZ!=0) {
		// build some more indices
		MALLOC_VECTOR(ginX,int,gX2+1,ALL);
		MALLOC_VECTOR(ginY,int,gY2+1,ALL);
		MALLOC_VECTOR(ginZ,int,locgZ2+1,ALL);
		for (i=0;i<=gX2;i++) ginX[i]=(int)ceil(x0+i*gdXh);
		id0=ginX[0];
		id1=ginX[gX2];
		for (i=0;i<=gY2;i++) ginY[i]=(int)ceil(y0+i*gdYh);
		jd0=ginY[0];
		jd1=ginY[gY2];
		for (i=0;i<=locgZ2;i++) ginZ[i]=(int)ceil(z0+(i+2*locz0)*gdZh);
		kd0=MAX(ginZ[0],local_z0);
		indZ=1;
		if (kd0>=ginZ[1]) indZ++;
		kd1=MIN(ginZ[locgZ2],local_z1_coer);
	}
	n=count=count_gr=false_count=0;
	nd=0;
	// crude estimate of the probability to place a small granule into domain
	if (sm_gr) overhead=((double)Ndip)/mat_count[gr_mat];
	else overhead=1;
	// main cycle
	D("Starting main iteration cycle");
	while (n<gr_N) {
		if (sm_gr) { // small granules
			// just generate granules
			if (IFROOT) {
				cur_Ngr=MIN(ceil((gr_N-n)*overhead),max_Ngr);
				// generate points and quick check
				ig=false_count=0;
				for (ui=0;ui<gr_gN;ui++) occup[ui]=MAX_GR_SET; // used as undefined
				while (ig<cur_Ngr) {
					count++;
					false_count++;
					fits=true;
					// random position in a grid
					gr[0]=genrand(0,gX);
					gr[1]=genrand(0,gY);
					gr[2]=genrand(0,gZ);
					// coordinates in a grid
					t1=floor(gr[0]);
					t2=floor(gr[1]);
					t3=floor(gr[2]);
					indX=(int)t1;
					indY=(int)t2;
					indZ=(int)t3;
					t1=gr[0]-t1; // t_i are distances to the edges
					t2=gr[1]-t2;
					t3=gr[2]-t3;
					// convert to usual coordinates (in dipole grid)
					gr[0]=gr[0]*gdX+x0;
					gr[1]=gr[1]*gdY+y0;
					gr[2]=gr[2]*gdZ+z0;
					index=indZ*gXY+indY*gX+indX;
					// 'last' is used only if fits, so when this test actually reaches last element
					last=CHECK_CELL(index);
					// weird construction (7-level nested 'ifs') but should be fast
					if (fits) {
						t1*=gdX; // transform shifts to usual coordinates; done only when needed
						sx=0;
						if (t1<Di) {
							if (indX!=0) sx=-1;
						}
						else if ((t1=gdX-t1)<Di && indX!=gX-1) sx=1;
						if (sx==0 || CHECK_CELL_TEST(index+sx)) { // test for x-neighbor
							t2*=gdY;
							sy=0;
							if (t2<Di) {
								if (indY!=0) sy=-gX;
							}
							else if ((t2=gdY-t2)<Di && indY!=gY-1) sy=gX;
							if (sy==0 || CHECK_CELL_TEST(index+sy)) { // test for y-neighbor
								t3*=gdZ;
								sz=0;
								if (t3<Di) {
									if (indZ!=0) sz=-(int)gXY;
								}
								else if ((t3=gdZ-t3)<Di && indZ!=gZ-1) sz=gXY;
								if (sz!=0) {
									if (CHECK_CELL_TEST(index+sz)) { // test for z-neighbor
										if (sy!=0) {
											tmp1=Di2-t2*t2-t3*t3;
											// test for yz- and xyz-neighbors
											if (tmp1>0 && CHECK_CELL_TEST(index+sy+sz )&& sx!=0 && t1*t1<tmp1)
												CHECK_CELL(index+sx+sy+sz);
										}
										else if (sx!= 0 && t1*t1+t3*t3<Di2) CHECK_CELL(index+sx+sz);
									}
								}
								// test for xy-neighbor
								else if (sx!=0 && sy!=0 && t1*t1+t2*t2<Di2) CHECK_CELL(index+sx+sy);
							}
						}
					}
					if (fits) {
						memcpy(vgran+3*ig,gr,3*sizeof(double));
						tree_index[ig]=MAX_GR_SET;
						if (last==MAX_GR_SET) occup[index]=(unsigned short)ig;
						else tree_index[last]=(unsigned short)ig;
						ig++;
						false_count=0;
					}
					if (false_count>MAX_FALSE_SKIP_SMALL) break;
				}
				// real number of placed granules for this set
				cur_Ngr=ig;
			}
		}
		else { // large granules
			// generate domain pattern
			if (locgZ!=0) {
				for (i=0;i<gr_locgN;i++) dom[i]=0;
				/* indices 'index' and 'dom_index' are build up gradually for optimization. Finally,
				 * index=(k-local_z0)*boxXY + j*boxX +i.
				 * TODO: ??? final formula for dom_index is unclear, moreover indZ seems to be not always initialized
				 */
				dom_index2=0;
				index2=(kd0-local_z0)*boxXY;
				bit=((indZ&1)^1)<<2;
				for (k=kd0;k<kd1;k++,index2+=boxXY) {
					index1=index2+jd0*boxX;
					dom_index1=dom_index2;
					indY=1;
					bit&=~2;
					for (j=jd0;j<jd1;j++,index1+=boxX) {
						index=index1+id0;
						dom_index=dom_index1;
						indX=1;
						bit&=~1;
						for (i=id0;i<id1;i++,index++) {
							if (material_tmp[index]!=gr_mat) dom[dom_index]|=(unsigned char)(1<<bit);
							if (i+1==ginX[indX]) {
								indX++;
								bit^=1;
								if (indX&1) dom_index++;
							}
						}
						if (j+1==ginY[indY]) {
							indY++;
							bit^=2;
							if (indY&1) dom_index1+=gX;
						}
					}
					if (k+1==ginZ[indZ]) {
						indZ++;
						bit^=4;
						if (indZ&1) dom_index2+=gXY;
					}
				}
			}
			D("Domain pattern generated");
			// send/collect domain pattern
			CollectDomainGranul(dom,gXY,locz0,locgZ,&Timing_GranulComm);
			if (IFROOT) {
				// analyze domain pattern
				avail=0;
				for (ui=0;ui<gr_gN;ui++) if (dom[ui]!=0xFF) avail++;
				cur_Ngr=MIN(avail,(size_t)max_Ngr);
				tmp1=(gr_N-n)*overhead;
				if (cur_Ngr>tmp1) cur_Ngr=(int)ceil(tmp1);
				// generate points and quick check
				ig=false_count=0;
				for (ui=0;ui<gr_gN;ui++) occup[ui]=MAX_GR_SET; // used as undefined
				while (ig<cur_Ngr) {
					count++;
					// random position in a double grid
					gr[0]=genrand(0,gX2);
					gr[1]=genrand(0,gY2);
					gr[2]=genrand(0,gZ2);
					// coordinates in doubled grid
					indX=(int)floor(gr[0]);
					indY=(int)floor(gr[1]);
					indZ=(int)floor(gr[2]);
					bit=1<<((indX&1)+((indY&1)<<1)+((indZ&1)<<2)); // position bit inside one cell
					// coordinates in usual grid
					indX/=2;
					indY/=2;
					indZ/=2;
					index=indZ*gXY+indY*gX+indX;
					// two simple checks
					if (!(dom[index]&bit) && occup[index]==MAX_GR_SET) {
						// convert to usual coordinates (in dipole grid)
						gr[0]=gr[0]*gdXh+x0;
						gr[1]=gr[1]*gdYh+y0;
						gr[2]=gr[2]*gdZh+z0;
						fits=true;
						false_count++;
						if ((i0=indX-sx)<0) i0=0;
						if ((i1=indX+sx+1)>gX) i1=gX;
						if ((j0=indY-sy)<0) j0=0;
						if ((j1=indY+sy+1)>gY) j1=gY;
						if ((k0=indZ-sz)<0) k0=0;
						if ((k1=indZ+sz+1)>gZ) k1=gZ;
						dom_index2=k0*gXY;
						for (k=k0;k<k1;k++,dom_index2+=gXY) {
							dom_index1=dom_index2+j0*gX;
							for (j=j0;j<j1;j++,dom_index1+=gX) {
								dom_index=dom_index1+i0;
								for (i=i0;i<i1;i++,dom_index++) if (occup[dom_index]!=MAX_GR_SET) {
									index1=3*occup[dom_index];
									t1=gr[0]-vgran[index1];
									t2=gr[1]-vgran[index1+1];
									t3=gr[2]-vgran[index1+2];
									if ((t1*t1+t2*t2+t3*t3)<Di2) {
										fits=false;
										break;
									}
								}
								if (!fits) break;
							}
							if (!fits) break;
						}
						if (fits) {
							memcpy(vgran+3*ig,gr,3*sizeof(double));
							occup[index]=(unsigned short)ig;
							ig++;
							false_count=0;
							/* Here it is possible to correct the domain pattern because of the presence of a new
							 * granule. However it probably will be useful only for large volume fractions
							 */
						}
						if (false_count>MAX_FALSE_SKIP) break;
					}
				}
				// real number of placed granules for this set
				cur_Ngr=ig;
			}
		} // end of large granules
		D("Set of possible granules produced");
		// cast to all processors
		MyBcast(&cur_Ngr,int_type,1,&Timing_GranulComm);
		MyBcast(vgran,double_type,3*cur_Ngr,&Timing_GranulComm);
		count_gr+=cur_Ngr;
		// final check if granules belong to the domain
		for (ig=0;ig<cur_Ngr;ig++) {
			memcpy(gr,vgran+3*ig,3*sizeof(double));
			k0=MAX((int)ceil(gr[2]-R),local_z0);
			k1=MIN((int)floor(gr[2]+R),local_z1_coer-1);
			fits=true;
			index2=(k0-local_z0)*boxXY;
			for (k=k0;k<=k1;k++,index2+=boxXY) {
				tmp1=R2-(gr[2]-k)*(gr[2]-k);
				tmp2=sqrt(tmp1);
				j0=(int)ceil(gr[1]-tmp2);
				j1=(int)floor(gr[1]+tmp2);
				index1=index2+j0*boxX;
				for (j=j0;j<=j1;j++,index1+=boxX) {
					tmp2=sqrt(tmp1-(gr[1]-j)*(gr[1]-j));
					i0=(int)ceil(gr[0]-tmp2);
					i1=(int)floor(gr[0]+tmp2);
					index=index1+i0;
					for (i=i0;i<=i1;i++,index++) {
						if (material_tmp[index]!=gr_mat) {
							fits=false;
							break;
						}
					}
					if (!fits) break;
				}
				if (!fits) break;
			}
			vfit[ig]=fits;
		}
		// collect fits
		ExchangeFits(vfit,cur_Ngr,&Timing_GranulComm);
		// fit dipole grid with successive granules
		Nfit=n;
		for (ig=0;ig<cur_Ngr && n<gr_N;ig++) {
			if (vfit[ig]) { // a successful granule
				n++;
				// fill dipoles in the sphere with granule material
				memcpy(gr,vgran+3*ig,3*sizeof(double));
				k0=MAX((int)ceil(gr[2]-R),local_z0);
				k1=MIN((int)floor(gr[2]+R),local_z1_coer-1);
				index2=(k0-local_z0)*boxXY;
				for (k=k0;k<=k1;k++,index2+=boxXY) {
					tmp1=R2-(gr[2]-k)*(gr[2]-k);
					tmp2=sqrt(tmp1);
					j0=(int)ceil(gr[1]-tmp2);
					j1=(int)floor(gr[1]+tmp2);
					index1=index2+j0*boxX;
					for (j=j0;j<=j1;j++,index1+=boxX) {
						tmp2=sqrt(tmp1-(gr[1]-j)*(gr[1]-j));
						i0=(int)ceil(gr[0]-tmp2);
						i1=(int)floor(gr[0]+tmp2);
						index=index1+i0;
						for (i=i0;i<=i1;i++,index++) {
							material_tmp[index]=(unsigned char)(Nmat-1);
							nd++;
						}
					}
				}
			}
		}
		cur_Ngr=ig; // this is non-trivial only if n=gr_N occurred above
		// save correct granule positions to file
		if (store_grans && IFROOT) for (ig=0;ig<cur_Ngr;ig++) if (vfit[ig]) fprintf(file,GFORM3L"\n",
			gridspace*(vgran[3*ig]-cX),gridspace*(vgran[3*ig+1]-cY),gridspace*(vgran[3*ig+2]-cZ));
		Nfit=n-Nfit;
		/* overhead is estimated based on the estimation of mean value - 1*standard deviation for the probability of
		 * fitting one granule. It is estimated from the Bernoulli statistics: k out of n successful hits.
		 * M(p)=(k+1)/(n+2); s^2(p)=(k+1)(n-k+1)/(n+3)(n+2)^2; M(p)-s(p)=[(k+1)/(n+2)]*[1-sqrt((n-k+1)/(k+1)(n+3))];
		 * overhead=1/latter.
		 */
		overhead=(cur_Ngr+2)/((1-sqrt((cur_Ngr-Nfit+1)/(double)((Nfit+1)*(cur_Ngr+3))))*(Nfit+1));
		if (Nfit!=0) zerofit=0;
		else {
			zerofit++;
			// check if taking too long
			if (zerofit>MAX_ZERO_FITS) {
				MyInnerProduct(&nd,sizet_type,1,&Timing_GranulComm);
				LogError(ONE_POS,"The granule generator failed to reach required volume fraction ("GFORMDEF") of "
					"granules. %zu granules were successfully placed up to a volume fraction of "GFORMDEF".",
					gr_vf,n,((double)nd)/mat_count[gr_mat]);
			}
		}
	}
	if (IFROOT) printf("Granule generator: total random placements= %zu (efficiency 1 = "GFORMDEF")\n"
		               "                   possible granules= %zu (efficiency 2 = "GFORMDEF")\n",
		               count,count_gr/(double)count,count_gr,gr_N/(double)count_gr);
	MyInnerProduct(&nd,sizet_type,1,&Timing_GranulComm);
	// free everything
	if (IFROOT) {
		Free_general(occup);
		if (sm_gr) Free_general(tree_index);
		else Free_general(dom);
	}
	else if (!sm_gr && locgZ!=0) Free_general(dom);
	FreeGranulComm(sm_gr);
	Free_general(vgran);
	Free_general(vfit);
	if (!sm_gr && locgZ!=0) {
		Free_general(ginX);
		Free_general(ginY);
		Free_general(ginZ);
	}
	// close granule file if needed and print info
	if (store_grans && IFROOT) {
		FCloseErr(file,fname,ONE_POS);
		printf("Granule coordinates saved to file\n");
	}
	return nd;
}
#undef KEY_LENGTH
#undef MAX_ZERO_FITS
#undef MAX_FALSE_SKIP
#undef MAX_FALSE_SKIP_SMALL
#undef MAX_GR_SET
#undef MIN_CELL_SIZE
#undef CHECK_CELL
#undef CHECK_CELL_TEST

//======================================================================================================================

#endif // !SPARSE

static void InitDipFile(const char * restrict fname,int *bX,int *bY,int *bZ,int *Nm,const char **rft)
/* read dipole file first to determine box sizes and Nmat; input is not checked for very large numbers (integer
 * overflows) to increase speed; this function opens file for reading, the file is closed in ReadDipFile.
 */
{
	int x,y,z,mat,scanned,mustbe;
	size_t line,skiplines,nd;
	double ds_Ndip; // number of dipoles declared in DDSCAT file
	bool anis_warned;
	bool ds_lf; // whether 6-th line matches pattern for DDSCAT 7 format
	int t2,t3; // dumb variables
	float td1,td2,td3; // dumb float variables
	int maxX,maxY,maxZ,maxN;
	char linebuf[BUF_LINE];
	const char *rf_text;

	TIME_TYPE tstart=GET_TIME();
	dipfile=FOpenErr(fname,"r",ALL_POS);

	// detect file format
	mustbe=UNDEF; // used as indicator whether shape was auto-detected
	/* test for DDSCAT format; in not-DDSCAT format, the line scanned below may be a long comment; therefore we first
	 * skip all comments
	 */
	line=SkipComments(dipfile);
	if (line<=1 // common extensive test for both DDSCAT formats
		&& (SkipNLines(dipfile,1-line),line=1, // to always skip first line
			FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL)
		&& sscanf(linebuf,"%lf",&ds_Ndip)==1 // %zu is not universally supported in scanf
		&& FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL
		&& sscanf(linebuf,"%f %f %f",&td1,&td2,&td3)==3
		&& FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL
		&& sscanf(linebuf,"%f %f %f",&td1,&td2,&td3)==3
		&& FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL
		&& sscanf(linebuf,"%f %f %f",&td1,&td2,&td3)==3
		&& FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL
		&& (ds_lf=(sscanf(linebuf,"%f %f %f",&td1,&td2,&td3)==3), // this line is new in DDSCAT7
			FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL) ) {
		/* Since the 6th or 7th line in DDSCAT format is arbitrary, it can match any test. Therefore, it is impossible
		 * to make rigorous automatic detection of version of DDSCAT format. The failure at this step will only be
		 * detected after scanning the whole file and comparing the number of data lines to ds_Ndip. Having said this,
		 * it is less probable for a comment line to match 7-element pattern than 3-element one, so the result of the
		 * following test takes precedence over ds_lf.
		 */
		if (sscanf(linebuf,ddscat_format_read1,&x,&y,&z,&mat,&t2,&t3)==6) {
			read_format=SF_DDSCAT6;
			rf_text="DDSCAT 6 format (FRMFIL)";
			mustbe=6;
			skiplines=line-1; // should be 6
		}
		else if (ds_lf && FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL
			&& sscanf(linebuf,ddscat_format_read1,&x,&y,&z,&mat,&t2,&t3)==6) {
			read_format=SF_DDSCAT7;
			rf_text="DDSCAT 7 format (FROM_FILE)";
			mustbe=6;
			skiplines=line-1; // should be 7
		}
	}
	if (mustbe==UNDEF) { // test for ADDA text formats, restarting the file
		fseek(dipfile,0,SEEK_SET);
		line=SkipComments(dipfile);
		/* scanf and analyze Nmat; if there is blank line between comments and Nmat, it fails later; the value of Nmat
		 * obtained here is not actually relevant, the main factor is maximum domain number among all dipoles.
		 */
		scanned=fscanf(dipfile,"Nmat=%d\n",Nm);
		if (scanned==EOF) LogError(ONE_POS,"No dipole positions are found in %s",fname);
		else if (scanned==0) { // no "Nmat=..."
			/* It is hard to perform rigorous test for SF_TEXT since any number of blank lines can be present before the
			 * data lines. Therefore this format is determined by exclusion of other possibilities. Hence, errors in
			 * files of other formats will most likely result in assignment of SF_TEXT format to it.
			 */
			read_format=SF_TEXT;
			rf_text="ADDA text format (single domain)";
			mustbe=3;
			skiplines=line;
		}
		else { // "Nmat=..." present
			read_format=SF_TEXT_EXT;
			rf_text="ADDA text format (multi-domain)";
			mustbe=4;
			skiplines=line+1;
		}
	}
	/* TO ADD NEW FORMAT OF SHAPE FILE
	 * Add code to detect file format to this if-else sequence in appropriate place. Currently SF_TEXT is very flexible
	 * - so it is not strictly checked against, effectively making it the default format. So test for new format should
	 * probably go before SF_TEXT. Failure of detection of other formats is indicated by mustbe==UNDEF. The test if
	 * successful should at least assign:
	 * 1) read_format to a handle from const.h;
	 * 2) rf_text to some literal descriptive string;
	 * 2) mustbe to a number of values to be scanned in each data line;
	 * 3) skiplines to a number of non-data files in the beginning of the file
	 */
	// currently this test is redundant
	if (mustbe==UNDEF) LogError(ONE_POS,"Failed to recognize the format of shape file %s",fname);

	// scan main part of the file
	fseek(dipfile,0,SEEK_SET);
	line=SkipNLines(dipfile,skiplines);
	maxX=maxY=maxZ=INT_MIN;
	minX=minY=minZ=INT_MAX;
	maxN=1;
	anis_warned=false;
	// reading is performed in lines
	nd=0;
	scanned=0; // redundant initialization to remove warnings
	while(FGetsError(dipfile,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL) {
		// scan numbers in a line
		switch (read_format) {
			case SF_TEXT: scanned=sscanf(linebuf,geom_format,&x,&y,&z); break;
			case SF_TEXT_EXT: scanned=sscanf(linebuf,geom_format_ext,&x,&y,&z,&mat); break;
			case SF_DDSCAT6:
			case SF_DDSCAT7: // for ddscat formats, only first material is used, other two are ignored
				scanned=sscanf(linebuf,ddscat_format_read1,&x,&y,&z,&mat,&t2,&t3);
				if (!anis_warned && (t2!=mat || t3!=mat)) {
					LogWarning(EC_WARN,ONE_POS,"Anisotropic dipoles are detected in file %s (first on line %zu). ADDA "
						"ignores this anisotropy, using only the identifier of x-component of refractive index as "
						"domain number",fname,line);
					anis_warned=true;
				}
				break;
		}
		/* TO ADD NEW FORMAT OF SHAPE FILE
		 * Add code to scan a single data line and perform consistency checks if necessary. The common variables to be
		 * scanned are 'x','y','z' (integer positions of dipoles) and possibly 'mat' - domain number. 'scanned' should
		 * be set to be tested below.
		 */
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			if (scanned!=mustbe) LogError(ONE_POS,"%s was detected, but error occurred during scanning of line %zu "
				"from dipole file %s",rf_text,line,fname); // this in most cases indicates wrong format
			nd++;
			if (read_format!=SF_TEXT) {
				if (mat<=0) LogError(ONE_POS,"%s was detected, but nonpositive material number (%d) encountered during "
					"scanning of line %zu from dipole file %s",rf_text,mat,line,fname);
				else if (mat>maxN) maxN=mat;
			}
			/* TO ADD NEW FORMAT OF SHAPE FILE
			 * The code above processes scanned mat. It is fairly general but may need modification for a new shape
			 * format. In particular the code should be executed if and only if the shape format is (potentially)
			 * multi-domain.
			 */
			// update maxima and minima
			if (x>maxX) maxX=x;
			if (x<minX) minX=x;
			if (y>maxY) maxY=y;
			if (y<minY) minY=y;
			if (z>maxZ) maxZ=z;
			if (z<minZ) minZ=z;
		}
	}

	// consistency checks and assignment of return values
	switch (read_format) {
		case SF_TEXT: break; // no specific tests
		case SF_TEXT_EXT:
			if (*Nm!=maxN) LogWarning(EC_WARN,ONE_POS,"Nmat (%d), as given in %s, is not equal to the maximum domain "
				"number (%d) among all specified dipoles; hence the former is ignored",*Nm,fname,maxN);
			break;
		case SF_DDSCAT6:
		case SF_DDSCAT7:
			if (nd!=ds_Ndip) LogWarning(EC_WARN,ONE_POS,"Number of dipoles (%.0f), as given in the beginning of %s, is "
				"not equal to the number of data lines actually present and scanned (%zu)",ds_Ndip,fname,nd);
			break;
	}
	nvoid_Ndip=nd*jagged*jagged*jagged;
	/* TO ADD NEW FORMAT OF SHAPE FILE
	 * If any consistency checks for the whole file can be performed, add them as a new case above.
	 */
	*rft=rf_text;
	*Nm=maxN;
	// set grid (box) sizes
	*bX=jagged*(maxX-minX+1);
	*bY=jagged*(maxY-minY+1);
	*bZ=jagged*(maxZ-minZ+1);
	/* return to first data line for ReadDipFile; not optimal way, but works more robustly when non-system EOL is used
	 * in data file
	 */
	fseek(dipfile,0,SEEK_SET);
	SkipNLines(dipfile,skiplines);
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

static void ReadDipFile(const char * restrict fname)
/* read dipole file; no consistency checks are made since they are made in InitDipFile. The file is opened in
 * InitDipFile; this function only closes the file.
 *
 * The operation is quite different in FFT and sparse modes. In FFT mode only material is set here, while position is
 * set in the main loop over shapes in MakeParticle(). By contrast, in sparse mode position is set here as well (since
 * the loop in MakeParticle() is skipped altogether.
 */
{
	int x,y,z,x0,y0,z0,mat,scanned;
	size_t index=0,line=0;
	char linebuf[BUF_LINE];	
#ifndef SPARSE
	// to remove possible overflows
	size_t boxX_l=(size_t)boxX;
#endif // !SPARSE

	TIME_TYPE tstart=GET_TIME();
	
	mat=1; // the default value for single-domain shape formats
	scanned=0; // redundant initialization to remove warnings
	while(fgets(linebuf,BUF_LINE,dipfile)!=NULL) {
		// scan numbers in a line
		switch (read_format) {
			case SF_TEXT: scanned=sscanf(linebuf,geom_format,&x0,&y0,&z0); break;
			case SF_TEXT_EXT: scanned=sscanf(linebuf,geom_format_ext,&x0,&y0,&z0,&mat); break;
			case SF_DDSCAT6:
			case SF_DDSCAT7: scanned=sscanf(linebuf,ddscat_format_read2,&x0,&y0,&z0,&mat); break;
		}
		/* TO ADD NEW FORMAT OF SHAPE FILE
		 * Add code to scan a single data line. The common variables to be scanned are 'x0','y0','z0' (integer positions
		 * of dipoles) and  possibly 'mat' - domain number. 'scanned' should be set to be tested against EOF below. The
		 * code is similar to the one in InitDipFile() but can be simpler because no consistency checks are performed.
		 */
		line++;
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			// shift dipole position to be nonnegative
			x0-=minX;
			y0-=minY;
			z0-=minZ;
			// initialize box jagged*jagged*jagged instead of one dipole
#ifndef SPARSE
			for (z=jagged*z0;z<jagged*(z0+1);z++) if (z>=local_z0 && z<local_z1_coer)
				for (x=jagged*x0;x<jagged*(x0+1);x++) for (y=jagged*y0;y<jagged*(y0+1);y++) {
					index=(z-local_z0)*boxXY+y*boxX_l+x;
					if (material_tmp[index]!=Nmat)
						LogError(ONE_POS,"Duplicate dipole was found at line %zu in dipole file %s",line,fname);
					material_tmp[index]=(unsigned char)(mat-1);
			}
#else
			for (z=0;z<jagged;z++) for (y=0;y<jagged;y++) for (x=0;x<jagged;x++) {
				if ((index >= local_nvoid_d0) && (index < local_nvoid_d1)) {
					material[index-local_nvoid_d0]=(unsigned char)(mat-1);
					position_full[3*index]=x0*jagged+x;
					position_full[3*index+1]=y0*jagged+y;
					position_full[3*index+2]=z0*jagged+z;
				}
				index++;
			}
#endif // SPARSE
		}
	}
	FCloseErr(dipfile,fname,ALL_POS);
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

static int FitBox(const int box)
// finds the smallest value for which program would work (should divide 2*jagged); the limit is also checked
{
	int res;

	res=2*jagged*((box+2*jagged-1)/(2*jagged));
	if (res>BOX_MAX) LogError(ONE_POS,"Derived grid size (%d) is too large (>%d)",res,BOX_MAX);
	return res;
}

//======================================================================================================================

static int FitBox_yz(const double size)
/* given the size of the particle in y or z direction (in units of dipoles), finds the grid size, which would satisfy
 * the FitBox function and so that all dipole centers (more precisely, centers of J^3 dipoles) would fall into the
 * particle (and increasing the number further will produce only the void dipoles).
 *
 * !!! It is still possible, however, that the shape will contain void layers of dipoles because the estimate do not
 * take into account the details of the shape e.g. its curvature. For instance, 'adda -grid 6 -ellipsoid 1 1.5' will
 * result in grid 6x6x10, but the layers z=0, z=10 will be void. This is because the dipole centers in this layers
 * always have non-zero x and y coordinates (at least half-dipole in absolute value) and do not fall inside the sphere,
 * also the points {+-4.5,0,0} do fall into it.
 */
{
	return (2*jagged*(int)floor((size+jagged)/(2*jagged)));
}

//======================================================================================================================

void InitShape(void)
/* perform of initialization of symmetries and boxY, boxZ. Estimate the volume of the particle, when not discretized.
 * Check whether enough refractive indices are specified.
 */
{
	int n_boxX,n_boxY,n_boxZ; // new values for dimensions
	double n_sizeX; // new value for size
	double yx_ratio,zx_ratio;
	double tmp1,tmp2;
#ifndef SPARSE
	double tmp3;
#endif
	TIME_TYPE tstart;
	int Nmat_need,i,temp;
	int small_Nmat=UNDEF; // is set to Nmat, when it is smaller than needed (during prognosis)
	bool size_given_cmd;  // if size is given in the command line
	const char *sizename; // type of input size, used in diagnostic messages
	/* TO ADD NEW SHAPE
	 * Add here all intermediate variables, which are used only inside this function. You may as well use 'tmp1'-'tmp3'
	 * variables defined above.
	 */

	tstart=GET_TIME();

	/* New boxes and sizes are defined only by some shapes, hence they are set to UNDEF here and can be tested against
	 * UNDEF afterwards.
	 */
	n_boxX=n_boxY=n_boxZ=UNDEF;
	n_sizeX=UNDEF;
	sh_form_str2=""; // this variable can be left uninitialized by some shapes

	size_given_cmd=(sizeX!=UNDEF || a_eq!=UNDEF);
	if (sizeX!=UNDEF) sizename="size";
	else if (a_eq!=UNDEF) sizename="eq_rad";
	else sizename=""; // redundant initialization to remove warnings
	// calculate default dpl - 10*sqrt(max(|m|)); for anisotropic each component is considered separately
	tmp2=0;
	for (i=0;i<Ncomp*Nmat;i++) {
		tmp1=cAbs2(ref_index[i]);
		if (tmp2<tmp1) tmp2=tmp1;
	}
	dpl_def=10*sqrt(tmp2);
	// initialization of global option index for error messages
	opt=opt_sh;
	// shape initialization
	switch (shape) {
#ifndef SPARSE //shapes other than "read" are disabled in sparse mode
		case SH_AXISYMMETRIC:
			/* Axisymmetric homogeneous shape, defined by its contour in ro-z plane of the cylindrical coordinate
			 * system. Its symmetry axis coincides with the z-axis, and the contour is read from file. Each line defines
			 * ro and z coordinates of a point, the first and the last points are connected automatically. Linear
			 * interpolation is used between the points.
			 */
			if (IFROOT) sh_form_str1=dyn_sprintf(
				"axisymmetric, defined by a contour in ro-z plane from file %s; diameter:",shape_fname);
			InitContour(shape_fname,&zx_ratio,&n_sizeX);
			yx_ratio=1;
			symZ=false; // input contour is assumed asymmetric over ro-axis
			/* TODO: volume_ratio can be determined from the contour. However, it is not trivial, especially when the
			 * contour intersects itself.
			 */
			volume_ratio=UNDEF;
			Nmat_need=1;
			break;
		case SH_BICOATED: { // based on code by Jin You Lu
			double diskratio,coat_ratio;

			diskratio=sh_pars[0];
			coat_ratio=sh_pars[1];
			coat_r2=0.25*coat_ratio*coat_ratio;
			TestNonNegative(diskratio,"center-to-center distance to diameter ratio");
			TestRangeII(coat_ratio,"inner/outer diameter ratio",0,1);
			if (IFROOT) {
				sh_form_str1="bicoated; diameter(d):";
				sh_form_str2=dyn_sprintf(", center-center distance R_cc/d="GFORM", inner diameter d_in/d="GFORM,
					diskratio,coat_ratio);
			}
			coat_r2=0.25*coat_ratio*coat_ratio;
			hdratio=diskratio/2.0;
			if (diskratio>=1) volume_ratio = 2*PI_OVER_SIX;
			else volume_ratio = PI_OVER_SIX*(2-diskratio)*(1+diskratio)*(1+diskratio)/2;
			yx_ratio=1;
			zx_ratio=diskratio+1;
			Nmat_need=2;
			break;
		}
		case SH_BIELLIPSOID: { // based on code by Alexander Moskalensky
			double aspectY,aspectZ,aspectY2,aspectZ2,aspectXs,aspectY2sc,aspectZ2sc,invmaxX;

			aspectY=sh_pars[0];
			TestPositive(aspectY,"aspect ratio y1/x1");
			aspectZ=sh_pars[1];
			TestPositive(aspectZ,"aspect ratio z1/x1");
			aspectXs=sh_pars[2];
			TestPositive(aspectXs,"aspect ratio x2/x1");
			aspectY2=sh_pars[3];
			TestPositive(aspectY2,"aspect ratio y2/x2");
			aspectZ2=sh_pars[4];
			TestPositive(aspectZ2,"aspect ratio z2/x2");
			// set descriptive string and symmetry
			if (IFROOT) {
				sh_form_str1="biellipsoid; size along x-axis:";
				sh_form_str2=dyn_sprintf("; aspect ratios: y1/x1="GFORM", z1/x1="GFORM", x2/x1="GFORM", y2/x2="GFORM", "
					"z2/x2="GFORM,aspectY,aspectZ,aspectXs,aspectY2,aspectZ2);
			}
			if (aspectY!=1 || aspectY2!=1) symR=false;
			symZ=false; // since upper and lower ellipsoids are generally different both in size and RI
			// set inverse squares of aspect ratios
			invsqY=1/(aspectY*aspectY);
			invsqZ=1/(aspectZ*aspectZ);
			invsqY2=1/(aspectY2*aspectY2);
			invsqZ2=1/(aspectZ2*aspectZ2);
			// determine scale to be used for variables below; and "radii" of ellipsoids
			invmaxX=1/MAX(aspectXs,1);
			ell_x1=0.5*invmaxX;
			ell_x2=ell_x1*aspectXs;
			ell_rsq1=ell_x1*ell_x1;
			ell_rsq2=ell_x2*ell_x2;
			// rescale aspect ratios for second ellipsoid
			aspectY2sc=aspectY2*aspectXs;
			aspectZ2sc=aspectZ2*aspectXs;
			// set z positions of centers and intersection
			zcenter1=-aspectZ2sc*invmaxX/2;
			zcenter2=aspectZ*invmaxX/2;
			boundZ=zcenter1+zcenter2;
			// set box and volume ratios
			volume_ratio=PI_OVER_SIX*(aspectY*aspectZ+aspectY2sc*aspectZ2sc*aspectXs)*invmaxX*invmaxX*invmaxX;
			yx_ratio=MAX(aspectY,aspectY2sc)*invmaxX;
			zx_ratio=(aspectZ+aspectZ2sc)*invmaxX;
			Nmat_need=2;
			break;
		}
		case SH_BISPHERE: { // based on code by Jin You Lu
			double diskratio;

			diskratio=sh_pars[0];
			TestNonNegative(diskratio,"center-to-center distance to diameter ratio");
			if (IFROOT) {
				sh_form_str1="bisphere; diameter(d):";
				sh_form_str2=dyn_sprintf(", center-center distance R_cc/d="GFORM,diskratio);
			}
			hdratio=diskratio/2.0;
			if (diskratio>=1) volume_ratio = 2*PI_OVER_SIX;
			else volume_ratio = PI_OVER_SIX*(2-diskratio)*(1+diskratio)*(1+diskratio)/2;
			yx_ratio=1;
			zx_ratio=diskratio+1;
			Nmat_need=1;
			break;
		}
		case SH_BOX: {
			double aspectY,aspectZ;

			if (sh_Npars==0) {
				if (IFROOT) sh_form_str1="cube; size of edge along x-axis:";
				aspectY=aspectZ=1;
			}
			else { // 2 parameters are given
				aspectY=sh_pars[0];
				TestPositive(aspectY,"aspect ratio y/x");
				aspectZ=sh_pars[1];
				TestPositive(aspectZ,"aspect ratio z/x");
				if (IFROOT) {
					sh_form_str1="rectangular parallelepiped; size along x-axis:";
					sh_form_str2=dyn_sprintf(", aspect ratio y/x="GFORM", z/x="GFORM,aspectY,aspectZ);
				}
			}
			if (aspectY!=1) symR=false;
			// set half-aspect ratios
			haspY=aspectY/2;
			haspZ=aspectZ/2;
			volume_ratio=aspectY*aspectZ;
			yx_ratio=aspectY;
			zx_ratio=aspectZ;
			Nmat_need=1;
			break;
		}
		case SH_CAPSULE: {
			double diskratio;

			diskratio=sh_pars[0];
			TestNonNegative(diskratio,"height to diameter ratio");
			if (IFROOT) {
				sh_form_str1="capsule; diameter(d):";
				sh_form_str2=dyn_sprintf(", cylinder height h/d="GFORM,diskratio);
			}
			hdratio=diskratio/2;
			volume_ratio = PI_OVER_FOUR*diskratio + PI_OVER_SIX;
			yx_ratio=1;
			zx_ratio=diskratio+1;
			Nmat_need=1;
			break;
		}
		case SH_CHEBYSHEV: {
			double Dx,Dz; // grid sizes along x and z (in units of scale parameter r0)
			double sz;    // coordinate of center of Dz relative to sphere origin (in units of r0)
			double ae;    // |eps|

			chebeps=sh_pars[0];
			TestRangeII(chebeps,"height to diameter ratio",-1,1);
			chebn=sh_pars[1];
			ConvertToInteger(sh_pars[1],"number of sides",&chebn);
			TestPositive_i(chebn,"order n");
			ChebyshevParams(chebeps,chebn,&Dx,&Dz,&sz,&volume_ratio);
			if (IFROOT) {
				sh_form_str1="axisymmetric chebyshev particle; size along x-axis (Dx):";
				sh_form_str2=dyn_sprintf(", amplitude eps="GFORM", order n=%d, initial radius r0/Dx="GFORM,
					chebeps,chebn,1/Dx);
			}
			yx_ratio=1;
			zx_ratio=Dz/Dx;
			symZ=(sz!=0);
			zcenter=-sz/Dx;
			r0_2=1/(Dx*Dx);
			ae=fabs(chebeps);
			rc_2=(1+ae)*(1+ae)*r0_2;
			ri_2=(1-ae)*(1-ae)*r0_2;
			Nmat_need=1;
			break;
		}
		case SH_COATED: {
			double coat_ratio;
			char *buf=NULL; // redundant initialization to remove warnings

			coat_ratio=sh_pars[0];
			TestRangeII(coat_ratio,"inner/outer diameter ratio",0,1);
			if (IFROOT) {
				sh_form_str1="coated sphere; diameter(d):";
				buf=dyn_sprintf(", inner diameter d_in/d="GFORM,coat_ratio);
			}
			if (sh_Npars==4) {
				coat_x=sh_pars[1];
				coat_y=sh_pars[2];
				coat_z=sh_pars[3];
				if (coat_x*coat_x+coat_y*coat_y+coat_z*coat_z>0.25*(1-coat_ratio)*(1-coat_ratio))
					PrintErrorHelp("Inner sphere is not fully inside the outer");
				if (IFROOT) buf=rea_sprintf(buf,"\n       position of inner sphere center r/d= "GFORM3V,
					coat_x,coat_y,coat_z);
			}
			else coat_x=coat_y=coat_z=0; // initialize default values
			if (IFROOT) sh_form_str2=buf;
			coat_r2=0.25*coat_ratio*coat_ratio;
			volume_ratio=PI_OVER_SIX;
			if (coat_x!=0) symX=symR=false;
			if (coat_y!=0) symY=symR=false;
			if (coat_z!=0) symZ=false;
			yx_ratio=zx_ratio=1;
			Nmat_need=2;
			break;
		}
		case SH_CYLINDER: {
			double diskratio;

			diskratio=sh_pars[0];
			TestPositive(diskratio,"height to diameter ratio");
			if (IFROOT) {
				sh_form_str1="cylinder; diameter(d):";
				sh_form_str2=dyn_sprintf(", height h/d="GFORM,diskratio);
			}
			hdratio=diskratio/2;
			volume_ratio=PI_OVER_FOUR*diskratio;
			yx_ratio=1;
			zx_ratio=diskratio;
			Nmat_need=1;
			break;
		}
		case SH_EGG: {
			/* determined by equation: (a/r)^2=1+nu*cos(theta)-(1-eps)cos^2(theta) or equivalently:
			 * a^2=r^2+nu*r*z-(1-eps)z^2. Parameters must be 0<eps<=1, 0<=nu<eps. This shape is proposed in:
			 * Hahn D.V., Limsui D., Joseph R.I., Baldwin K.C., Boggs N.T., Carr A.K., Carter C.C., Han T.S., and
			 * Thomas M.E. "Shape characteristics of biological spores", paper 6954-31 to be presented at
			 * "SPIE Defence + Security", March 2008
			 */
			double ad;
			double ct,ct2; // cos(theta0) and its square

			egeps=sh_pars[0];
			TestRangeNI(egeps,"egg parameter epsilon",0,1);
			egnu=sh_pars[1];
			TestRangeIN(egnu,"egg parameter nu",0,egeps);
			// egg shape is symmetric about z-axis (xz and yz planes, but generally NOT xy plane)
			if (egnu!=0) symZ=false;
			/* cos(theta0): ct=-nu/[eps+sqrt(eps^2-nu^2)]; this expression for root of the quadratic equation is used
			 * for numerical stability (i.e. when nu=0); at this theta0 the diameter (maximum width perpendicular to
			 * z-axis) d=Dx is located
			 */
			ct=-egnu/(egeps+sqrt(egeps*egeps-egnu*egnu));
			ct2=ct*ct;
			// Determine ad=(a/d) and its square
			ad2=(1+egnu*ct-(1-egeps)*ct2)/(4*(1-ct2));
			ad=sqrt(ad2);
			tmp1=1/sqrt(egeps+egnu);
			tmp2=1/sqrt(egeps-egnu);
			tmp3=2*(1-egeps);
			/* Center of the computational box (z coordinate): z0=(a/d)*[1/sqrt(eps+nu)+1/sqrt(eps-nu)]/2; but more
			 * numerically stable expression is used (for nu->0). Although it may overflow faster for nu->eps,
			 * volume_ratio (below) will overflow even faster. It is used to shift coordinates from the computational
			 * reference frame (centered at z0) to the natural one
			 */
			zcenter=ad*egnu*(tmp1*tmp1*tmp2*tmp2)/(tmp1+tmp2);
			// (V/d^3)=(4*pi/3)*(a/d)^3*{[2(1-eps)-nu]/sqrt(eps+nu)+[2(1-eps)+nu]/sqrt(eps-nu)}/[nu^2+4(1-eps)]
			volume_ratio=FOUR_PI_OVER_THREE*ad2*ad*((tmp3-egnu)*tmp1+(tmp3+egnu)*tmp2)/(egnu*egnu+2*tmp3);
			if (IFROOT) {
				sh_form_str1="egg; diameter(d):";
				sh_form_str2=dyn_sprintf(", epsilon="GFORM", nu="GFORM", a/d="GFORM,egeps,egnu,ad);
			}
			Nmat_need=1;
			yx_ratio=1;
			zx_ratio=ad*(tmp1+tmp2); // (a/d)*[1/sqrt(eps+nu)+1/sqrt(eps-nu)]
			break;
		}
		case SH_ELLIPSOID: {
			double aspectY,aspectZ;

			aspectY=sh_pars[0];
			TestPositive(aspectY,"aspect ratio y/x");
			aspectZ=sh_pars[1];
			TestPositive(aspectZ,"aspect ratio z/x");
			if (IFROOT) {
				sh_form_str1="ellipsoid; size along x-axis:";
				sh_form_str2=dyn_sprintf(", aspect ratios y/x="GFORM", z/x="GFORM,aspectY,aspectZ);
			}
			if (aspectY!=1) symR=false;
			// set inverse squares of aspect ratios
			invsqY=1/(aspectY*aspectY);
			invsqZ=1/(aspectZ*aspectZ);
			volume_ratio=PI_OVER_SIX*aspectY*aspectZ;
			yx_ratio=aspectY;
			zx_ratio=aspectZ;
			Nmat_need=1;
			break;
		}
		case SH_LINE:
			if (IFROOT) sh_form_str1="line; length:";
			symY=symZ=symR=false;
			n_boxY=n_boxZ=jagged;
			yx_ratio=zx_ratio=UNDEF;
			volume_ratio=UNDEF;
			Nmat_need=1;
			break;
		case SH_PLATE: {
			double diskratio; // ratio of height to diameter

			diskratio=sh_pars[0];
			TestRangeNI(diskratio, "height to diameter ratio",0,1);
			if (IFROOT) {
				sh_form_str1="plate; full diameter(d):";
				sh_form_str2=dyn_sprintf(", height h/d="GFORM,diskratio);
			}
			volume_ratio=PI*diskratio*(6+3*(PI-4)*diskratio+(10-3*PI)*diskratio*diskratio)/24;
			yx_ratio=1;
			zx_ratio=diskratio;
			Nmat_need=1;
			hdratio=zx_ratio/2;
			ri_2=(0.5-hdratio)*(0.5-hdratio);
			break;
		}
		case SH_PRISM: {
			double Dx,Dy; // grid sizes along x and y (in units of a - side length)
			double Ri,Rc; // radii of inscribed and circumscribed circles (in units of a)
			double sx;    // distance between center of circles and center of Dx (in units of a)
			double diskratio; // ratio of height to Dx
			int Nsides; // number of sides

			ConvertToInteger(sh_pars[0],"number of sides",&Nsides);
			TestGreaterThan_i(Nsides,"number of sides",2);
			diskratio=sh_pars[1];
			TestPositive(diskratio,"height to x-size ratio");
			prang=PI/Nsides;
			Rc=0.5/sin(prang);
			Ri=Rc*cos(prang);
			// define grid sizes along x and y, shift of the center along x, and symmetries. Here we assume a=1.
			if (Nsides&1) { // N=2k+1
				Dx=Rc+Ri;
				symX=symR=false;
				sx=(Rc-Ri)/2;
				Dy=2*Rc*sin((Nsides-1)*prang/2);
			}
			else { // N=2k
				Dx=2*Ri;
				sx=0;
				if (Nsides&2) { // N=4k+2
					Dy=2*Rc;
					symR=false;
				}
				else Dy=Dx; // N=4k
			}
			if (IFROOT) {
				sh_form_str1=dyn_sprintf("%d-sided regular prism; size along x-axis (Dx):",Nsides);
				sh_form_str2=dyn_sprintf(", height h/Dx="GFORM", base side a/Dx="GFORM,diskratio,1/Dx);
			}
			xcenter=sx/Dx;
			yx_ratio=Dy/Dx;
			zx_ratio=diskratio;
			hdratio=zx_ratio/2;
			rc_2=Rc*Rc/(Dx*Dx);
			ri_2=Ri*Ri/(Dx*Dx);
			volume_ratio=hdratio*Nsides*Ri/(Dx*Dx);
			Nmat_need=1;
			break;
		}
		case SH_RBC: {
			/* three-parameter shape; developed by K.A.Semyanov,P.A.Tarasov,P.A.Avrorov based on work
			 * P.W.Kuchel and E.D.Fackerell, "Parametric-equation representation of biconcave erythrocytes,"
			 * Bulletin of Mathematical Biology 61, 209-220 (1999).
			 * ro^4+2S*ro^2*z^2+z^4+P*ro^2+Q*z^2+R=0,
			 * ro^2=x^2+y^2, P,Q,R,S are determined by d,h,b,c given in the command line.
			 */
			double h_d,b_d,c_d,h2,b2,c2;

			h_d=sh_pars[0];
			TestPositive(h_d,"ratio of maximum width to diameter");
			b_d=sh_pars[1];
			TestNonNegative(b_d,"ratio of minimum width to diameter");
			if (h_d<=b_d) PrintErrorHelp("given RBC is not biconcave; maximum width is in the center");
			c_d=sh_pars[2];
			TestRangeII(c_d,"relative diameter of maximum width",0,1);
			if (IFROOT) {
				sh_form_str1="red blood cell; diameter(d):";
				sh_form_str2=dyn_sprintf(", maximum and minimum width h/d="GFORM", b/d="GFORM", "
					"diameter of maximum width c/d="GFORM,h_d,b_d,c_d);
			}
			// calculate shape parameters
			h2=h_d*h_d;
			b2=b_d*b_d;
			c2=c_d*c_d;
			/* P={(b/d)^2*[c^4/(h^2-b^2)-h^2]-d^2}/4; Q=(d/b)^2*(P+d^2/4)-b^2/4; R=-d^2*(P+d^2/4)/4; S=-(2P+c^2)/h^2;
			 * here P,Q,R,S are made dimensionless dividing by respective powers of d. Calculation is performed so that
			 * Q is well defined even for b=0. Parameter names are prefixed by 'rbc'.
			 */
			tmp1=((c2*c2/(h2-b2))-h2)/4;
			rbcP=b2*tmp1-0.25;
			rbcQ=tmp1-(b2/4);
			rbcR=-b2*tmp1/4;
			rbcS=-(2*rbcP+c2)/h2;
			yx_ratio=1;
			zx_ratio=h_d;
			volume_ratio=UNDEF;
			Nmat_need=1;
			break;
		}
		case SH_SPHERE:
			if (IFROOT) sh_form_str1="sphere; diameter:";
			volume_ratio=PI_OVER_SIX;
			yx_ratio=zx_ratio=1;
			Nmat_need=1;
			break;
		case SH_SPHEREBOX: {
			double coat_ratio;

			coat_ratio=sh_pars[0];
			TestRangeII(coat_ratio,"sphere diameter/cube edge ratio",0,1);
			if (IFROOT) {
				sh_form_str1="sphere in cube; size of cube edge(a):";
				sh_form_str2=dyn_sprintf(", diameter of sphere d/a="GFORM,coat_ratio);
			}
			coat_r2=0.25*coat_ratio*coat_ratio;
			yx_ratio=zx_ratio=1;
			volume_ratio=1;
			Nmat_need=2;
			break;
		}
#endif // !SPARSE
		case SH_READ: { // it is first, because always relevant
			symX=symY=symZ=symR=false; // input file is assumed fully asymmetric
			const char *rf_text;
			InitDipFile(shape_fname,&n_boxX,&n_boxY,&n_boxZ,&Nmat_need,&rf_text);
			if (IFROOT) sh_form_str1=dyn_sprintf("specified by file %s - %s; size along x-axis:",shape_fname,rf_text);
			yx_ratio=zx_ratio=UNDEF;
			volume_ratio=UNDEF;
			break;
		}
		default: LogError(ONE_POS,"Unknown shape"); // this is mainly to remove 'uninitialized' warnings
			// no break
	}
	/* TO ADD NEW SHAPE
	 * add a case above in alphabetical order (before #endif). Identifier ('SH_...') should be defined inside 'enum sh'
	 * in const.h. This case should
	 * 1) save all the input parameters from array 'sh_pars' to local variables (defined in the beginning of this source
	 *    files)
	 * 2) test all input parameters (for that you're encouraged to use functions from param.h since they would
	 *    automatically produce informative output in case of error).
	 * 3) if shape breaks any symmetry, corresponding variable should be set to false. Do not set any of them to true,
	 *    as they can be set to false by some other factors.
	 *    symX, symY, symZ - symmetries of reflection over planes YZ, XZ, XY respectively.
	 *    symR - symmetry of rotation for 90 degrees over the Z axis
	 * 4) initialize the following:
	 * sh_form_str - descriptive string, should contain format descriptor (like "%g", but using macro GFORM is
	 *               recommended) - it would be replaced by box size along the x-axis afterwards (in param.c).
	 * Either yx_ratio (preferably) or n_boxY. The former is a ratio of particle sizes along y and x axes. Initialize
	 *     n_boxY directly only if it is not proportional to boxX, like in shape LINE above, since boxX is not
	 *     initialized at this moment. If yx_ratio is not initialized, set it explicitly to UNDEF.
	 * Analogously either zx_ratio (preferably) or n_boxZ.
	 * Nmat_need - number of different domains in this shape (void is not included)
	 * volume_ratio - ratio of particle volume to (boxX)^3. Initialize it if it can be calculated analytically or set to
	 *                UNDEF otherwise. This parameter is crucial if one wants to initialize computational grid from
	 *                '-eq_rad' and '-dpl'.
	 * n_boxX - grid size for the particle, defined by shape; initialize only when relevant, e.g. for shapes such as
	 *          'read'.
	 * n_sizeX - absolute size of the particle, defined by shape; initialize only when relevant, e.g. for shapes such as
	 *           'axisymmetric'.
	 *
	 * All other auxiliary variables, which are used in shape generation (MakeParticle(), see below), should be defined
	 * in the beginning of this file. If you need temporary local variables (which are used only in this part of the
	 * code), either use 'tmp1'-'tmp3' or define your own (with more informative names) in the beginning of this
	 * function.
	 */

	// check for redundancy of input data
	if (dpl!=UNDEF) {
		if (size_given_cmd) {
			if (boxX!=UNDEF) PrintError("Too much information is given by setting '-dpl', '-grid', and '-%s'",sizename);
			else if (n_boxX!=UNDEF) PrintError("Too much information is given by setting both '-dpl' and '-%s', while "
				"shape '%s' sets the size of the grid",sizename,shapename);
		}
		else if (n_sizeX!=UNDEF) {
			if (boxX!=UNDEF) PrintError("Too much information is given by setting '-dpl' and '-grid', while shape '%s' "
				"sets the particle size",shapename);
			// currently this can't happen, but may become relevant in the future
			else if (n_boxX!=UNDEF) PrintError("Too much information is given by setting '-dpl', while shape '%s' sets "
				"both the particle size and the size of the grid",shapename);
		}
	}
#ifndef SPARSE //granulation disabled in sparse mode
	// initialize domain granulation
	if (sh_granul) {
		symX=symY=symZ=symR=false; // no symmetry with granules
		if (gr_mat+1>Nmat_need) PrintError("Specified domain number to be granulated (%d) is larger than total number "
			"of domains (%d) for the given shape (%s)",gr_mat+1,Nmat_need,shapename);
		else Nmat_need++;
		// update shapename; use new storage
		shapename=dyn_sprintf("%s_gran",shapename);
	}
#endif
	// check if enough refractive indices or extra
	if (Nmat<Nmat_need) {
		if (prognosis) small_Nmat=Nmat;
		else PrintError("Only %d refractive indices are given. %d are required",Nmat,Nmat_need);
	}
	else if (Nmat>Nmat_need)
		LogWarning(EC_INFO,ONE_POS,"More refractive indices are given (%d) than actually used (%d)",Nmat,Nmat_need);
	Nmat=Nmat_need;

	// check anisotropic refractive indices for symmetries
	if (anisotropy) for (i=0;i<Nmat;i++) symR=symR && ref_index[3*i]==ref_index[3*i+1];

	// determine which size to use
	if (n_sizeX!=UNDEF) {
		if (size_given_cmd) LogWarning(EC_INFO,ONE_POS,"Particle size specified by command line option '-%s' overrides "
			"the internal specification of the shape '%s'. The particle will be scaled accordingly.",
			sizename,shapename);
		else sizeX=n_sizeX;
	}
	// use analytic connection between sizeX and a_eq if available
	if (a_eq!=UNDEF && volume_ratio!=UNDEF) sizeX=pow(FOUR_PI_OVER_THREE/volume_ratio,ONE_THIRD)*a_eq;
	/* Initialization of boxX;
	 * if boxX is not defined by command line, it is either set by shape itself or
	 *   if sizeX is set, boxX is initialized to default
	 *   else dpl is initialized to default (if undefined) and boxX is calculated from sizeX and dpl
	 * else adjust boxX if needed.
	 */
	if (boxX==UNDEF && n_boxX==UNDEF) {
		if (sizeX==UNDEF) {
			// if a_eq is set, but sizeX was not initialized before - error
			if (a_eq!=UNDEF) PrintError("Grid size can not be automatically determined from equivalent radius and dpl "
				"for shape '%s', because its volume is not known analytically. Either use '-size' instead of '-eq_rad' "
				"or specify grid size manually by '-grid'.",shapename);
			// default value for boxX; FitBox is redundant but safer for future changes
			boxX=FitBox(DEF_GRID);
		}
		else {
			if (dpl==UNDEF) {
				// use default dpl, but make sure that it does not produce too small grid (e.g. for nanoparticles).
				temp=(int)ceil(sizeX*dpl_def/lambda);
				boxX=FitBox(MAX(temp,MIN_AUTO_GRID));
				if (small_Nmat!=UNDEF) PrintError("Given number of refractive indices (%d) is less than number of "
					"domains (%d). Since computational grid is initialized based on the default dpl, it may change "
					"depending on the actual refractive indices.",small_Nmat,Nmat_need);
			}
			else { // if dpl is given in the command line; then believe it
				boxX=FitBox((int)ceil(sizeX*dpl/lambda));
				dpl=UNDEF; // dpl is given correct value in make_particle()
			}
		}
	}
	else {
		// warnings are issued if specified boxX need to be adjusted, especially when '-size' is used
		if (boxX!=UNDEF) temp=boxX;
		else temp=n_boxX; // this happens only if n_boxX!=UNDEF
		if ((boxX=FitBox(temp))!=temp) {
			if (sizeX==UNDEF) LogWarning(EC_WARN,ONE_POS,"boxX has been adjusted from %i to %i. Size along X-axis in "
				"the shape description is the size of new (adjusted) computational grid.",temp,boxX);
			else LogWarning(EC_WARN,ONE_POS,"boxX has been adjusted from %i to %i. Size specified by the command line "
				"option '-size' is used for the new (adjusted) computational grid.",temp,boxX);
		}
		if (n_boxX!=UNDEF && n_boxX>boxX)
			PrintError("Particle (boxX=%d) does not fit into specified boxX=%d",n_boxX,boxX);
	}
	/* If shape is determined by ratios, calculate proposed grid sizes along y and z axes. Either ratios or n_box should
	 * necessarily be defined.
	 */
	if (yx_ratio!=UNDEF) n_boxY=FitBox_yz(yx_ratio*boxX);
	else if (n_boxY==UNDEF) LogError(ONE_POS,"Both yx_ratio and n_boxY are undefined");
	if (zx_ratio!=UNDEF) n_boxZ=FitBox_yz(zx_ratio*boxX);
	else if (n_boxZ==UNDEF) LogError(ONE_POS,"Both zx_ratio and n_boxZ are undefined");
	// set boxY and boxZ
	if (boxY==UNDEF) { // assumed that boxY and boxZ are either both defined or both not defined
		boxY=FitBox(n_boxY);
		boxZ=FitBox(n_boxZ);
	}
	else {
		temp=boxY;
		if ((boxY=FitBox(boxY))!=temp) LogWarning(EC_WARN,ONE_POS,"boxY has been adjusted from %i to %i",temp,boxY);
		temp=boxZ;
		if ((boxZ=FitBox(boxZ))!=temp) LogWarning(EC_WARN,ONE_POS,"boxZ has been adjusted from %i to %i",temp,boxZ);
		// this error is not duplicated in the log file since it does not yet exist
		if (n_boxY>boxY || n_boxZ>boxZ)
			PrintError("Particle (boxY,Z={%d,%d}) does not fit into specified boxY,Z={%d,%d}",n_boxY,n_boxZ,boxY,boxZ);
	}
#ifndef SPARSE //this check is not needed in sparse mode
	// initialize number of dipoles; first check that it fits into size_t type
	double tmp=((double)boxX)*((double)boxY)*((double)boxZ);
	if (tmp > SIZE_MAX) LogError(ONE_POS,"Total number of dipoles in the circumscribing box (%.0f) is larger than "
		"supported by size_t type on this system (%zu). If possible, please recompile ADDA in 64-bit mode.",
		tmp,SIZE_MAX);
#endif // !SPARSE
	Ndip=boxX*boxY*boxZ;
	// initialize maxiter; not very realistic
	if (maxiter==UNDEF) maxiter=MIN(INT_MAX,3*Ndip);
	// some old, not really logical heuristics for Ntheta, but better than constant value
	if (nTheta==UNDEF) {
		if (Ndip<1000) nTheta=91;
		else if (Ndip<10000) nTheta=181;
		else if (Ndip<100000) nTheta=361;
		else nTheta=721;
	}
	Timing_Particle = GET_TIME() - tstart;
}

//======================================================================================================================

void MakeParticle(void)
// creates a particle; initializes all dipoles counts, dpl, gridspace
{
	
	size_t index,dip,i3;
	TIME_TYPE tstart;
	int i;
#ifndef SPARSE
	double jcX,jcY,jcZ; // center for jagged
	size_t local_nRows_tmp;
	int j,k,ns;	
	double tmp1,tmp2,tmp3;
	double xr,yr,zr,xcoat,ycoat,zcoat,r2,ro2,z2,zshift,xshift;
	int local_z0_unif; // should be global or semi-global
	int largerZ,smallerZ; // number of larger and smaller z in intersections with contours
	int xj,yj,zj;
	int mat;
	unsigned short us_tmp;
	TIME_TYPE tgran;
#endif // !SPARSE

	/* TO ADD NEW SHAPE
	 * Add here all intermediate variables, which are used only inside this function. You may as well use 'tmp1'-'tmp3'
	 * variables defined above.
	 */
	tstart=GET_TIME();
	
	cX=(boxX-1)/2.0;
	cY=(boxY-1)/2.0;
	cZ=(boxZ-1)/2.0;
		
#ifndef SPARSE //shapes other than "read" are disabled in sparse mode
	index=0;	
	// assumed that box's are even
	jcX=jcY=jcZ=jagged/2.0;

	local_nRows_tmp=MultOverflow(3,local_Ndip,ALL_POS,"local_nRows_tmp");
	
	/* allocate temporary memory; even if prognosis, since they are needed for exact estimation; they will be
	 * reallocated afterwards (when local_nRows is known).
	 */
	MALLOC_VECTOR(material_tmp,uchar,local_Ndip,ALL);
	MALLOC_VECTOR(position_tmp,ushort,local_nRows_tmp,ALL);

	for(k=local_z0;k<local_z1_coer;k++) for(j=0;j<boxY;j++) for(i=0;i<boxX;i++) {
		xj=jagged*(i/jagged)-boxX/2;
		yj=jagged*(j/jagged)-boxY/2;
		zj=jagged*(k/jagged)-boxZ/2;
		/* all coordinates are scaled by the same box size (boxX), so yr and zr are not necessarily in fixed
		 * ranges (like from -1/2 to 1/2). This is done to treat adequately cases when particle dimensions are
		 * the same (along different axes), but e.g. boxY!=boxX (so there are some extra void dipoles). All
		 * anisotropies in the particle itself are treated in the specific shape modules below (see e.g.
		 * ELLIPSOID).
		 */
		xr=(xj+jcX)/(boxX);
		yr=(yj+jcY)/(boxX);
		zr=(zj+jcZ)/(boxX);

		mat=Nmat; // corresponds to void

		switch (shape) {
			case SH_AXISYMMETRIC:
				ro2=xr*xr+yr*yr;
				if (ro2>=ri_2 && ro2<=0.25) {
					largerZ=smallerZ=0;
					contCurRo=sqrt(ro2);
					contCurZ=zr;
					for (ns=0;ns<contNseg;ns++) if (contCurRo>=contSegRoMin[ns] && contCurRo<=contSegRoMax[ns])
						CheckContourSegment(contSeg+ns) ? largerZ++ : smallerZ++;
					// check for consistency; if the code is perfect, this is not needed
					if (!IS_EVEN(largerZ+smallerZ)) LogError(ALL_POS,"Point (ro,z)=("GFORMDEF","GFORMDEF") "
						"produced weird result when checking whether it lies inside the contour. Larger than z %d "
						"intersections, smaller - %d.",contCurRo,contCurZ,largerZ,smallerZ);
					if (!IS_EVEN(largerZ)) mat=0;
				}
				break;
			case SH_BICOATED:
				ro2=xr*xr+yr*yr;
				if (ro2<=0.25) {
					tmp1=fabs(zr)-hdratio;
					if (tmp1*tmp1+ro2<=0.25) {
						if (tmp1*tmp1+ro2<=coat_r2) mat=1;
						else mat=0;
					}
				}
				break;
			case SH_BIELLIPSOID:
				if (zr<=boundZ) { // lower ellipsoid
					if (fabs(xr)<=ell_x1) {
						zshift=zr-zcenter1;
						if (xr*xr+yr*yr*invsqY+zshift*zshift*invsqZ<=ell_rsq1) mat=0;
					}
				}
				else { // upper ellipsoid
					if (fabs(xr)<=ell_x2) {
						zshift=zr-zcenter2;
						if (xr*xr+yr*yr*invsqY2+zshift*zshift*invsqZ2<=ell_rsq2) mat=1;
					}
				}
				break;
			case SH_BISPHERE:
				ro2=xr*xr+yr*yr;
				if (ro2<=0.25) {
					tmp1=fabs(zr)-hdratio;
					if (tmp1*tmp1+ro2<=0.25) mat=0;
				}
				break;
			case SH_BOX:
				if (fabs(yr)<=haspY && fabs(zr)<=haspZ) mat=0;
				break;
			case SH_CAPSULE:
				ro2=xr*xr+yr*yr;
				if (ro2<=0.25) {
					tmp1=fabs(zr)-hdratio;
					if (tmp1<=0 || tmp1*tmp1+ro2<=0.25) mat=0;
				}
				break;
			case SH_CHEBYSHEV:
				ro2=xr*xr+yr*yr;
				zshift=zr-zcenter;
				r2=ro2+zshift*zshift;
				if (r2<=ri_2) mat=0;
				else if (r2<=rc_2) {
					/* This can be optimized using Chebyshev polynomials, but would probably be efficient only for
					 * relatively small n.
					 */
					tmp1=1+chebeps*cos(chebn*atan2(sqrt(ro2),zshift));
					if (r2 <= r0_2*tmp1*tmp1) mat=0;
				}
				break;
			case SH_COATED:
				if (xr*xr+yr*yr+zr*zr<=0.25) { // first test to skip some dipoles immediately)
					xcoat=xr-coat_x;
					ycoat=yr-coat_y;
					zcoat=zr-coat_z;
					if (xcoat*xcoat+ycoat*ycoat+zcoat*zcoat<=coat_r2) mat=1;
					else mat=0;
				}
				break;
			case SH_CYLINDER:
				if(xr*xr+yr*yr<=0.25 && fabs(zr)<=hdratio) mat=0;
				break;
			case SH_EGG:
				ro2=xr*xr+yr*yr;
				zshift=zr-zcenter;
				z2=zshift*zshift;
				if (ro2+egeps*z2+egnu*zshift*sqrt(ro2+z2)<=ad2) mat=0;
				break;
			case SH_ELLIPSOID:
				if (xr*xr+yr*yr*invsqY+zr*zr*invsqZ<=0.25) mat=0;
				break;
			case SH_LINE:
				if (yj==0 && zj==0) mat=0;
				break;
			case SH_PLATE:
				ro2=xr*xr+yr*yr;
				if(ro2<=0.25 && fabs(zr)<=hdratio) {
					if (ro2<=ri_2) mat=0;
					else {
						tmp1=sqrt(ro2)-0.5+hdratio; // ro-ri
						if (tmp1*tmp1+zr*zr<=hdratio*hdratio) mat=0;
					}
				}
				break;
			case SH_PRISM:
				xshift=xr-xcenter;
				ro2=xshift*xshift+yr*yr;
				if(ro2<=rc_2 && fabs(zr)<=hdratio) {
					if (ro2<=ri_2) mat=0;
					/* this can be optimized considering special cases for small N. For larger N the relevant
					 * fraction of dipoles decrease as N^-2, so this part is less problematic.
					 */
					else {
						tmp1=cos(fmod(fabs(atan2(yr,xshift))+prang,2*prang)-prang);
						if (tmp1*tmp1*ro2<=ri_2) mat=0;
					}
				}
				break;
			case SH_RBC:
				ro2=xr*xr+yr*yr;
				z2=zr*zr;
				if (ro2*ro2+2*rbcS*ro2*z2+z2*z2+rbcP*ro2+rbcQ*z2+rbcR<=0) mat=0;
				break;
			case SH_READ: break; // just to have a complete set of cases; this cases is treated separately below
			case SH_SPHERE:
				if (xr*xr+yr*yr+zr*zr<=0.25) mat=0;
				break;
			case SH_SPHEREBOX:
				if (xr*xr+yr*yr+zr*zr<=coat_r2) mat=1;
				else if (fabs(yr)<=0.5 && fabs(zr)<=0.5) mat=0;
				break;
		}
		/* TO ADD NEW SHAPE
		 * add a case above (in alphabetical order). Identifier ('SH_...') should be defined inside 'enum sh' in
		 * const.h. This option should set 'mat' - index of domain for a point, specified by {xr,yr,zr} - coordinates
		 * divided by grid size along X (xr from -0.5 to 0.5, others - depending on aspect ratios). C array indexing
		 * used: mat=0 - first domain, etc. If point corresponds to void, do not set 'mat'. If you need temporary local
		 * variables (which are used* only in this part of the code), either use 'tmp1'-'tmp3' or define your own (with
		 * more informative names) in the beginning of this function.
		 */
		position_tmp[3*index]=(unsigned short)i;
		position_tmp[3*index+1]=(unsigned short)j;
		position_tmp[3*index+2]=(unsigned short)k;
		// afterwards multiplied by gridspace
		material_tmp[index]=(unsigned char)mat;
		index++;
	} // End box loop
#else // SPARSE
	// local_nvoid_d0 and local_nvoid_d1 are set earlier in ParSetup()
	local_nvoid_Ndip=local_nvoid_d1-local_nvoid_d0;
	local_Ndip=local_nvoid_Ndip;
	MALLOC_VECTOR(material,uchar,local_nvoid_Ndip,ALL);
	// check if 3*nvoid_Ndip can be computed; check is redundant for sequential mode
	size_t nRows=MultOverflow(3,nvoid_Ndip,ONE_POS_FUNC);
	MALLOC_VECTOR(position_full,int,nRows,ALL);
	memory+=3*sizeof(int)*nvoid_Ndip+sizeof(char)*local_nvoid_Ndip;
#endif // SPARSE
	if (shape==SH_READ) ReadDipFile(shape_fname);
	// initialization of mat_count and dipoles counts
	for(i=0;i<=Nmat;i++) mat_count[i]=0;
#ifdef SPARSE
	for(dip=0;dip<local_Ndip;dip++) mat_count[material[dip]]++;
	MyInnerProduct(mat_count,sizet_type,Nmat+1,NULL);
#else
	for(dip=0;dip<local_Ndip;dip++) mat_count[material_tmp[dip]]++;
	local_nvoid_Ndip=local_Ndip-mat_count[Nmat];
	SetupLocalD();
	MyInnerProduct(mat_count,sizet_type,Nmat+1,NULL);
	nvoid_Ndip=Ndip-mat_count[Nmat];
#endif // !SPARSE
	if (nvoid_Ndip==0) LogError(ONE_POS,"All dipoles of the scatterer are void");
	local_nRows=3*local_nvoid_Ndip;
	// initialize dpl and gridspace
	volcor_used=(volcor && (volume_ratio!=UNDEF));
	if (sizeX==UNDEF) {
		if (a_eq!=UNDEF) dpl=lambda*pow(nvoid_Ndip*THREE_OVER_FOUR_PI,ONE_THIRD)/a_eq;
		else if (dpl==UNDEF) dpl=dpl_def; // default value of dpl
		// sizeX is determined to give correct volume
		if (volcor_used) sizeX=lambda*pow(nvoid_Ndip/volume_ratio,ONE_THIRD)/dpl;
		else sizeX=lambda*boxX/dpl;
	}
	else {
		// dpl is determined to give correct volume
		if (volcor_used) dpl=lambda*pow(nvoid_Ndip/volume_ratio,ONE_THIRD)/sizeX;
		else dpl=lambda*boxX/sizeX;
	}
	// Check consistency for FCD
	if ((IntRelation==G_FCD || PolRelation==POL_FCD) && dpl<=2)
		LogError(ONE_POS,"Too small dpl for FCD formulation, should be at least 2");
	// initialize gridspace and dipvol
	gridspace=lambda/dpl;
	dipvol=gridspace*gridspace*gridspace;
	// initialize equivalent size parameter and cross section
	kd = TWO_PI/dpl;
	/* from this moment on a_eq and all derived quantities are based on the real a_eq, which can in several cases be
	 * slightly different from the one given by '-eq_rad' option.
	 */
	a_eq = pow(THREE_OVER_FOUR_PI*nvoid_Ndip,ONE_THIRD)*gridspace;
	ka_eq = WaveNum*a_eq;
	inv_G = 1/(PI*a_eq*a_eq);
	
#ifndef SPARSE
	// granulate one domain, if needed
	if (sh_granul) {
		tgran=GET_TIME();
		Timing_GranulComm=0;
		// calculate number of granules
		if (mat_count[gr_mat]==0) LogError(ONE_POS,"Domain to be granulated does not contain any dipoles");
		tmp1=gridspace/gr_d;
		tmp2=mat_count[gr_mat]*gr_vf*SIX_OVER_PI;
		tmp3=tmp2*tmp1*tmp1*tmp1;
		CheckOverflow(tmp3,ONE_POS,"gr_N");
		gr_N=(size_t)ceil(tmp3);
		// correct granules diameter to get exact volume fraction (if volume correction is used)
		if (volcor_used) {
			tmp1=gridspace*pow(tmp2/gr_N,ONE_THIRD);
			tmp3=100*fabs((tmp1/gr_d)-1);
			if (tmp3>10) LogWarning(EC_WARN,ONE_POS,
				"Granule size was adjusted by %.0f%% (to "GFORMDEF") to satisfy volume correction",tmp3,tmp1);
			gr_d=tmp1;
		}
		// actually place granules
		mat_count[Nmat-1]=PlaceGranules();
		// calculate exact volume fraction
		gr_vf_real=((double)mat_count[Nmat-1])/mat_count[gr_mat];
		mat_count[gr_mat]-=mat_count[Nmat-1];
		Timing_Granul=GET_TIME()-tgran;
	}
	/* allocate main particle arrays, using precise local_nRows even when prognosis is used to enable save_geom
	 * afterwards.
	 */
	MALLOC_VECTOR(material,uchar,local_nvoid_Ndip,ALL);
	MALLOC_VECTOR(position,ushort,local_nRows,ALL);
	memory+=(3*sizeof(short int)+sizeof(char))*local_nvoid_Ndip;
	// copy nontrivial part of arrays
	index=0;
	for (dip=0;dip<local_Ndip;dip++) if (material_tmp[dip]<Nmat) {
		material[index]=material_tmp[dip];
		memcpy(position+3*index,position_tmp+3*dip,3*sizeof(short int));
		index++;
	}

	// free temporary memory
	Free_general(material_tmp);
	Free_general(position_tmp);
	if (shape==SH_AXISYMMETRIC) {
		for (ns=0;ns<contNseg;ns++) FreeContourSegment(contSeg+ns);
		Free_general(contSegRoMin);
		Free_general(contSegRoMax);
	}
#else
	position=position_full + 3*local_nvoid_d0;
#endif // SPARSE
	// initialize DipoleCoord
	MALLOC_VECTOR(DipoleCoord,double,local_nRows,ALL);
	memory+=3*sizeof(double)*local_nvoid_Ndip;
	double minZco=0; // minimum Z coordinates of dipoles
	for (index=0; index<local_nvoid_Ndip; index++) {
		i3=3*index;
		DipoleCoord[i3] = (position[i3]-cX)*gridspace;
		DipoleCoord[i3+1] = (position[i3+1]-cY)*gridspace;
		DipoleCoord[i3+2] = (position[i3+2]-cZ)*gridspace;
		if (minZco>DipoleCoord[i3+2]) minZco=DipoleCoord[i3+2]; // crude way to find the minimum on the way
	}
	/* test that particle is wholly above the substrate; strictly speaking, we test dipole centers to be above the
	 * substrate - hsub+minZco>0, while the geometric boundary of the particle may still intersect with the substrate.
	 * However, the current test is sufficient to ensure that corresponding routines to calculate reflected Green's
	 * tensor do not fail. And accuracy of the DDA itself is anyway questionable when some of the dipoles are very close
	 * to the substrate (whether they cross it or not).
	 */
	if (surface && hsub<=-minZco) LogError(ALL_POS,"The particle must be entirely above the substrate. There exist a "
		"dipole with z="GFORMDEF" (relative to the center), making specified height of the center ("GFORMDEF") too "
		"small",minZco,hsub);
	// save geometry
	if (save_geom)
#ifndef SPARSE 
		SaveGeometry();
#else
		LogError(ONE_POS,"Saving the geometry is not supported in sparse mode");
#endif // SPARSE

#ifndef SPARSE
	/* adjust z-axis of position vector, to speed-up matrix-vector multiplication a little bit; after this point
	 * 'position(z)' is taken relative to the local_z0.
	 */
	if (local_z0!=0) {
		us_tmp=(unsigned short)local_z0;
		for (dip=2;dip<3*local_nvoid_Ndip;dip+=3) position[dip]-=us_tmp;
	}
	local_Nz_unif=position[3*local_nvoid_Ndip-1]+1;
	local_z0_unif=local_z0; // TODO: should be changed afterwards
#endif // !SPARSE

	box_origin_unif[0]=-gridspace*cX;
	box_origin_unif[1]=-gridspace*cY;
#ifndef SPARSE
	box_origin_unif[2]=gridspace*(local_z0_unif-cZ);
	if (surface) ZsumShift=2*((hsub/gridspace)-cZ+local_z0);
#else
	box_origin_unif[2]=-gridspace*cZ;
	if (surface) ZsumShift=2*((hsub/gridspace)-cZ);
#	ifdef PARALLEL
	AllGather(NULL,position_full,int3_type,NULL);
#	endif
#endif // SPARSE
	
	Timing_Particle += GET_TIME() - tstart;
}
