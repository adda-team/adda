/* File: make_particlce.c
 * $Author$
 * $Date::                            $
 * Descr: this module initializes the dipole set, either using predefined shapes or reading from a
 *        file; includes granule generator
 *
 *        Up to 2004 it was maintained by Alfons Hoekstra.
 *        -----------------------------------------------------------
 *        rewritten,
 *        Michel Grimminck 1995
 *        -----------------------------------------------------------
 *        included ellipsoidal particle for work with Victor Babenko
 *        September 2002
 *        -----------------------------------------------------------
 *        included many more new particles:
 *        leucocyte, stick, rotatable oblate spheroid, lymphocyte,
 *        rotatable RBC, etc etc etc
 *        December 2003 - February 2004, by Konstantin Semyanov
 *          (not used now)
 *        -----------------------------------------------------------
 *        Shapes 'capsule' and 'egg' are implemented by Daniel Hahn and Richard Joseph.
 *        -----------------------------------------------------------
 *        Shape 'axisymmetric' is based on the code by Konstantin Gilev
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h> // for time and clock (used for random seed)
#include <limits.h>
#include <stdbool.h>
#include "vars.h"
#include "const.h"
#include "cmplx.h"
#include "types.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"
#include "io.h"
#include "param.h"
#include "timing.h"
#include "mt19937ar.h"

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const int shape,sh_Npars;
extern const double sh_pars[];
extern const int sym_type;
extern const double lambda;
extern double sizeX,dpl,a_eq;
extern const int jagged;
extern const char shape_fname[];
extern char shapename[];
extern char save_geom_fname[];
extern const int volcor,save_geom;
extern opt_index opt_sh;
extern const double gr_vf;
extern double gr_d;
extern const int gr_mat;
extern int sg_format;
extern int store_grans;

// defined and initialized in timing.c
extern TIME_TYPE Timing_Particle,Timing_Granul,Timing_Granul_comm;

// used in param.c
int volcor_used;                 // volume correction was actually employed
char sh_form_str[MAX_PARAGRAPH]; // string for log file with shape parameters
size_t gr_N;                     // number of granules
double gr_vf_real;               // actual granules volume fraction
double mat_count[MAX_NMAT+1];    // number of dipoles in each domain

// LOCAL VARIABLES

static const char geom_format[]="%d %d %d\n";              // format of the geom file
static const char geom_format_ext[]="%d %d %d %d\n";       // extended format of the geom file
/* C99 allows use of %zu for size_t variables, but this is not supported by MinGW due to dependence
 * on Microsoft libraries
 */
static const char ddscat_format[]="%ld %d %d %d %d %d %d\n";// DDSCAT shape format (FRMFIL)
// ratio of scatterer volume to enclosing cube; used for dpl correction and initialization by a_eq
static double volume_ratio;
static double Ndip;            // total number of dipoles (in a circumscribing cube)
static double dpl_def;         // default value of dpl
static int minX,minY,minZ;     // minimum values of dipole positions in dipole file
static FILE *dipfile;          // handle of dipole file
static int read_format;        // format of dipole file, which is read
static char linebuf[BUF_LINE]; // buffer for reading lines from dipole file
double cX,cY,cZ;               // center for DipoleCoord, it is sometimes used in PlaceGranules
// shape parameters
static double coat_ratio,coat_x,coat_y,coat_z,coat_r2;
static double ad2,egnu,egeps,egz0; // for egg
static double hdratio,invsqY,invsqZ,haspY,haspZ;
static double P,Q,R,S; // for RBC
// for axisymmetric; all coordinates defined here are relative
static double *contSegRoMin,*contSegRoMax,*contRo,*contZ;
static double contCurRo, contCurZ, contRoSqMin;
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
struct segment *contSeg;
/* TO ADD NEW SHAPE
 * Add here all internal variables (aspect ratios, etc.), which you initialize in InitShape()
 * and use in MakeParticle() afterwards. If you need local, intermediate variables, put them into
 * the beginning of the corresponding function.
 * Add descriptive comments, use 'static'.
 */

// temporary arrays before their real counterparts are allocated
static unsigned char *material_tmp;
static double *DipoleCoord_tmp;
static unsigned short *position_tmp;

//============================================================

static void SaveGeometry(void)
// saves dipole configuration to .geom file
{
	char fname[MAX_FNAME];
	FILE *geom;
	size_t i,j;
	int mat;

	// create save_geom_fname if not specified
	if (save_geom_fname[0]==0) sprintf(save_geom_fname,"%s.geom",shapename);
	// automatically change format if needed
	if (sg_format==SF_TEXT && Nmat>1) sg_format=SF_TEXT_EXT;
	// choose filename
#ifdef PARALLEL
	sprintf(fname,"%s/" F_GEOM_TMP,directory,ringid);
#else
	sprintf(fname,"%s/%s",directory,save_geom_fname);
#endif
	geom=FOpenErr(fname,"w",ALL_POS);
	// print head of file
#ifdef PARALLEL
	if (ringid==0) { // this condition can be different from being ROOT
#endif
		if (sg_format==SF_TEXT || sg_format==SF_TEXT_EXT) {
			fprintf(geom,"#generated by ADDA v." ADDA_VERSION "\n"
			             "#shape: '%s'\n"
			             "#box size: %dx%dx%d\n",shapename,boxX,boxY,boxZ);
			if (sg_format==SF_TEXT_EXT) fprintf(geom,"Nmat=%d\n",Nmat);
		}
		else if (sg_format==SF_DDSCAT)
			fprintf(geom,"shape: '%s'; box size: %dx%dx%d; generated by ADDA v." ADDA_VERSION "\n"
			             "%0.f = NAT\n"
			             "1 0 0 = A_1 vector\n"
			             "0 1 0 = A_2 vector\n"
			             "1 1 1 = lattice spacings (d_x,d_y,d_z)/d\n"
			             "JA  IX  IY  IZ ICOMP(x,y,z)\n",shapename,boxX,boxY,boxZ,nvoid_Ndip);
#ifdef PARALLEL
	} // end of if
#endif
	// save geometry
	if (sg_format==SF_TEXT) for(i=0;i<local_nvoid_Ndip;i++) {
		j=3*i;
		fprintf(geom,geom_format,position[j],position[j+1],position[j+2]);
	}
	else if (sg_format==SF_TEXT_EXT) for(i=0;i<local_nvoid_Ndip;i++) {
		j=3*i;
		fprintf(geom,geom_format_ext,position[j],position[j+1],position[j+2],material[i]+1);
	}
	else if (sg_format==SF_DDSCAT) for(i=0;i<local_nvoid_Ndip;i++) {
		j=3*i;
		mat=material[i]+1;
		fprintf(geom,ddscat_format,(long)(i+1),position[j],position[j+1],position[j+2],mat,mat,mat);
		/* conversion to long is needed (to remove warnings) because %z printf
		 * argument is not yet supported by all target compiler environments
		 */
	}
	FCloseErr(geom,fname,ALL_POS);
#ifdef PARALLEL
	// wait for all processes to save their part of geometry
	Synchronize();
	// combine all files into one and clean
	if (ringid==ROOT) CatNFiles(directory,F_GEOM_TMP,save_geom_fname);
#endif
	PRINTZ("Geometry saved to file\n");
}

//===========================================================

INLINE void SkipFullLine(FILE* file)
// skips full line in the file, starting from current position; it uses predefined buffer 'linebuf'
{
	do fgets(linebuf,BUF_LINE,file); while (strchr(linebuf,'\n')==NULL && !feof(file));
}

//===========================================================

INLINE char *FgetsError(FILE* file,const char *fname,int *line,const char *s_fname,const int s_line)
/* calls fgets, checks for errors and increments line number; s_fname and s_line are source fname
 * and line number to be shown in error message; result is stored in predefined buffer 'linebuf'.
 */
{
	char *res;

	res=fgets(linebuf,BUF_LINE,file);
	if (res!=NULL) {
		(*line)++;
		if (strchr(linebuf,'\n')==NULL && !feof(file)) LogError(EC_ERROR,ONE,s_fname,s_line,
			"Buffer overflow while scanning lines in file '%s' (size of line %d > %d)",
			fname,*line,BUF_LINE-1);
	}
	return res;
}

//===========================================================

INLINE void SkipNLines(FILE *file,int n)
// skips n lines from the file starting from current position in a file
{
	while (n>0) {
		SkipFullLine(file);
		n--;
	}
}

//===========================================================

static int SkipComments(FILE *file)
/* skips comments (#...), all lines, starting from current position in a file.
 * returns number of lines skipped
 */
{
	int lines=0,ch;

	while ((ch=fgetc(file))=='#') {
		SkipFullLine(file);
		lines++;
	}
	if (ch!=EOF) ungetc(ch,file);

	return lines;
}

//===========================================================
#define DDSCAT_HL 6 // number of header lines in DDSCAT format

static void InitDipFile(const char *fname,int *bX,int *bY,int *bZ,int *Nm)
/* read dipole file first to determine box sizes and Nmat; input is not checked for very large
 * numbers (integer overflows) to increase speed; this function opens file for reading, the file is
 * closed in ReadDipFile.
 */
{
	int x,y,z,mat,line,scanned,mustbe,skiplines,anis_warned;
	long tl; // dumb variable
	int t2,t3; // dumb variables
	int maxX,maxY,maxZ,maxN;
	char formtext[MAX_LINE];

	dipfile=FOpenErr(fname,"r",ALL_POS);
	read_format=UNDEF;
	/* test for DDSCAT format; in not-DDSCAT format, the line scanned below may be a long comment;
	 * therefore we first skip all comments
	 */
	line=SkipComments(dipfile);
	if (line<=DDSCAT_HL) {
		SkipNLines(dipfile,DDSCAT_HL-line);
		if (FgetsError(dipfile,fname,&line,POSIT)!=NULL
			&& sscanf(linebuf,ddscat_format,&tl,&x,&y,&z,&mat,&t2,&t3)==7) {
				read_format=SF_DDSCAT;
				strcpy(formtext,"DDSCAT format (FRMFIL)");
				mustbe=7;
				line=DDSCAT_HL;
				fseek(dipfile,0,SEEK_SET);
				SkipNLines(dipfile,line);
		}
	}
	// if format is not yet determined, test for ADDA text formats
	if (read_format==UNDEF) {
		fseek(dipfile,0,SEEK_SET);
		line=SkipComments(dipfile);
		/* scanf and analyze Nmat; if there is blank line between comments and Nmat, it fails later;
		 * the value of Nmat obtained here is not actually relevant, the main factor is maximum
		 * domain number among all dipoles.
		 */
		scanned=fscanf(dipfile,"Nmat=%d\n",Nm);
		if (scanned==EOF) LogError(EC_ERROR,ONE_POS,"No dipole positions are found in %s",fname);
		else if (scanned==0) { // no "Nmat=..."
			read_format=SF_TEXT;
			strcpy(formtext,"ADDA text format (single domain)");
			*Nm=1;
			mustbe=3;
		}
		else { // "Nmat=..." present
			read_format=SF_TEXT_EXT;
			strcpy(formtext,"ADDA text format (multi-domain)");
			mustbe=4;
			line++;
		}
	}
	// scan main part of the file
	skiplines=line;
	maxX=maxY=maxZ=INT_MIN;
	minX=minY=minZ=INT_MAX;
	maxN=1;
	anis_warned=FALSE;
	// reading is performed in lines
	while(FgetsError(dipfile,fname,&line,POSIT)!=NULL) {
		// scan numbers in a line
		if (read_format==SF_TEXT) scanned=sscanf(linebuf,geom_format,&x,&y,&z);
		else if (read_format==SF_TEXT_EXT) scanned=sscanf(linebuf,geom_format_ext,&x,&y,&z,&mat);
		// for ddscat format, only first material is used, other two are ignored
		else if (read_format==SF_DDSCAT) {
			scanned=sscanf(linebuf,ddscat_format,&tl,&x,&y,&z,&mat,&t2,&t3);
			if (!anis_warned && (t2!=mat || t3!=mat)) {
				LogError(EC_WARN,ONE_POS,"Anisotropic dipoles are detected in file %s (first on "
					"line %d). ADDA ignores this anisotropy, using only the identifier of "
					"x-component of refractive index as domain number",fname,line);
				anis_warned=TRUE;
			}
		}
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			if (scanned!=mustbe) // this in most cases indicates wrong format
				LogError(EC_ERROR,ONE_POS,"%s was detected, but error occurred during scanning "
					"of line %d from dipole file %s",formtext,line,fname);
			if (read_format!=SF_TEXT) {
				if (mat<=0) LogError(EC_ERROR,ONE_POS,"%s was detected, but nonpositive material "
					"number (%d) encountered during scanning of line %d from dipole file %s",
					formtext,mat,line,fname);
				else if (mat>maxN) maxN=mat;
			}
			// update maxima and minima
			if (x>maxX) maxX=x;
			if (x<minX) minX=x;
			if (y>maxY) maxY=y;
			if (y<minY) minY=y;
			if (z>maxZ) maxZ=z;
			if (z<minZ) minZ=z;
		}
	}
	if (read_format==SF_TEXT_EXT) {
		if (*Nm!=maxN) LogError(EC_WARN,ONE_POS,"Nmat (%d), as given in %s, is not equal to the "
				"maximum domain number (%d) among all specified dipoles; hence the former is "
				"ignored",*Nm,fname,maxN);
	}
	*Nm=maxN;
	// set grid (box) sizes
	*bX=jagged*(maxX-minX+1);
	*bY=jagged*(maxY-minY+1);
	*bZ=jagged*(maxZ-minZ+1);
	// not optimal way, but works more robustly when non-system EOL is used in data file
	fseek(dipfile,0,SEEK_SET);
	SkipNLines(dipfile,skiplines);
}
#undef DDSCAT_HL
//===========================================================

static void ReadDipFile(const char *fname)
/* read dipole file; no consistency checks are made since they are made in InitDipFile.
 * the file is opened in InitDipFile; this function only closes the file.
 */
{
	int x,y,z,x0,y0,z0,mat,scanned;
	long tl; // dumb variable
	int t2,t3; // dumb variables
	int index;
	size_t boxXY,boxX_l;

	// to remove possible overflows
	boxX_l=(size_t)boxX;
	boxXY=boxX_l*boxY;

	mat=1;
	while(fgets(linebuf,BUF_LINE,dipfile)!=NULL) {
		// scan numbers in a line
		if (read_format==SF_TEXT) scanned=sscanf(linebuf,geom_format,&x0,&y0,&z0);
		else if (read_format==SF_TEXT_EXT) scanned=sscanf(linebuf,geom_format_ext,&x0,&y0,&z0,&mat);
		else if (read_format==SF_DDSCAT)
			scanned=sscanf(linebuf,ddscat_format,&tl,&x0,&y0,&z0,&mat,&t2,&t3);
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			// shift dipole position to be nonnegative
			x0-=minX;
			y0-=minY;
			z0-=minZ;
			// initialize box jagged*jagged*jagged instead of one dipole
			for (z=jagged*z0;z<jagged*(z0+1);z++) if (z>=local_z0 && z<local_z1_coer)
				for (x=jagged*x0;x<jagged*(x0+1);x++) for (y=jagged*y0;y<jagged*(y0+1);y++) {
					index=(z-local_z0)*boxXY+y*boxX_l+x;
					material_tmp[index]=(unsigned char)(mat-1);
			}
		}
	}
	FCloseErr(dipfile,fname,ALL_POS);
}

//==========================================================
#define ALLOCATE_SEGMENTS(N) (struct segment *)voidVector((N)*sizeof(struct segment),ALL_POS,\
	"contour segment");

void InitContourSegment(struct segment *seg,const bool increasing)
/* recursively initialize a segment of a contour: allocates memory and calculates all elements
 * some elements are calculated during the forward sweep (from long segments to short), others -
 * during the backward sweep.
 * Recursive function calls incurs certain overhead, however here it is not critical.
 */
{
	int i;
	struct segment *s1,*s2;

	/* Remove constant parts in the beginning and end of segment, if present.
	 * After this procedure first is guaranteed to be less than last by definition of the segment
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

//===========================================================
#define CHUNK_SIZE 128 // how many numbers are allocated at once for adjustable arrays

static void InitContour(const char *fname,double *ratio,double *shSize)
/* Reads a contour from the file, rotates it so that it starts from a local minimum in ro,
 * then divides it into monotonic (over ro) segments. It produces data, which are later used to
 * test each dipole for being inside the contour. Segments are either increasing or non-decreasing.
 */
{
	int line; // current line number
	int nr;   // number of contour points read from the file
	int size; // current size of the allocated memory for contour
	int i,j,scanned;
	double *bufRo,*bufZ; // temporary buffers
	int *index;
	double ro,z,romin,romax,zmin,zmax,mult,zmid;
	FILE* file;
	bool increasing;

	D("InitContour has started");
	// Read contour from file
	file=FOpenErr(fname,"r",ALL_POS);
	line=SkipComments(file);
	size=CHUNK_SIZE;
	MALLOC_VECTOR(bufRo,double,size,ALL);
	MALLOC_VECTOR(bufZ,double,size,ALL);
	nr=0;
	// reading is performed in lines
	while(FgetsError(file,fname,&line,POSIT)!=NULL) {
		// scan numbers in a line
		scanned=sscanf(linebuf,"%lf %lf",&ro,&z);
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			if (scanned!=2) // this in most cases indicates wrong format
				LogError(EC_ERROR,ONE_POS,"Error occurred during scanning of line %d from contour "
					"file %s",line,fname);
			// check for consistency of input
			if (ro<0) LogError(EC_ERROR,ONE_POS,"Negative ro-coordinate is found on line %d in "
				"contour file %s",line,fname);
			// update extreme values
			if (nr==0) {
				zmax=zmin=z;
				romax=romin=ro;
			}
			else {
				if (z>zmax) zmax=z;
				if (z<zmin) zmin=z;
				if (ro>romax) romax=ro;
				if (ro<romin) romin=ro;
			}
			// add allocated memory to buf, if needed
			if (nr >= size) {
				size+=CHUNK_SIZE;
				REALLOC_DVECTOR(bufRo,size,ALL);
				REALLOC_DVECTOR(bufZ,size,ALL);
			}
			bufRo[nr]=ro;
			bufZ[nr]=z;
			nr++;
		}
	}
	FCloseErr(file,fname,ALL_POS);
	// Check number of points read
	if (nr<3) LogError(EC_ERROR,ONE_POS,
		"Contour from file %s contains less than three points",fname);

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
	else if (i==nr-1 && bufRo[nr-1]==bufRo[0]) LogError(EC_ERROR,ONE_POS,
		"Contour from file %s has zero area. Hence the scatterer is void",fname);
	/* Construct working contour so that its first point = last and is a local minimum. It is done
	 * by rotating buf and adding one extra point. Then free the buffer.
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
	contRoSqMin=romin*romin*mult*mult;
	for (i=0;i<=nr;i++) {
		contRo[i]*=mult;
		contZ[i]=(contZ[i]-zmid)*mult;
	}

	/* divide the contour into the segments; actually only the index is constructed marking end
	 * points of the segments
	 */
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
	/* Calculate maximum and minimum ro for segments;
	 * We implicitly use that first segment is increasing, second - decreasing, and so on.
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
	D2("Nseg=%d",contNseg);
	D2("minroSq=%g",contRoSqMin);
}
#undef CHUNK_SIZE
#undef ALLOCATE_SEGMENTS

//==========================================================

bool CheckContourSegment(struct segment *seg)
/* Checks, whether point is under or above the segment, by traversing the tree of segments.
 * It returns true, if intersecting z value is larger than given z, and false otherwise. Point is
 * defined by local variables contCurRo and contCurZ.
 */
{
	while (true) {
		if (contCurZ < seg->zmin) return true;
		else if (contCurZ > seg->zmax) return false;
		else if (seg->single) return (contCurZ < seg->add + contCurRo*seg->slope);
		else seg=(contCurRo<seg->romid ? seg->left : seg->right);
	}
}

//==========================================================

void FreeContourSegment(struct segment *seg)
/* recursively frees memory allocated for contour segments
 * Recursive function calls incurs certain overhead, however here it is not critical.
 */
{
	if (!(seg->single)) {
		FreeContourSegment(seg->left);
		FreeContourSegment(seg->right);
	}
}

//==========================================================

#define KEY_LENGTH 2            // length of key for initialization of random generator
#define MAX_ZERO_FITS 1E4       // maximum number of zero fits in a row (each - many granules)
#define MAX_FALSE_SKIP 10       // number of false skips in granule placement to complete the set
#define MAX_FALSE_SKIP_SMALL 10 // the same for small granules
#define MAX_GR_SET USHRT_MAX    // maximum size of granule set
#define MIN_CELL_SIZE 4.0       // minimum cell size for small granules

INLINE int CheckCell(const double *gr,const double *vgran,const unsigned short *tree_index,
                     const double Di2,const int start,int *fits)
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
		if ((t1*t1+t2*t2+t3*t3)<Di2) *fits=FALSE;
		index=tree_index[index];
	}
	return last;
}
//==========================================================

static double PlaceGranules(void)
/* Randomly places granules inside the specified domain;
 * Mersenne Twister is used for generating random numbers
 *
 * A simplest algorithm is used: to place randomly a sphere, and see if it overlaps with any
 * dipoles (more exactly: centers of dipoles) of not correct domain; if not, accept it and
 * fill all this dipoles with granules' domain. Optimized to perform in two steps:
 * First it places of set of not-intersecting granules and do only a quick check against the
 * "domain pattern" - coarse representation of the domain. On the second step granules of the
 * whole set are thoroughly checked against the whole domain on each processor. When small
 * granules are used, no domain pattern is used - makes it simpler.
 * Intersection of two granules between the sets is checked only through dipoles, which is not
 * exact, however it allows considering arbitrary complex domains, which is described only by
 * a set of occupied dipoles.
 * This algorithm is unsuitable for high volume fractions, it becomes very slow and for some
 * volume fractions may fail at all (Metropolis algorithm should be more suitable, however it
 * is hard to code for arbitrary domains). Moreover, statistical properties of the obtained
 * granules distribution may be not perfect, however it seems good enough for our applications.
 *
 * Currently it is not working with jagged. That should be improved, by rewriting the jagged
 * calculation throughout the program
 */

{
	int i,j,k,zerofit,last;
	size_t n,count,count_gr,false_count,ui;
	size_t boxXY;
	double nd;                           // number of dipoles occupied by granules
	int index,index1,index2;             // indices for dipole grid
	int dom_index,dom_index1,dom_index2; // indices for auxiliary grid
	int gX,gY,gZ;                        // auxiliary grid dimensions
	size_t gXY,gr_gN;                    // ... and their products
	size_t avail;                        // number of available (free) domain cells
	int gX2,gY2,gZ2,locgZ2;
	int i0,i1,j0,j1,k0,k1;
	int fits;
	int cur_Ngr,ig,max_Ngr; // number of granules in a current set, index, and maximum set size
	double gdX,gdY,gdZ,gdXh,gdYh,gdZh; // auxiliary grid cell sizes and their halfs (h)
	int locz0,locz1,locgZ,gr_locgN;
	double R,R2,Di,Di2;          // radius and diameter of granule, and their squares
	double x0,x1,y0,y1,z0,z1;    // where to put random number (inner box)
	int id0,id1,jd0,jd1,kd0,kd1; // dipoles limit that fall inside inner box
	int Nfit;        // number of successfully placed granules in a current set
	double overhead; // estimate of the overhead needed to have exactly needed N of granules
	double tmp1,tmp2,tmp3,t1,t2,t3;
	int sx,sy,sz; /* maximum shifts for checks of neighboring cells in auxiliary grid
	               *  for 'small' it is the shift in index
	               */
	unsigned long key[KEY_LENGTH]; // key to initialize random number generator
	unsigned char *dom;            // information about the domain on a granule grid
	unsigned short *occup;         // information about the occupied cells
	int sm_gr;                     // whether granules are small (then simpler algorithm is used)
	unsigned short *tree_index;    // index for traversing granules inside one cell (for small)
	double *vgran;                 // coordinates of a set of granules
	char *vfit;                    // results of granule fitting on the grid (boolean)
	int *ginX,*ginY,*ginZ;         // indices to find dipoles inside auxiliary grid
	int indX,indY,indZ;            // indices for doubled auxiliary grid
	int bit;                       // bit position in char of 'dom'
	double gr[3];                  // coordinates of a single granule
	FILE *file;                    // file for saving granule positions
	char fname[MAX_FNAME];         // filename of file

	// prepare granule file for saving if needed
	if (store_grans && ringid==ROOT) {
		sprintf(fname,"%s/" F_GRANS,directory);
		file=FOpenErr(fname,"w",ONE_POS);
		fprintf(file,"#generated by ADDA v." ADDA_VERSION "\n"
		             "#granule diameter = %.10g\n",gr_d);
	}
	// set variables; consider jagged
	Di=gr_d/(gridspace*jagged);
	if (Di<1) LogError(EC_WARN,ONE_POS,"Granule diameter is smaller than dipole size. It is "
		"recommended to increase resolution");
	R=Di/2;
	R2=R*R;
	Di2=4*R2;
	boxXY=boxX*(size_t)boxY;
	// inner box
	if (Di>MIN(boxX,MIN(boxY,boxZ))) LogError(EC_WARN,ONE_POS,
		"Granule size is larger than minimum particle dimension");
	x0=R-0.5;
	x1=boxX-R-0.5;
	y0=R-0.5;
	y1=boxY-R-0.5;
	z0=R-0.5;
	z1=boxZ-R-0.5;
	// initialize auxiliary grid
	CheckOverflow(MAX(boxX,MAX(boxY,boxZ))*10/Di,ONE_POS,"PlaceGranules()");
	tmp1=sqrt(3)/Di;
	gX=(int)ceil((x1-x0)*tmp1);
	gdX=(x1-x0)/gX;
	gY=(int)ceil((y1-y0)*tmp1);
	gdY=(y1-y0)/gY;
	gZ=(int)ceil((z1-z0)*tmp1);
	gdZ=(z1-z0)/gZ;
	sm_gr=(gdX<2 || gdY<2 || gdZ<2); // sets the discrimination for small or large granules
	if (sm_gr) {
		PRINTZ("Using algorithm for small granules\n");
		// redefine auxiliary grid
		tmp1=1/MAX(2*Di,MIN_CELL_SIZE);
		gX=(int)floor((x1-x0)*tmp1);
		gdX=(x1-x0)/gX;
		gY=(int)floor((y1-y0)*tmp1);
		gdY=(y1-y0)/gY;
		gZ=(int)floor((z1-z0)*tmp1);
		gdZ=(z1-z0)/gZ;
	}
	else {
		PRINTZ("Using algorithm for large granules\n");
		gX2=2*gX;
		gdXh=gdX/2;
		gY2=2*gY;
		gdYh=gdY/2;
		gZ2=2*gZ;
		gdZh=gdZ/2;
		/* this sets maximum distance of neighboring cells to check; condition gdX<R can only occur
		 * if gX<=7, which is quite rare, so no optimization is performed. sx>3 can only occur if
		 * gX<=2 and then it doesn't make sense to take bigger sx. Absolutely analogous for y, z.
		 */
		if (gdX<R) sx=3;
		else sx=2;
		if (gdY<R) sy=3;
		else sy=2;
		if (gdZ<R) sz=3;
		else sz=2;
	}
	gXY=MultOverflow(gX,gY,ONE_POS,"PlaceGranules()");
	gr_gN=MultOverflow(gXY,gZ,ONE_POS,"PlaceGranules()");
	// calculate maximum number of granules in a grid; crude estimate
	tmp2=(ceil((x1-x0)/Di)+1)*(ceil((y1-y0)/Di)+1)*(ceil((z1-z0)/Di)+1);
	max_Ngr=MIN(MAX_GR_SET,tmp2);
	// local z grid + initialize communications
	SetGranulComm(z0,z1,gdZ,gZ,gXY,max_Ngr,&locz0,&locz1,sm_gr);
	if (!sm_gr) {
		locgZ=locz1-locz0;
		locgZ2=2*locgZ;
		gr_locgN=gXY*locgZ;
	}
	if (ringid==ROOT) {
		// initialize random generator
		key[0]=(unsigned long)time(NULL);
		key[1]=(unsigned long)(clock()-wt_start);
		init_by_array(key,KEY_LENGTH);
		// allocate memory
		MALLOC_VECTOR(occup,ushort,gr_gN,ONE);
		if (sm_gr) MALLOC_VECTOR(tree_index,ushort,max_Ngr,ONE);
		else MALLOC_VECTOR(dom,uchar,gr_gN,ALL);
	}
	else if (!sm_gr && locgZ!=0) MALLOC_VECTOR(dom,uchar,gr_locgN,ALL);
	MALLOC_VECTOR(vgran,double,3*max_Ngr,ALL);
	MALLOC_VECTOR(vfit,char,max_Ngr,ALL);
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
		jd1=ginY[gZ2];
		for (i=0;i<=locgZ2;i++) ginZ[i]=(int)ceil(z0+(i+2*locz0)*gdZh);
		kd0=MAX(ginZ[0],local_z0);
		indZ=1;
		if (kd0>=ginZ[1]) indZ++;
		kd1=MIN(ginZ[locgZ2],local_z1_coer);
	}
	n=count=count_gr=false_count=0;
	nd=0;
	// crude estimate of the probability to place a small granule into domain
	if (sm_gr) overhead=Ndip/mat_count[gr_mat];
	else overhead=1;
	// main cycle
	while (n<gr_N) {
		if (sm_gr) { // small granules
			// just generate granules
			if (ringid==ROOT) {
				cur_Ngr=MIN(ceil((gr_N-n)*overhead),max_Ngr);
				// generate points and quick check
				ig=false_count=0;
				for (ui=0;ui<gr_gN;ui++) occup[ui]=MAX_GR_SET; // used as undefined
				while (ig<cur_Ngr) {
					count++;
					false_count++;
					fits=TRUE;
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
					last=CheckCell(gr,vgran,tree_index,Di2,occup[index],&fits);
					// weird construction (7 inclosed 'if' structures) but should be fast
					if (fits) {
						// possible x-neighbor
						t1*=gdX; // transform shifts to usual coordinates; done only when needed
						sx=0;
						if (t1<Di) {
							if (indX!=0) sx=-1;
						}
						else if ((t1=gdX-t1)<Di && indX!=gX-1) sx=1;
						if (sx!=0) CheckCell(gr,vgran,tree_index,Di2,occup[index+sx],&fits);
						if (fits) {
							// possible y-neighbor
							t2*=gdY;
							sy=0;
							if (t2<Di) {
								if (indY!=0) sy=-gX;
							}
							else if ((t2=gdY-t2)<Di && indY!=gY-1) sy=gX;
							if (sy!=0) CheckCell(gr,vgran,tree_index,Di2,occup[index+sy],&fits);
							if (fits) {
								// possible z-neighbor
								t3*=gdZ;
								sz=0;
								if (t3<Di) {
									if (indZ!=0) sz=-(int)gXY;
								}
								else if ((t3=gdZ-t3)<Di && indZ!=gZ-1) sz=gXY;
								if (sz!=0) CheckCell(gr,vgran,tree_index,Di2,occup[index+sz],&fits);
								if (fits) {
									// possible xy-neighbor
									if (sx!=0 && sy!=0 && ((tmp1=t1*t1)+(tmp2=t2*t2)<Di2))
										CheckCell(gr,vgran,tree_index,Di2,occup[index+sx+sy],&fits);
									if (fits) {
										// possible xz-neighbor
										if (sx!=0 && sz!=0 && ((tmp1+(tmp3=t3*t3))<Di2))
											CheckCell(gr,vgran,tree_index,Di2,
												occup[index+sx+sz],&fits);
										if (fits) {
											// possible yz-neighbor & xyz-neighbor
											if (sy!=0 && sz!=0 && (tmp2+tmp3<Di2)) {
												CheckCell(gr,vgran,tree_index,Di2,
													occup[index+sy+sz],&fits);
												if (fits && sx!=0 && (tmp1+tmp2+tmp3<Di2))
													CheckCell(gr,vgran,tree_index,Di2,
														occup[index+sx+sy+sz],&fits);
											}
										}
									}
								}
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
							if (material_tmp[index]!=gr_mat)
								dom[dom_index]|=(unsigned char)(1<<bit);
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
			// send/collect domain pattern
			CollectDomainGranul(dom,gXY,locz0,locgZ,&Timing_Granul_comm);
			if (ringid==ROOT) {
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
					indZ=(int)floor(gr[2]); // position bit inside one cell
					bit=1<<((indX&1)+((indY&1)<<1)+((indZ&1)<<2));
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
						fits=TRUE;
						false_count++;
						if ((i0=indX-sx)<0) i0=0;
						if ((i1=indX+sx+1)>gZ) i1=gX;
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
										fits=FALSE;
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
							/* Here it is possible to correct the domain pattern because of the
							 * presence of a new granule. However it probably will be useful only
							 * for large volume fractions
							 */
						}
						if (false_count>MAX_FALSE_SKIP) break;
					}
				}
				// real number of placed granules for this set
				cur_Ngr=ig;
			}
		} // end of large granules
		// cast to all processors
		MyBcast(&cur_Ngr,int_type,1,&Timing_Granul_comm);
		MyBcast(vgran,double_type,3*cur_Ngr,&Timing_Granul_comm);
		count_gr+=cur_Ngr;
		// final check if granules belong to the domain
		for (ig=0;ig<cur_Ngr;ig++) {
			memcpy(gr,vgran+3*ig,3*sizeof(double));
			k0=MAX((int)ceil(gr[2]-R),local_z0);
			k1=MIN((int)floor(gr[2]+R),local_z1_coer-1);
			fits=TRUE;
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
							fits=FALSE;
							break;
						}
					}
					if (!fits) break;
				}
				if (!fits) break;
			}
			vfit[ig]=(char)fits;
		}
		// collect fits
		ExchangeFits(vfit,cur_Ngr,&Timing_Granul_comm);
		// fit dipole grid with successive granules
		Nfit=n;
		for (ig=0;ig<cur_Ngr;ig++) {
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
				// if the allocation was too optimistic
				if (n>=gr_N) break;
			}
		}
		// save correct granule positions to file
		if (store_grans && ringid==ROOT) for (ig=0;ig<cur_Ngr;ig++) if (vfit[ig])
			fprintf(file,"%.10g %.10g %.10g\n",gridspace*(vgran[3*ig]-cX),
				gridspace*(vgran[3*ig+1]-cY),gridspace*(vgran[3*ig+2]-cZ));
		Nfit=n-Nfit;
		/* overhead is estimated based on the estimation of mean value - 1*standard deviation
		 * for the probability of fitting one granule. It is estimated from the Bernoulli statistics
		 * k out of n successful hits. M(p)=(k+1)/(n+2); s^2(p)=(k+1)(n-k+1)/(n+3)(n+2)^2
		 * M(p)-s(p)=[(k+1)/(n+2)]*[1-sqrt((n-k+1)/(k+1)(n+3))];
		 * overhead=1/latter.
		 */
		overhead=(cur_Ngr+2)/((1-sqrt((cur_Ngr-Nfit+1)/(double)((Nfit+1)*(cur_Ngr+3))))*(Nfit+1));
		if (Nfit!=0) zerofit=0;
		else {
			zerofit++;
			// check if taking too long
			if (zerofit>MAX_ZERO_FITS) {
				MyInnerProduct(&nd,double_type,1,&Timing_Granul_comm);
				LogError(EC_ERROR,ONE_POS,"The granule generator failed to reach required volume "
					"fraction (%g) of granules. %zu granules were successfully placed up to a "
					"volume fraction of %g.",gr_vf,n,nd/mat_count[gr_mat]);
			}
		}
	}
	/* conversions to (unsigned long) are needed (to remove warnings) because %z printf argument is
	 * not yet supported by all target compiler environments
	 */
	PRINTZ("Granule generator: total random placements= %lu (efficiency 1 = %g)\n"
	       "                   possible granules= %lu (efficiency 2 = %g)\n",
	       (unsigned long)count,count_gr/(double)count,(unsigned long)count_gr,
	       gr_N/(double)count_gr);
	MyInnerProduct(&nd,double_type,1,&Timing_Granul_comm);
	// free everything
	if (ringid==ROOT) {
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
	if (store_grans && ringid==ROOT) {
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
//==========================================================

static int FitBox(const int box)
/* finds the smallest value for which program would work (should be even and divide jagged);
 * the limit is also checked
 */
{
	int res;

	if (IS_EVEN(jagged)) res=jagged*((box+jagged-1)/jagged);
	else res=2*jagged*((box+2*jagged-1)/(2*jagged));
	if (res>BOX_MAX) LogError(EC_ERROR,ONE_POS,
		"Derived grid size (%d) is too large (>%d)",res,BOX_MAX);
	return res;
}

//==========================================================

void InitShape(void)
/* perform of initialization of symmetries and boxY, boxZ. Estimate the volume of the particle, when
 * not discretized. Check whether enough refractive indices are specified.
 */
{
	int n_boxX,n_boxY,n_boxZ; // new values for dimensions
	double n_sizeX; // new value for size
	double h_d,b_d,c_d,h2,b2,c2;
	double yx_ratio,zx_ratio,tmp1,tmp2,tmp3;
	double diskratio,aspectY,aspectZ;
	double ad,ct,ct2; // cos(theta0) and its square
	TIME_TYPE tstart;
	int Nmat_need,i,temp;
	int dpl_def_used; // if default dpl is used for grid initialization
	bool box_det_sh; // if boxX is determined by shape itself
	bool size_det_sh; // if size is determined by shape itself
	bool size_given_cmd; // if size is given in the command line
	char sizename[MAX_LINE]; // type of input size, used in diagnostic messages
	/* TO ADD NEW SHAPE
	 * Add here all intermediate variables, which are used only inside this function. You may as
	 * well use 'tmp1'-'tmp3' variables defined above.
	 */

	tstart=GET_TIME();

	box_det_sh=(shape==SH_READ);
	size_det_sh=(shape==SH_AXISYMMETRIC);
	/* TO ADD NEW SHAPE
	 * If new shape defines dimension of the computational grid or absolute size of the particle,
	 * change corresponding definition in one of two lines above. In many cases this is not
	 * relevant.
	 */

	size_given_cmd=(sizeX!=UNDEF || a_eq!=UNDEF);
	if (sizeX!=UNDEF) strcpy(sizename,"size");
	else if (a_eq!=UNDEF) strcpy(sizename,"eq_rad");
	// check for redundancy of input data
	if (dpl!=UNDEF) {
		if (size_given_cmd) {
			if (boxX!=UNDEF) PrintError("Extra information is given by setting '-dpl', '-grid', "
				"and '-%s'",sizename);
			else if (box_det_sh) PrintError("Extra information is given by setting both '-dpl' and "
				"'-%s', while shape '%s' sets the size of the grid",sizename,shapename);
		}
		else if (size_det_sh) {
			if (boxX!=UNDEF) PrintError("Extra information is given by setting '-dpl' and '-grid', "
				"while shape '%s' sets the particle size",shapename);
			// currently this can't happen, but may become relevant in the future
			else if (box_det_sh) PrintError("Extra information is given by setting '-dpl', while "
				"shape '%s' sets both the particle size and the size of the grid",shapename);
		}
	}
	/* calculate default dpl - 10*sqrt(max(|m|));
	 * for anisotropic each component is considered separately
	 */
	tmp2=0;
	for (i=0;i<Ncomp*Nmat;i++) {
		tmp1=cAbs2(ref_index[i]);
		if (tmp2<tmp1) tmp2=tmp1;
	}
	dpl_def=10*sqrt(tmp2);
	dpl_def_used=FALSE;
	// initialization of global option index for error messages
	opt=opt_sh;
	// shape initialization
	if (shape==SH_AXISYMMETRIC) {
		/* Axisymmetric homogeneous shape, defined by its contour in ro-z plane of the cylindrical
		 * coordinate system. Its symmetry axis coincides with the z-axis, and the contour is read
		 * from file. Each line defines ro and z coordinates of a point, the first and the last
		 * points are connected automatically. Linear interpolation is used between the points.
		 */
		SPRINTZ(sh_form_str,"axisymmetric, defined by a contour in ro-z plane from file %s;"
			" diameter:%%.10g",shape_fname);
		InitContour(shape_fname,&zx_ratio,&n_sizeX);
		yx_ratio=1;
		symZ=FALSE; // input contour is assumed asymmetric over ro-axis
		/* TODO: this can be determined from the contour. However, it is not trivial, especially
		 * when the contour intersects itself.
		 */
		volume_ratio=UNDEF;
		Nmat_need=1;
	}
	else if (shape==SH_BOX) {
		if (sh_Npars==0) {
			STRCPYZ(sh_form_str,"cube; size of edge along x-axis:%.10g");
			aspectY=aspectZ=1;
		}
		else { // 2 parameters are given
			aspectY=sh_pars[0];
			TestPositive(aspectY,"aspect ratio y/x");
			aspectZ=sh_pars[1];
			TestPositive(aspectZ,"aspect ratio z/x");
			SPRINTZ(sh_form_str,"rectangular parallelepiped; size along x-axis:%%.10g, aspect "
				"ratios y/x=%.10g, z/x=%.10g",aspectY,aspectZ);
		}
		if (aspectY!=1) symR=FALSE;
		// set half-aspect ratios
		haspY=aspectY/2;
		haspZ=aspectZ/2;
		volume_ratio=aspectY*aspectZ;
		yx_ratio=aspectY;
		zx_ratio=aspectZ;
		Nmat_need=1;
	}
	else if(shape==SH_CAPSULE) {
		diskratio=sh_pars[0];
		TestNonNegative(diskratio,"height to diameter ratio");
		SPRINTZ(sh_form_str,"capsule; diameter(d):%%.10g, cylinder height h/d=%.10g",diskratio);
		hdratio=diskratio/2;
		volume_ratio = PI_OVER_FOUR*diskratio + PI_OVER_SIX;
		yx_ratio=1;
		zx_ratio=diskratio+1;
		Nmat_need=1;
	}
	else if (shape==SH_COATED) {
		coat_ratio=sh_pars[0];
		TestRangeII(coat_ratio,"inner/outer diameter ratio",0,1);
		SPRINTZ(sh_form_str,"coated sphere; diameter(d):%%.10g, inner diameter d_in/d=%.10g",
			coat_ratio);
		if (sh_Npars==4) {
			coat_x=sh_pars[1];
			coat_y=sh_pars[2];
			coat_z=sh_pars[3];
			if (coat_x*coat_x+coat_y*coat_y+coat_z*coat_z>0.25*(1-coat_ratio)*(1-coat_ratio))
				PrintErrorHelp("Inner sphere is not fully inside the outer");
			SPRINTZ(sh_form_str+strlen(sh_form_str),
			"\n       position of inner sphere center r/d= {%.10g,%.10g,%.10g}",
				coat_x,coat_y,coat_z);
		}
		else coat_x=coat_y=coat_z=0; // initialize default values
		coat_r2=0.25*coat_ratio*coat_ratio;
		volume_ratio=PI_OVER_SIX;
		if (coat_x!=0) symX=symR=FALSE;
		if (coat_y!=0) symY=symR=FALSE;
		if (coat_z!=0) symZ=FALSE;
		yx_ratio=zx_ratio=1;
		Nmat_need=2;
	}
	else if(shape==SH_CYLINDER) {
		diskratio=sh_pars[0];
		TestPositive(diskratio,"height to diameter ratio");
		SPRINTZ(sh_form_str,"cylinder; diameter(d):%%.10g, height h/d=%.10g",diskratio);
		hdratio=diskratio/2;
		volume_ratio=PI_OVER_FOUR*diskratio;
		yx_ratio=1;
		zx_ratio=diskratio;
		Nmat_need=1;
	}
	else if (shape==SH_EGG) {
		/* determined by equation: (a/r)^2=1+nu*cos(theta)-(1-eps)cos^2(theta)
		 * or equivalently: a^2=r^2+nu*r*z-(1-eps)z^2. Parameters must be 0<eps<=1, 0<=nu<eps.
		 * This shape is proposed in: Hahn D.V., Limsui D., Joseph R.I., Baldwin K.C., Boggs N.T.,
		 * Carr A.K., Carter C.C., Han T.S., and Thomas M.E. "Shape characteristics of biological
		 * spores", paper 6954-31 to be presented at "SPIE Defence + Security", March 2008
		 */
		egeps=sh_pars[0];
		TestRangeNI(egeps,"egg parameter epsilon",0,1);
		egnu=sh_pars[1];
		TestRangeIN(egnu,"egg parameter nu",0,egeps);
		// egg shape is symmetric about z-axis (xz and yz planes, but generally NOT xy plane)
		if (egnu!=0) symZ=FALSE;
		/* cos(theta0): ct=-nu/[eps+sqrt(eps^2-nu^2)]; this expression for root of the quadratic
		 * equation is used for numerical stability (i.e. when nu=0); at this theta0 the diameter
		 * (maximum width perpendicular to z-axis) d=Dx is located
		 */
		ct=-egnu/(egeps+sqrt(egeps*egeps-egnu*egnu));
		ct2=ct*ct;
		// Determine ad=(a/d) and its square
		ad2=(1+egnu*ct-(1-egeps)*ct2)/(4*(1-ct2));
		ad=sqrt(ad2);
		tmp1=1/sqrt(egeps+egnu);
		tmp2=1/sqrt(egeps-egnu);
		tmp3=2*(1-egeps);
		/* Center of the computational box (z coordinate):
		 * z0=(a/d)*[1/sqrt(eps+nu)+1/sqrt(eps-nu)]/2; but more numerically stable expression is
		 * used (for nu->0). Although it may overflow faster for nu->eps, volume_ratio (below) will
		 * overflow even faster. It is used to shift coordinates from the computational reference
		 * frame (centered at z0) to the natural one
		 */
		egz0=-ad*egnu*(tmp1*tmp1*tmp2*tmp2)/(tmp1+tmp2);
		/* (V/d^3)=(4*pi/3)*(a/d)^3*{[2(1-eps)-nu]/sqrt(eps+nu)+[2(1-eps)+nu]/sqrt(eps-nu)}/
		 *        /[nu^2+4(1-eps)]
		 */
		volume_ratio=FOUR_PI_OVER_THREE*ad2*ad*((tmp3-egnu)*tmp1+(tmp3+egnu)*tmp2)
		            /(egnu*egnu+2*tmp3);
		SPRINTZ(sh_form_str,"egg; diameter(d):%%.10g, epsilon=%.10g, nu=%.10g, a/d=%.10g",
			egeps,egnu,ad);
		Nmat_need=1;
		yx_ratio=1;
		zx_ratio=ad*(tmp1+tmp2); // (a/d)*[1/sqrt(eps+nu)+1/sqrt(eps-nu)]
	}
	else if (shape==SH_ELLIPSOID) {
		aspectY=sh_pars[0];
		TestPositive(aspectY,"aspect ratio y/x");
		aspectZ=sh_pars[1];
		TestPositive(aspectZ,"aspect ratio z/x");
		SPRINTZ(sh_form_str,"ellipsoid; size along x-axis:%%.10g, aspect ratios y/x=%.10g, "
			"z/x=%.10g",aspectY,aspectZ);
		if (aspectY!=1) symR=FALSE;
		// set inverse squares of aspect ratios
		invsqY=1/(aspectY*aspectY);
		invsqZ=1/(aspectZ*aspectZ);
		volume_ratio=PI_OVER_SIX*aspectY*aspectZ;
		yx_ratio=aspectY;
		zx_ratio=aspectZ;
		Nmat_need=1;
	}
	else if (shape==SH_LINE) {
		STRCPYZ(sh_form_str,"line; length:%g");
		symY=symZ=symR=FALSE;
		n_boxY=n_boxZ=jagged;
		yx_ratio=zx_ratio=UNDEF;
		volume_ratio=UNDEF;
		Nmat_need=1;
	}
	else if(shape==SH_RBC) {
		/* three-parameter shape; developed by K.A.Semyanov,P.A.Tarasov,P.A.Avrorov
		 * based on work by P.W.Kuchel and E.D.Fackerell, "Parametric-equation representation
		 * of biconcave erythrocytes," Bulletin of Mathematical Biology 61, 209-220 (1999).
		 * ro^4+2S*ro^2*z^2+z^4+P*ro^2+Q*z^2+R=0, ro^2=x^2+y^2, P,Q,R,S are determined by d,h,b,c
		 * given in the command line.
		 */
		h_d=sh_pars[0];
		TestPositive(h_d,"ratio of maximum width to diameter");
		b_d=sh_pars[1];
		TestNonNegative(b_d,"ratio of minimum width to diameter");
		if (h_d<=b_d) PrintErrorHelp("given RBC is not biconcave; maximum width is in the center");
		c_d=sh_pars[2];
		TestRangeII(c_d,"relative diameter of maximum width",0,1);
		SPRINTZ(sh_form_str,
			"red blood cell; diameter(d):%%.10g, maximum and minimum width h/d=%.10g, b/d=%.10g\n"
			"       diameter of maximum width c/d=%.10g",h_d,b_d,c_d);
		// calculate shape parameters
		h2=h_d*h_d;
		b2=b_d*b_d;
		c2=c_d*c_d;
		/* P={(b/d)^2*[c^4/(h^2-b^2)-h^2]-d^2}/4; Q=(d/b)^2*(P+d^2/4)-b^2/4; R=-d^2*(P+d^2/4)/4;
		 * S=-(2P+c^2)/h^2;  here P,Q,R,S are made dimensionless dividing by respective powers of d
		 * Calculation is performed so that Q is well defined even for b=0.
		 */
		tmp1=((c2*c2/(h2-b2))-h2)/4;
		P=b2*tmp1-0.25;
		Q=tmp1-(b2/4);
		R=-b2*tmp1/4;
		S=-(2*P+c2)/h2;
		yx_ratio=1;
		zx_ratio=h_d;
		volume_ratio=UNDEF;
		Nmat_need=1;
	}
	else if (shape==SH_READ) {
		SPRINTZ(sh_form_str,"specified by file %s; size along x-axis:%%.10g",shape_fname);
		symX=symY=symZ=symR=FALSE; // input file is assumed asymmetric
		InitDipFile(shape_fname,&n_boxX,&n_boxY,&n_boxZ,&Nmat_need);
		yx_ratio=zx_ratio=UNDEF;
		volume_ratio=UNDEF;
	}
	else if (shape==SH_SPHERE) {
		STRCPYZ(sh_form_str,"sphere; diameter:%.10g");
		volume_ratio=PI_OVER_SIX;
		yx_ratio=zx_ratio=1;
		Nmat_need=1;
	}
	else if (shape==SH_SPHEREBOX) {
		coat_ratio=sh_pars[0];
		TestRangeII(coat_ratio,"sphere diameter/cube edge ratio",0,1);
		SPRINTZ(sh_form_str,
			"sphere in cube; size of cube edge(a):%%.10g, diameter of sphere d/a=%.10g",coat_ratio);
		coat_r2=0.25*coat_ratio*coat_ratio;
		yx_ratio=zx_ratio=1;
		volume_ratio=1;
		Nmat_need=2;
	}
	/* TO ADD NEW SHAPE
	 * add an option here (in the end of 'else if' sequence). Identifier ('SH_...') should be
	 * defined in const.h. The option should
	 * 1) save all the input parameters from array 'sh_pars' to local variables
	 *    (defined in the beginning of this source files)
	 * 2) test all input parameters (for that you're encouraged to use functions from param.h since
	 *    they would automatically produce informative output in case of error). If the shape can
	 *    accept different number of parameters (UNDEF was set in array shape_opt) then also test
	 *    the number of parameters.
	 * 3) if shape breaks any symmetry, corresponding variable should be set to FALSE. Do not set
	 *    any of them to TRUE, as they can be set to FALSE by some other factors.
	 *    symX, symY, symZ - symmetries of reflection over planes YZ, XZ, XY respectively.
	 *    symR - symmetry of rotation for 90 degrees over the Z axis
	 * 4) initialize the following:
	 * sh_form_str - descriptive string, should contain %g - it would be replaced by box size along
	 *               x-axis afterwards (in param.c).
	 * Either yx_ratio (preferably) or n_boxY. The former is a ratio of particle sizes along y and x
	 *          axes. Initialize n_boxY directly only if it is not proportional to boxX, like in
	 *          shape LINE above, since boxX is not initialized at this moment. If yx_ratio is not
	 *          initialized, set it explicitly to UNDEF.
	 * Analogously either zx_ratio (preferably) or n_boxZ.
	 * Nmat_need - number of different domains in this shape (void is not included)
	 * volume_ratio - ratio of particle volume to (boxX)^3. Initialize it if it can be calculated
	 *                analytically or set to UNDEF otherwise. This parameter is crucial if one wants
	 *                to initialize computational grid from '-eq_rad' and '-dpl'.
	 * n_sizeX - absolute size of the particle, defined by shape; initialize only when relevant,
	 *           e.g. for shapes such as 'axisymmetric'.
	 * all other auxiliary variables, which are used in shape generation (MakeParticle(), see
	 *   below), should be defined in the beginning of this file. If you need temporary local
	 *   variables (which are used only in this part of the code), either use 'tmp1'-'tmp3' or
	 *   define your own (with more informative names) in the beginning of this function.
	 * Also (rarely) if the shape defines dimension of the computational grid or absolute size of
	 * the particle, correct values of box_det_sh and size_det_sh in the beginning of this function.
	 */

	// initialize domain granulation
	if (sh_granul) {
		symX=symY=symZ=symR=FALSE;  // no symmetry with granules
		if (gr_mat+1>Nmat_need)
			PrintError("Specified domain number to be granulated (%d) is larger than total number "
				"of domains (%d) for the given shape (%s)",gr_mat+1,Nmat_need,shapename);
		else Nmat_need++;
		strcat(shapename,"_gran");
	}
	// check if enough refractive indices or extra
	if (Nmat<Nmat_need) {
		if (prognose) {
			if (dpl_def_used) PrintError("Given number of refractive indices (%d) is less "
				"than number of domains (%d). Since computational grid is initialized based on "
				"the default dpl, it may change depending on the actual refractive indices.",
				Nmat,Nmat_need);
		}
		else PrintError("Only %d refractive indices are given. %d are required",Nmat,Nmat_need);
	}
	else if (Nmat>Nmat_need) LogError(EC_INFO,ONE_POS,
		"More refractive indices are given (%d) than actually used (%d)",Nmat,Nmat_need);
	Nmat=Nmat_need;

	// check anisotropic refractive indices for symmetries
	if (anisotropy) for (i=0;i<Nmat;i++) symR=symR && ref_index[3*i][RE]==ref_index[3*i+1][RE]
	                                         && ref_index[3*i][IM]==ref_index[3*i+1][IM];

	if (sym_type==SYM_NO) symX=symY=symZ=symR=FALSE;
	else if (sym_type==SYM_ENF) symX=symY=symZ=symR=TRUE;

	// determine which size to use
	if (size_det_sh) {
		if (size_given_cmd) LogError(EC_INFO,ONE_POS,"Particle size specified by command line "
			"option '-%s' overrides the internal specification of the shape '%s'. The particle "
			"will be scaled accordingly.",sizename,shapename);
		else sizeX=n_sizeX;
	}
	// use analytic connection between sizeX and a_eq if available
	if (a_eq!=UNDEF && volume_ratio!=UNDEF)
		sizeX=pow(FOUR_PI_OVER_THREE/volume_ratio,ONE_THIRD)*a_eq;
	/* Initialization of boxX;
	 * if boxX is not defined by command line, it is either set by shape itself or
	 *   if sizeX is set, boxX is initialized to default
	 *   else dpl is initialized to default (if undefined) and boxX is calculated from sizeX and dpl
	 * else adjust boxX if needed.
	 */
	if (boxX==UNDEF && !box_det_sh) {
		if (sizeX==UNDEF) {
			// if a_eq is set, but sizeX was not initialized before - error
			if (a_eq!=UNDEF) PrintError("Grid size can not be automatically determined from "
				"equivalent radius and dpl for shape '%s', because its volume is not known "
				"analytically. Either use '-size' instead of '-eq_rad' or specify grid size "
				"manually by '-grid'.",shapename);
			// default value for boxX; FitBox is redundant but safer for future changes
			boxX=FitBox(DEF_GRID);
		}
		else {
			if (dpl==UNDEF) {
				/* use default dpl, but make sure that it does not produce too small grid
				 * (e.g. for nanoparticles).
				 */
				temp=(int)ceil(sizeX*dpl_def/lambda);
				boxX=FitBox(MAX(temp,MIN_AUTO_GRID));
				dpl_def_used=TRUE;
			}
			else { // if dpl is given in the command line; then believe it
				boxX=FitBox((int)ceil(sizeX*dpl/lambda));
				dpl=UNDEF; // dpl is given correct value in make_particle()
			}
		}
	}
	else {
		/* warnings are issued if specified boxX need to be adjusted,
		 * especially when '-size' is used
		 */
		if (boxX!=UNDEF) temp=boxX;
		else temp=n_boxX;
		if ((boxX=FitBox(temp))!=temp) {
			if (sizeX==UNDEF) LogError(EC_WARN,ONE_POS,"boxX has been adjusted from %i to %i. "
				"Size along X-axis in the shape description is the size of new (adjusted) "
				"computational grid.",temp,boxX);
			else LogError(EC_WARN,ONE_POS,
				"boxX has been adjusted from %i to %i. Size specified by the command line option "
				"'-size' is used for the new (adjusted) computational grid.",temp,boxX);
		}
		if (box_det_sh && n_boxX>boxX)
			PrintError("Particle (boxX=%d) does not fit into specified boxX=%d",n_boxX,boxX);
	}
	// if shape is determined by ratios, calculate proposed grid sizes along y and z axes
	if (yx_ratio!=UNDEF) n_boxY=(int)ceil(yx_ratio*boxX);
	if (zx_ratio!=UNDEF) n_boxZ=(int)ceil(zx_ratio*boxX);
	// set boxY and boxZ
	if (boxY==UNDEF) { // assumed that boxY and boxZ are either both defined or both not defined
		boxY=FitBox(n_boxY);
		boxZ=FitBox(n_boxZ);
	}
	else {
		temp=boxY;
		if ((boxY=FitBox(boxY))!=temp)
			LogError(EC_WARN,ONE_POS,"boxY has been adjusted from %i to %i",temp,boxY);
		temp=boxZ;
		if ((boxZ=FitBox(boxZ))!=temp)
			LogError(EC_WARN,ONE_POS,"boxZ has been adjusted from %i to %i",temp,boxZ);
		// this error is not duplicated in the log file since it does not yet exist
		if (n_boxY>boxY || n_boxZ>boxZ)
			PrintError("Particle (boxY,Z={%d,%d}) does not fit into specified boxY,Z={%d,%d}",
				n_boxY,n_boxZ,boxY,boxZ);
	}
	// initialize number of dipoles
	Ndip=boxX*((double)boxY)*boxZ;
	// initialize maxiter; not very realistic
	if (maxiter==UNDEF) maxiter=MIN(INT_MAX,3*Ndip);
	// some old, not really logical heuristics for Ntheta, but better than constant value
	if (nTheta==UNDEF) {
		if (Ndip<1000) nTheta=91;
		else if (Ndip<10000) nTheta=181;
		else if (Ndip<100000) nTheta=361;
		else nTheta=721;
	}
	// this limitation should be removed in the future
	if (chp_type!=CHP_NONE && (!symR || scat_grid)) LogError(EC_ERROR,ONE_POS,
		"Currently checkpoints can be used when internal fields are calculated only once,"
		"i.e. for a single incident polarization.");
	Timing_Particle = GET_TIME() - tstart;
}

//==========================================================

void MakeParticle(void)
// creates a particle; initializes all dipoles counts, dpl, gridspace
{
	int i,j,k,ns;
	size_t index,dip,nlocalRows_tmp;
	double tmp1,tmp2,tmp3;
	double xr,yr,zr,xcoat,ycoat,zcoat,r2,z2,zshift;
	double jcX,jcY,jcZ; // center for jagged
	int local_z0_unif; // should be global or semi-global
	int largerZ,smallerZ; // number of larger and smaller z in intersections with contours
	int xj,yj,zj;
	int mat;
	unsigned short us_tmp;
	TIME_TYPE tstart,tgran;
	/* TO ADD NEW SHAPE
	 * Add here all intermediate variables, which are used only inside this function. You may as
	 * well use 'tmp1'-'tmp3' variables defined above.
	 */

	tstart=GET_TIME();

	index=0;
	// assumed that box's are even
	jcX=jcY=jcZ=jagged/2.0;
	cX=(boxX-1)/2.0;
	cY=(boxY-1)/2.0;
	cZ=(boxZ-1)/2.0;
	nlocalRows_tmp=MultOverflow(3,local_Ndip,ALL_POS,"nlocalRows_tmp");
	/* allocate temporary memory; even if prognosis, since they are needed for exact estimation
	 * they will be reallocated afterwards (when nlocalRows is known).
	 */
	MALLOC_VECTOR(material_tmp,uchar,local_Ndip,ALL);
	MALLOC_VECTOR(DipoleCoord_tmp,double,nlocalRows_tmp,ALL);
	MALLOC_VECTOR(position_tmp,ushort,nlocalRows_tmp,ALL);

	for(k=local_z0;k<local_z1_coer;k++)
		for(j=0;j<boxY;j++)
			for(i=0;i<boxX;i++) {
				xj=jagged*(i/jagged)-boxX/2;
				yj=jagged*(j/jagged)-boxY/2;
				zj=jagged*(k/jagged)-boxZ/2;

				xr=(xj+jcX)/(boxX);
				yr=(yj+jcY)/(boxX);
				zr=(zj+jcZ)/(boxX);

				mat=Nmat; // corresponds to void

				if (shape==SH_AXISYMMETRIC) {
					r2=xr*xr+yr*yr;
					if (r2>=contRoSqMin && r2<=0.25) {
						largerZ=smallerZ=0;
						contCurRo=sqrt(r2);
						contCurZ=zr;
						for (ns=0;ns<contNseg;ns++)
							if (contCurRo>=contSegRoMin[ns] && contCurRo<=contSegRoMax[ns])
								CheckContourSegment(contSeg+ns) ? largerZ++ : smallerZ++;
						// check for consistency; if the code is perfect, this is not needed
						if (!IS_EVEN(largerZ+smallerZ)) LogError(EC_ERROR,ALL_POS,
							"Point (ro,z)=(%g,%g) produced weird result when checking whether it lies "
							"inside the contour. Larger than z %d intersections, smaller - %d.",
							contCurRo,contCurZ,largerZ,smallerZ);
						if (!IS_EVEN(largerZ)) mat=0;
					}
				}
				else if (shape==SH_BOX) {
					if (fabs(yr)<=haspY && fabs(zr)<=haspZ) mat=0;
				}
				else if (shape==SH_CAPSULE) {
					r2=xr*xr+yr*yr;
					if (r2<=0.25) {
						tmp1=fabs(zr)-hdratio;
						if (tmp1<=0 || tmp1*tmp1+r2<=0.25) mat=0;
					}
				}
				else if (shape==SH_COATED) {
					if (xr*xr+yr*yr+zr*zr<=0.25) { // first test to skip some dipoles immediately)
						xcoat=xr-coat_x;
						ycoat=yr-coat_y;
						zcoat=zr-coat_z;
						if (xcoat*xcoat+ycoat*ycoat+zcoat*zcoat<=coat_r2) mat=1;
						else mat=0;
					}
				}
				else if (shape==SH_CYLINDER) {
					if(xr*xr+yr*yr<=0.25 && fabs(zr)<=hdratio) mat=0;
				}
				else if (shape==SH_EGG) {
					r2=xr*xr+yr*yr;
					zshift=zr+egz0;
					z2=zshift*zshift;
					if (r2+egeps*z2+egnu*zshift*sqrt(r2+z2)<=ad2) mat=0;
				}
				else if (shape==SH_ELLIPSOID) {
					if (xr*xr+yr*yr*invsqY+zr*zr*invsqZ<=0.25) mat=0;
				}
				else if (shape==SH_LINE) {
					if (yj==0 && zj==0) mat=0;
				}
				else if (shape==SH_RBC) {
					r2=xr*xr+yr*yr;
					z2=zr*zr;
					if (r2*r2+2*S*r2*z2+z2*z2+P*r2+Q*z2+R<=0) mat=0;
				}
				else if (shape==SH_SPHERE) {
					if (xr*xr+yr*yr+zr*zr<=0.25) mat=0;
				}
				else if (shape==SH_SPHEREBOX) {
					if (xr*xr+yr*yr+zr*zr<=coat_r2) mat=1;
					else if (fabs(yr)<=0.5 && fabs(zr)<=0.5) mat=0;
				}
				/* TO ADD NEW SHAPE
				 * add an option here (in the end of 'else if' sequence). Identifier ('SH_...')
				 * should be defined in const.h. This option should set 'mat' - index of domain for
				 * a point, specified by {xr,yr,zr} - coordinates divided by grid size along X (xr
				 * from -0.5 to 0.5, others - depending on aspect ratios). C array indexing used:
				 * mat=0 - first domain, etc. If point corresponds to void, do not set 'mat'. If you
				 * need temporary local variables (which are used only in this part of the code),
				 * either use 'tmp1'-'tmp3' or define your own (with more informative names) in the
				 * beginning of this function.
				 */

				position_tmp[3*index]=(unsigned short)i;
				position_tmp[3*index+1]=(unsigned short)j;
				position_tmp[3*index+2]=(unsigned short)k;
				// afterwards multiplied by gridspace
				DipoleCoord_tmp[3*index]=i-cX;
				DipoleCoord_tmp[3*index+1]=j-cY;
				DipoleCoord_tmp[3*index+2]=k-cZ;
				material_tmp[index]=(unsigned char)mat;
				index++;
			} // End box loop
	if (shape==SH_READ) ReadDipFile(shape_fname);
	// initialization of mat_count and dipoles counts
	for(i=0;i<=Nmat;i++) mat_count[i]=0;
	for(dip=0;dip<local_Ndip;dip++) mat_count[material_tmp[dip]]++;
	local_nvoid_Ndip=local_Ndip-mat_count[Nmat];
	MyInnerProduct(mat_count,double_type,Nmat+1,NULL);
	if ((nvoid_Ndip=Ndip-mat_count[Nmat])==0)
		LogError(EC_ERROR,ONE_POS,"All dipoles of the scatterer are void");
	nlocalRows=3*local_nvoid_Ndip;
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
		LogError(EC_ERROR,ONE_POS,"Too small dpl for FCD formulation, should be at least 2");
	gridspace=lambda/dpl;
	// initialize equivalent size parameter and cross section
	kd = TWO_PI/dpl;
	/* from this moment on a_eq and all derived quantities are based on the real a_eq, which can
	 * in several cases be slightly different from the one given by '-eq_rad' option.
	 */
	a_eq = pow(THREE_OVER_FOUR_PI*nvoid_Ndip,ONE_THIRD)*gridspace;
	ka_eq = WaveNum*a_eq;
	inv_G = 1/(PI*a_eq*a_eq);
	// granulate one domain, if needed
	if (sh_granul) {
		tgran=GET_TIME();
		Timing_Granul_comm=0;
		// calculate number of granules
		if (mat_count[gr_mat]==0) LogError(EC_ERROR,ONE_POS,
			"Domain to be granulated does not contain any dipoles");
		tmp1=gridspace/gr_d;
		tmp2=mat_count[gr_mat]*gr_vf*SIX_OVER_PI;
		tmp3=tmp2*tmp1*tmp1*tmp1;
		CheckOverflow(tmp3,ONE_POS,"Make_Particle()");
		gr_N=(size_t)ceil(tmp3);
		// correct granules diameter to get exact volume fraction (if volume correction is used)
		if (volcor_used) gr_d=gridspace*pow(tmp2/gr_N,ONE_THIRD);
		// actually place granules
		mat_count[Nmat-1]=PlaceGranules();
		// calculate exact volume fraction
		gr_vf_real=mat_count[Nmat-1]/mat_count[gr_mat];
		mat_count[gr_mat]-=mat_count[Nmat-1];
		Timing_Granul=GET_TIME()-tgran;
	}
	/* allocate main particle arrays, using precise nlocalRows even when prognosis is used to enable
	 * save_geom afterwards.
	 */
	MALLOC_VECTOR(material,uchar,local_nvoid_Ndip,ALL);
	MALLOC_VECTOR(DipoleCoord,double,nlocalRows,ALL);
	MALLOC_VECTOR(position,ushort,nlocalRows,ALL);

	memory+=(3*(sizeof(short int)+sizeof(double))+sizeof(char))*local_nvoid_Ndip;
	// copy nontrivial part of arrays
	index=0;
	for (dip=0;dip<local_Ndip;dip++) if (material_tmp[dip]<Nmat) {
		material[index]=material_tmp[dip];
		// DipoleCoord=gridspace*DipoleCoord_tmp
		MultScal(gridspace,DipoleCoord_tmp+3*dip,DipoleCoord+3*index);
		memcpy(position+3*index,position_tmp+3*dip,3*sizeof(short int));
		index++;
	}
	// free temporary memory
	Free_general(material_tmp);
	Free_general(DipoleCoord_tmp);
	Free_general(position_tmp);
	if (shape==SH_AXISYMMETRIC) {
		for (ns=0;ns<contNseg;ns++) FreeContourSegment(contSeg+ns);
		Free_general(contSegRoMin);
		Free_general(contSegRoMax);
	}

	// save geometry
	if (save_geom) SaveGeometry();

	/* adjust z-axis of position vector, to speed-up matrix-vector multiplication a little bit;
	 * after this point 'position(z)' is taken relative to the local_z0.
	 */
	if (local_z0!=0) {
		us_tmp=(unsigned short)local_z0;
		for (dip=2;dip<3*local_nvoid_Ndip;dip+=3) position[dip]-=us_tmp;
	}
	local_Nz_unif=position[3*local_nvoid_Ndip-1]+1;
	local_z0_unif=local_z0; // TODO: should be changed afterwards
	box_origin_unif[0]=-gridspace*cX;
	box_origin_unif[1]=-gridspace*cY;
	box_origin_unif[2]=gridspace*(local_z0_unif-cZ);

	Timing_Particle += GET_TIME() - tstart;
}
