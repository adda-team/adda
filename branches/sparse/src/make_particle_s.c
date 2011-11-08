/* FILE : make_particle_s.c
 * Descr: this module initializes the dipole set, either using predefined shapes or reading from a
 *        file; a specialized version for the sparse mode ADDA
 *
 * Copyright (C) 2006-2011 ADDA contributors
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

/*
	A version of the functions in make_particle.h designed for the data model of
	the sparse mode of ADDA. Currently supports only reading the input from a file.
*/

//this should not be included in non-sparse builds, but making sure with #ifdef
#ifdef ADDA_SPARSE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "const.h"
#include "param.h"
#include "vars.h"
#include "memory.h"
#include "timing.h"
#include "cmplx.h"
#include "comm.h"
#include "debug.h"

bool volcor_used;                // volume correction was actually employed
char sh_form_str[MAX_PARAGRAPH]; // string for log file with shape parameters
double mat_count[MAX_NMAT+1];

static int minX,minY,minZ;      // minimum values of dipole positions in dipole file
static double cX,cY,cZ;         // center for DipoleCoord
double gridspace;

// defined and initialized in timing.c
extern TIME_TYPE Timing_Particle;

static double dpl_def;          // default value of dpl

extern const int jagged;
extern const char shape_fname[];
extern double sizeX;                    // size of particle along x-axis
extern double dpl;                      // number of dipoles per lambda (wavelength)
extern double a_eq;
extern const double lambda;

static const char geom_format[]="%d %d %d\n";              // format of the geom file
static const char geom_format_ext[]="%d %d %d %d\n";       // extended format of the geom file
/* DDSCAT shape formats; several format are used, since first variable is unpredictable and last two
 * are not actually used (only to produce warnings
 */
static const char ddscat_format_read1[]="%*f %d %d %d %d %d %d\n";
static const char ddscat_format_read2[]="%*f %d %d %d %d";
static const char ddscat_format_write[]="%zd %d %d %d %d %d %d\n";
static FILE * restrict dipfile; // handle of dipole file
static enum shform read_format; // format of dipole file, which is read
static char linebuf[BUF_LINE];  // buffer for reading lines from dipole file

static void InitDipFile(const char * restrict fname, int *Nm, int *Nd);
static void ReadDipFile(const char * restrict fname);

void InitShape(void) 
{		
	double tmp1,tmp2;
	time_t tstart=GET_TIME();	
	int n_nonvoid;
	
	symX=symY=symZ=symR=false;		
	
	InitDipFile(shape_fname, &Nmat, &n_nonvoid);
	nvoid_Ndip = n_nonvoid;
	
	tmp2=0;
	for (int i=0;i<Ncomp*Nmat;i++) {
		tmp1=cAbs2(ref_index[i]);
		if (tmp2<tmp1) tmp2=tmp1;
	}
	dpl_def=10*sqrt(tmp2);
	
	int Ndip = boxX*((double)boxY)*boxZ; 
	if (maxiter==UNDEF) maxiter=MIN(INT_MAX,3*Ndip);
	if (nTheta==UNDEF) {
		if (Ndip<1000) nTheta=91;
		else if (Ndip<10000) nTheta=181;
		else if (Ndip<100000) nTheta=361;
		else nTheta=721;
	}
	
	if (IFROOT) SnprintfErr(ONE_POS,sh_form_str,MAX_PARAGRAPH,
			"specified by file %s; size along x-axis:%s",shape_fname,GFORM);
	volcor_used=false;	

	// this limitation should be removed in the future
	if (chp_type!=CHP_NONE && (!symR || scat_grid)) LogError(ONE_POS,"Currently checkpoints can be "
		"used when internal fields are calculated only once,i.e. for a single incident "
		"polarization.");
	
	Timing_Particle = GET_TIME() - tstart;
}

void MakeParticle(void)
{	
	const int jagged_cu = jagged*jagged*jagged;
	size_t index;
	
	time_t tstart=GET_TIME();
	
	nvoid_Ndip *= jagged_cu;
	MALLOC_VECTOR(material_full,uchar,nvoid_Ndip,ALL);	
	MALLOC_VECTOR(position_full,int,nvoid_Ndip*3,ALL);
	memory += nvoid_Ndip*sizeof(char) + nvoid_Ndip*3*sizeof(int);
	
	local_d0 = local_f0*nvoid_Ndip;
	local_d1 = local_f1*nvoid_Ndip;
	local_nvoid_Ndip=local_d1-local_d0;
	nlocalRows = 3*local_nvoid_Ndip;
		
	ReadDipFile(shape_fname);
	
	// initialize dpl and gridspace
	if (sizeX==UNDEF) {
		if (a_eq!=UNDEF) dpl=lambda*pow(nvoid_Ndip*THREE_OVER_FOUR_PI,ONE_THIRD)/a_eq;
		else if (dpl==UNDEF) dpl=dpl_def; // default value of dpl
		// sizeX is determined to give correct volume
		sizeX=lambda*boxX/dpl;
	}
	else {
		// dpl is determined to give correct volume
		dpl=lambda*boxX/sizeX;
	}
	// Check consistency for FCD
	if ((IntRelation==G_FCD || PolRelation==POL_FCD) && dpl<=2)
		LogError(ONE_POS,"Too small dpl for FCD formulation, should be at least 2");
	// initialize gridspace and dipvol
	gridspace=lambda/dpl;
	dipvol=gridspace*gridspace*gridspace;
	// initialize equivalent size parameter and cross section
	kd = TWO_PI/dpl;
	/* from this moment on a_eq and all derived quantities are based on the real a_eq, which can
	 * in several cases be slightly different from the one given by '-eq_rad' option.
	 */
	a_eq = pow(THREE_OVER_FOUR_PI*nvoid_Ndip,ONE_THIRD)*gridspace;
	ka_eq = WaveNum*a_eq;
	inv_G = 1/(PI*a_eq*a_eq);

	index=nvoid_Ndip/jagged_cu;
	if (jagged>1) {
		//multiply dipoles according to jagged
		for (unsigned long i=0; i<nvoid_Ndip/jagged_cu; i++) {
			position[3*i] *= jagged;
			position[3*i+1] *= jagged;
			position[3*i+2] *= jagged;
			for (int dz=0; dz<jagged; dz++)
				for (int dy=0; dy<jagged; dy++)
					for (int dx=0; dx<jagged; dx++) {
						position[3*index] = position[3*i]+dx;
						position[3*index+1] = position[3*i+1]+dy;
						position[3*index+2] = position[3*i+2]+dz;
						index++;
					}
		} 				
	}	
		
	MALLOC_VECTOR(material,uchar,local_nvoid_Ndip,ALL);
	MALLOC_VECTOR(DipoleCoord,double,nlocalRows,ALL);
	MALLOC_VECTOR(position,int,nlocalRows,ALL);
	memory+=(3*(sizeof(short int)+sizeof(double))+sizeof(char))*local_nvoid_Ndip;
	memcpy(material, &(material_full[local_d0]), local_nvoid_Ndip*sizeof(*material_full));
	memcpy(position, &(position_full[3*local_d0]), nlocalRows*sizeof(*position_full));	
	
	for (index=0; index<local_nvoid_Ndip; index++) {
		DipoleCoord[3*index] = (position[3*index] - cX) * gridspace;
		DipoleCoord[3*index+1] = (position[3*index+1] - cY) * gridspace;
		DipoleCoord[3*index+2] = (position[3*index+2] - cZ) * gridspace;
	}
	
	/*
	D("cX: %f, cY: %f, cZ: %f\n", cX, cY, cZ);
	printf("%u:position\n", ringid);
	for (unsigned int i=0; i<3*local_nvoid_Ndip; i++)
		printf("%u-%u:%d ", ringid, i, position[i]);
	printf("\n");
	printf("%u:DipoleCoord\n", ringid);
	for (unsigned int i=0; i<3*local_nvoid_Ndip; i++)
		printf("%u-%u:%.4f ", ringid, i, DipoleCoord[i]);
	printf("\n");
	*/
	
	box_origin_unif[0]=-gridspace*cX;
	box_origin_unif[1]=-gridspace*cY;
	box_origin_unif[2]=-gridspace*cZ;

	D("local_d0: %u, local_d1: %u, local_f0: %f, local_f1: %f,\n nlocalRows: %d, local_nvoid_Ndip: %d", 
	  local_d0, local_d1, local_f0, local_f1, nlocalRows, local_nvoid_Ndip);

	InitArraySync();	
		
	/*
	for (unsigned int i=0; i<nlocalRows; i++)
		printf("%d ", position[i]);
	printf("\n");
	*/
	SyncPosition(position_full);
	SyncMaterial(material_full);

	Timing_Particle += GET_TIME() - tstart;	
	 
}

INLINE void SkipFullLine(FILE * restrict file)
// skips full line in the file, starting from current position; it uses predefined buffer 'linebuf'
{
	do fgets(linebuf,BUF_LINE,file); while (strchr(linebuf,'\n')==NULL && !feof(file));
}

//===========================================================

INLINE char *FgetsError(FILE * restrict file,const char * restrict fname,size_t *line,ERR_LOC_DECL)
/* calls fgets, checks for errors and increments line number; s_fname and s_line are source fname
 * and line number to be shown in error message; result is stored in predefined buffer 'linebuf'.
 */
{
	char *res;

	res=fgets(linebuf,BUF_LINE,file);
	if (res!=NULL) {
		(*line)++;
		if (strchr(linebuf,'\n')==NULL && !feof(file)) LogError(ERR_LOC_CALL,
			"Buffer overflow while scanning lines in file '%s' (size of line %zu > %d)",
			fname,*line,BUF_LINE-1);
	}
	return res;
}

//===========================================================

INLINE void SkipNLines(FILE * restrict file,size_t n)
// skips n lines from the file starting from current position in a file
{
	while (n>0) {
		SkipFullLine(file);
		n--;
	}
}

//===========================================================

static size_t SkipComments(FILE * restrict file)
/* skips comments (#...), all lines, starting from current position in a file.
 * returns number of lines skipped
 */
{
	int ch;
	size_t lines=0;

	while ((ch=fgetc(file))=='#') {
		SkipFullLine(file);
		lines++;
	}
	if (ch!=EOF) ungetc(ch,file);

	return lines;
}

//===========================================================
#define DDSCAT_HL 6 // number of header lines in DDSCAT format

static void InitDipFile(const char * restrict fname, int *Nm, int *Nd)
/* read dipole file first to determine box sizes, Nmat and dipole number; input is not checked for very large
 * numbers (integer overflows) to increase speed; this function opens file for reading, the file is
 * closed in ReadDipFile.
 */
{
	int x,y,z,mat,scanned,mustbe;
	size_t line,skiplines;
	bool anis_warned;
	int t2,t3; // dumb variables
	int maxX,maxY,maxZ,maxN;
	char formtext[MAX_LINE];
	
	*Nd=0;

	dipfile=FOpenErr(fname,"r",ALL_POS);
	/* test for DDSCAT format; in not-DDSCAT format, the line scanned below may be a long comment;
	 * therefore we first skip all comments
	 */
	line=SkipComments(dipfile);
	if (line<=DDSCAT_HL
		&& (SkipNLines(dipfile,DDSCAT_HL-line),FgetsError(dipfile,fname,&line,ONE_POS)!=NULL)
		&& sscanf(linebuf,ddscat_format_read1,&x,&y,&z,&mat,&t2,&t3)==6) {

		read_format=SF_DDSCAT;
		strcpy(formtext,"DDSCAT format (FRMFIL)");
		mustbe=6;
		line=DDSCAT_HL;
		fseek(dipfile,0,SEEK_SET);
		SkipNLines(dipfile,line);
	}
	else { // test for ADDA text formats
		fseek(dipfile,0,SEEK_SET);
		line=SkipComments(dipfile);
		/* scanf and analyze Nmat; if there is blank line between comments and Nmat, it fails later;
		 * the value of Nmat obtained here is not actually relevant, the main factor is maximum
		 * domain number among all dipoles.
		 */
		scanned=fscanf(dipfile,"Nmat=%d\n",Nm);
		if (scanned==EOF) LogError(ONE_POS,"No dipole positions are found in %s",fname);
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
	D("%s was detected",formtext);
	// scan main part of the file
	skiplines=line;
	maxX=maxY=maxZ=INT_MIN;
	minX=minY=minZ=INT_MAX;
	maxN=1;
	anis_warned=false;
	// reading is performed in lines
	while(FgetsError(dipfile,fname,&line,ONE_POS)!=NULL) {
		// scan numbers in a line
		if (read_format==SF_TEXT) scanned=sscanf(linebuf,geom_format,&x,&y,&z);
		else if (read_format==SF_TEXT_EXT) scanned=sscanf(linebuf,geom_format_ext,&x,&y,&z,&mat);
		// for ddscat format, only first material is used, other two are ignored
		else { // read_format==SF_DDSCAT
			scanned=sscanf(linebuf,ddscat_format_read1,&x,&y,&z,&mat,&t2,&t3);
			if (!anis_warned && (t2!=mat || t3!=mat)) {
				LogWarning(EC_WARN,ONE_POS,"Anisotropic dipoles are detected in file %s (first on "
					"line %zu). ADDA ignores this anisotropy, using only the identifier of "
					"x-component of refractive index as domain number",fname,line);
				anis_warned=true;
			}
		}
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			if (scanned!=mustbe) // this in most cases indicates wrong format
				LogError(ONE_POS,"%s was detected, but error occurred during scanning of line %zu "
					"from dipole file %s",formtext,line,fname);
			if (read_format!=SF_TEXT) {
				if (mat<=0) LogError(ONE_POS,"%s was detected, but nonpositive material number "
					"(%d) encountered during scanning of line %zu from dipole file %s",
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
			// increase dipole count
			(*Nd)++;
		}
	}
	if (read_format==SF_TEXT_EXT) {
		if (*Nm!=maxN) LogWarning(EC_WARN,ONE_POS,"Nmat (%d), as given in %s, is not equal to the "
				"maximum domain number (%d) among all specified dipoles; hence the former is "
				"ignored",*Nm,fname,maxN);
	}
	
	boxX = (maxX-minX+1)*jagged;
	boxY = (maxY-minY+1)*jagged;
	boxZ = (maxZ-minZ+1)*jagged;	
	boxXY = boxX*boxY;	
	cX = (boxX-jagged)/2.0;
	cY = (boxY-jagged)/2.0;
	cZ = (boxZ-jagged)/2.0;
	
	// not optimal way, but works more robustly when non-system EOL is used in data file
	fseek(dipfile,0,SEEK_SET);
	SkipNLines(dipfile,skiplines);
}
#undef DDSCAT_HL

//===========================================================

static void ReadDipFile(const char * restrict fname) 
/* read dipole file; no consistency checks are made since they are made in InitDipFile.
 * the file is opened in InitDipFile; this function only closes the file.
 */
{
	int x0,y0,z0,mat,scanned;
	unsigned int index=0;
	size_t boxX_l;
	const int jagged_cu = jagged*jagged*jagged;

	// to remove possible overflows
	boxX_l=(size_t)boxX;	

	mat=1;
	while(fgets(linebuf,BUF_LINE,dipfile)!=NULL) {
		// scan numbers in a line
		if (read_format==SF_TEXT) scanned=sscanf(linebuf,geom_format,&x0,&y0,&z0);
		else if (read_format==SF_TEXT_EXT) scanned=sscanf(linebuf,geom_format_ext,&x0,&y0,&z0,&mat);
		else scanned=sscanf(linebuf,ddscat_format_read2,&x0,&y0,&z0,&mat); // read_format==SF_DDSCAT
		// if sscanf returns EOF, that is a blank line -> just skip
		if (scanned!=EOF) {
			if ((jagged_cu*index >= local_d0) && (jagged_cu*index < local_d1)) {
				x0-=minX;
				y0-=minY;
				z0-=minZ;
				material_full[index]=(unsigned char)(mat-1);
				position_full[3*index]=x0;
				position_full[3*index+1]=y0;
				position_full[3*index+2]=z0;
			}
			index++;		
		}
	}
	FCloseErr(dipfile,fname,ALL_POS);
}

#endif //ADDA_SPARSE
