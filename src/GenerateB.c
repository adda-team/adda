/* File: GenerateB.c
 * $Date::                            $
 * Descr: generate a incident beam
 *
 *        Lminus beam is based on: G. Gouesbet, B. Maheu, G. Grehan, "Light scattering from a sphere arbitrary located
 *        in a Gaussian beam, using a Bromwhich formulation", J.Opt.Soc.Am.A 5,1427-1443 (1988).
 *        Eq.(22) - complex conjugate.
 *
 *        Davis beam is based on: L. W. Davis, "Theory of electromagnetic beams," Phys.Rev.A 19, 1177-1179 (1979).
 *        Eqs.(15a),(15b) - complex conjugate; in (15a) "Q" changed to "Q^2" (typo).
 *
 *        Barton beam is based on: J. P. Barton and D. R. Alexander, "Fifth-order corrected electromagnetic-field
 *        components for a fundamental Gaussian-beam," J.Appl.Phys. 66,2800-2802 (1989).
 *        Eqs.(25)-(28) - complex conjugate.
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
#include "io.h"
#include "param.h"
#include "vars.h"
// system headers
#include <stdio.h>
#include <string.h>

// SEMI-GLOBAL VARIABLES

// defined and initialized in param.c
extern const int beam_Npars;
extern const double beam_pars[];
extern const char *beam_fnameY;
extern const char *beam_fnameX;
extern opt_index opt_beam;

// used in crosssec.c
double beam_center_0[3]; // position of the beam center in laboratory reference frame
// used in param.c
char beam_descr[MAX_MESSAGE2]; // string for log file with beam parameters

// LOCAL VARIABLES
static double s,s2;            // beam confinement factor and its square
static double scale_x,scale_z; // multipliers for scaling coordinates
/* TO ADD NEW BEAM
 * Add here all internal variables (beam parameters), which you initialize in InitBeam() and use in GenerateB()
 * afterwards. If you need local, intermediate variables, put them into the beginning of the corresponding function.
 * Add descriptive comments, use 'static'.
 */

//======================================================================================================================

void InitBeam(void)
// initialize beam; produce description string
{
	double w0; // beam width
	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function.
	 */

	// initialization of global option index for error messages
	opt=opt_beam;
	// beam initialization
	switch (beamtype) {
		case B_PLANE:
			if (IFROOT) strcpy(beam_descr,"plane wave");
			beam_asym=false;
			return;
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			// initialize parameters
			w0=beam_pars[0];
			TestPositive(w0,"beam width");
			beam_asym=(beam_Npars==4 && (beam_pars[1]!=0 || beam_pars[2]!=0 || beam_pars[3]!=0));
			if (beam_asym) {
				memcpy(beam_center_0,beam_pars+1,3*sizeof(double));
				// if necessary break the symmetry of the problem
				if (beam_center_0[0]!=0) symX=symR=false;
				if (beam_center_0[1]!=0) symY=symR=false;
				if (beam_center_0[2]!=0) symZ=false;
			}
			else beam_center[0]=beam_center[1]=beam_center[2]=0;
			s=1/(WaveNum*w0);
			s2=s*s;
			scale_x=1/w0;
			scale_z=s*scale_x; // 1/(k*w0^2)
			// beam info
			if (IFROOT) {
				strcpy(beam_descr,"Gaussian beam (");
				switch (beamtype) {
					case B_LMINUS:
						strcat(beam_descr,"L- approximation)\n");
						break;
					case B_DAVIS3:
						strcat(beam_descr,"3rd order approximation, by Davis)\n");
						break;
					case B_BARTON5:
						strcat(beam_descr,"5th order approximation, by Barton)\n");
						break;
					default: break;
				}
				sprintf(beam_descr+strlen(beam_descr),"\tWidth="GFORMDEF" (confinement factor s="GFORMDEF")\n",w0,s);
				if (beam_asym) sprintf(beam_descr+strlen(beam_descr),"\tCenter position: "GFORMDEF3V,beam_center_0[0],
					beam_center_0[1],beam_center_0[2]);
				else strcat(beam_descr,"\tCenter is in the origin");
			}
			return;
		case B_READ:
			// the safest is to assume cancellation of all symmetries
			symX=symY=symZ=symR=false;
			if (IFROOT) {
				if (beam_Npars==1) sprintf(beam_descr,"specified by file '%s'",beam_fnameY);
				else sprintf(beam_descr,"specified by files '%s' and '%s'",beam_fnameY,beam_fnameX);
			}
			// we do not define beam_asym here, because beam_center is not defined anyway
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. The case should
	 * 1) save all the input parameters from array 'beam_pars' to local variables (defined in the beginning of this
	 *    source files)
	 * 2) test all input parameters (for that you're encouraged to use functions from param.h since they would
	 *    automatically produce informative output in case of error).
	 * 3) if shape breaks any symmetry, corresponding variable should be set to false. Do not set any of them to true,
	 *    as they can be set to false by other factors.
	 *    symX, symY, symZ - symmetries of reflection over planes YZ, XZ, XY respectively.
	 *    symR - symmetry of rotation for 90 degrees over the Z axis
	 * 4) initialize the following:
	 *    beam_descr - descriptive string, which will appear in log file.
	 *    beam_asym - whether beam center does not coincide with the reference frame origin. If it is set to true, then
	 *                set also beam_center_0 - 3D radius-vector of beam center in the laboratory reference frame (it
	 *                will be then automatically transformed to particle reference frame, if required).
	 * All other auxiliary variables, which are used in beam generation (GenerateB(), see below), should be defined in
	 * the beginning of this file. If you need temporary local variables (which are used only in this part of the code),
	 * define them in the beginning of this function.
	 */
}
//======================================================================================================================

static void ReadBeam(const enum incpol which,doublecomplex *restrict b)
// Reads beam from file, a separate file is used for each incident polarization
{
	char linebuf[BUF_LINE];
	const char *fname;
	if (which==INCPOL_Y) fname=beam_fnameY;
	else { // which==INCPOL_X
		if (beam_Npars==1)
			LogError(ONE_POS,"Only one beam file is specified, while X polarization need to be calculated");
		fname=beam_fnameX;
	}
	TIME_TYPE tstart=GET_TIME();
	FILE *file=FOpenErr(fname,"r",ALL_POS);
	// the same format as used for saving the beam by StoreFields(...) in make_particle.c
	const char beam_format[]="%*f %*f %*f %*f %lf %lf %lf %lf %lf %lf";
	const int mustbe=6;
	int scanned;
	size_t i,j;
	// skips first line with headers and any comments, if present
	size_t line=SkipNLines(file,1);
	line+=SkipComments(file);

	i=j=0;
	while(FGetsError(file,fname,&line,linebuf,BUF_LINE,ONE_POS)!=NULL) {
		// scan numbers in a line
		scanned=sscanf(linebuf,beam_format,&(b[j][RE]),&(b[j][IM]),&(b[j+1][RE]),&(b[j+1][IM]),&(b[j+2][RE]),
			&(b[j+2][IM]));
		if (scanned!=EOF) { // if sscanf returns EOF, that is a blank line -> just skip
			if (scanned!=mustbe) // this in most cases indicates wrong format
				LogError(ALL_POS,"Error occurred during scanning of line %zu from beam file %s",line,fname);
			if (i==nvoid_Ndip) LogError(ALL_POS,"Beam file %s contains more data rows than number of dipoles (%zu) in "
				"the particle",fname,nvoid_Ndip);
			if (i>=local_nvoid_d0) j+=3;
			i++;
			/* all processors stop reading file as soon as possible, but the last processor reads one more line to test
			 * (above) for extra strings in the file
			 */
			if(i==local_nvoid_d1 && i!=nvoid_Ndip) break;
		}
		Timing_FileIO+=GET_TIME()-tstart;
	}
}

//======================================================================================================================

void GenerateB (const enum incpol which,   // x - or y polarized incident light
                doublecomplex *restrict b) // the b vector for the incident field
// generates incident beam at every dipole
{
	size_t i,j;
	doublecomplex psi0,Q,Q2;
	doublecomplex v1[3],v2[3],v3[3];
	double ro2,ro4;
	double x,y,z,x2_s,xy_s;
	doublecomplex t1,t2,t3,t4,t5,t6,t7,t8,t0,ctemp;
	const double *ex; // coordinate axis of the beam reference frame
	double ey[3];
	double r1[3];
	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function. You may as well use 't1'-'t8'
	 * variables defined above.
	 */

	// set reference frame of the beam; ez=prop, ex - incident polarization
	if (which==INCPOL_Y) {
		ex=incPolY;
		vMultScal(-1,incPolX,ey);
	}
	else { // which==INCPOL_X
		ex=incPolX;
		memcpy(ey,incPolY,3*sizeof(double));
	}

	switch (beamtype) {
		case B_PLANE: // plane is separate to be fast
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				imExp(WaveNum*DotProd(DipoleCoord+j,prop),ctemp); // ctemp=exp(ik*r.a)
				cScalMultRVec(ex,ctemp,b+j); // b[i]=ctemp*ex
			}
			return;
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				// set relative coordinates (in beam's coordinate system)
				LinComb(DipoleCoord+j,beam_center,1,-1,r1);
				x=DotProd(r1,ex)*scale_x;
				y=DotProd(r1,ey)*scale_x;
				z=DotProd(r1,prop)*scale_z;
				ro2=x*x+y*y;
				// calculate Q=1/(2z-i)
				Q[IM]=1/(1+4*z*z);
				Q[RE]=2*z*Q[IM];
				// calculate psi0=-iQexp(iQro^2)
				cMult_i2(Q,t1);
				cMultReal(ro2,t1,t2);
				cExpSelf(t2);
				cMult(t1,t2,psi0);
				cInvSign(psi0);
				// ctemp=exp(ik*z0)*psi0, z0 - non-scaled coordinate (z/scale_z)
				imExp(WaveNum*z/scale_z,ctemp);
				cMultSelf(ctemp,psi0);
				// the following logic (if-else-if...) is hard to replace by a simple switch
				if (beamtype==B_LMINUS) {
					cScalMultRVec(ex,ctemp,b+j); // b[i]=ctemp*ex
				}
				else {
					x2_s=x*x/ro2;
					cSquare(Q,Q2);
					ro4=ro2*ro2;
					// some combinations that are used more than once
					cMultReal(s2*ro2,Q2,t4); // t4=(s*ro*Q)^2
					cMultReal(ro2,Q,t5);
					cMult_i(t5);             // t5=i*Q*ro^2
					cMultReal(ro4,Q2,t6);    // t6=ro^4*Q^2
					cMultReal(x*s,Q,t7);     // t7=x*s*Q
					if (beamtype==B_DAVIS3) {
						// t1=1+s^2(-4Q^2*x^2-iQ^3*ro^4)=1-t4(4x2_s+t5)
						cEqual(t5,t1);
						t1[RE]+=4*x2_s;
						cMultSelf(t1,t4);
						cMultReal(-1,t1,t1);
						t1[RE]+=1;
						// t2=0
						t2[RE]=t2[IM]=0;
						// t3=-s(2Qx)+s^3(8Q^3*ro^2*x+2iQ^4*ro^4*x-4iQ^2x)=2t7[-1+iQ*s2*(-4t5+t6-2)]
						cMultReal(-4,t5,t3);
						cAdd(t3,t6,t3);
						t3[RE]-=2;
						cMultReal(s2,t3,t3);
						cMultSelf(t3,Q);
						cMult_i(t3);
						t3[RE]-=1;
						cMultSelf(t3,t7);
						cMultReal(2,t3,t3);
					}
					else if (beamtype==B_BARTON5) {
						xy_s=x*y/ro2;
						cMultReal(2,t5,t8);
						t8[RE]+=8; // t8=8+2i*Q*ro^2
						/* t1 = 1 + s^2(-ro^2*Q^2-i*ro^4*Q^3-2Q^2*x^2)
						 *    + s^4[2ro^4*Q^4+3iro^6*Q^5-0.5ro^8*Q^6+x^2(8ro^2*Q^4+2iro^4*Q^5)]
						 *    = 1 + t4*{-1-2xs2-t5+t4*[2+3t5-0.5t6+x2_s*t8]}
						 */
						cMultReal(x2_s,t8,t1);
						cMultReal(-0.5,t6,t0);
						cAdd(t1,t0,t1);
						cMultReal(3,t5,t0);
						cAdd(t1,t0,t1);
						t1[RE]+=2;
						cMultSelf(t1,t4);
						cSubtr(t1,t5,t1);
						t1[RE]-=1+2*x2_s;
						cMultSelf(t1,t4);
						t1[RE]+=1;
						// t2=s^2(-2Q^2*xy)+s^4[xy(8ro^2*Q^4+2iro^4*Q^5)]=xy_s*t4(-2+t4*t8)
						cMult(t4,t8,t2);
						t2[RE]-=2;
						cMultSelf(t2,t4);
						cMultReal(xy_s,t2,t2);
						/* t3 = s(-2Qx) + s^3[(6ro^2*Q^3+2iro^4*Q^4)x] + s^5[(-20ro^4*Q^5-10iro^6*Q^6+ro^8*Q^7)x]
						 *    = t7{-2+t4[6+2t5+t4(-20-10t5+t6)]}
						 */
						cMultReal(-10,t5,t3);
						cAdd(t3,t6,t3);
						t3[RE]-=20;
						cMultSelf(t3,t4);
						cMultReal(2,t5,t0);
						cAdd(t3,t0,t3);
						t3[RE]+=6;
						cMultSelf(t3,t4);
						t3[RE]-=2;
						cMultSelf(t3,t7);
					} else LogError(ONE_POS,"Inconsistency in beam definition"); // to remove warnings
					// b[i]=ctemp(ex*t1+ey*t2+ez*t3)
					cScalMultRVec(ex,t1,v1);
					cScalMultRVec(ey,t2,v2);
					cScalMultRVec(prop,t3,v3);
					cvAdd2Self(v1,v2,v3);
					cvMultScal_cmplx(ctemp,v1,b+j);
				}
			}
			return;
		case B_READ:
			ReadBeam(which,b);
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. This case should set
	 * complex vector 'b', describing the incident field in the particle reference frame. It is set inside the cycle for
	 * each dipole of the particle and is calculated using
	 * 1) 'DipoleCoord' – array of dipole coordinates;
	 * 2) 'prop' – propagation direction of the incident field;
	 * 3) 'ex' – direction of incident polarization;
	 * 4) 'ey' – complementary unity vector of polarization (orthogonal to both 'prop' and 'ex');
	 * 5) 'beam_center' – beam center in the particle reference frame (automatically calculated from 'beam_center_0'
	 *                    defined in InitBeam).
	 * If you need temporary local variables (which are used only in this part of the code), either use 't1'-'t8' or
	 * define your own (with more informative names) in the beginning of this function.
	 */
}
