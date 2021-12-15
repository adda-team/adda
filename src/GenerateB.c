/* Generates the incident beam
 *
 * Lminus beam is based on: G. Gouesbet, B. Maheu, G. Grehan, "Light scattering from a sphere arbitrary located
 * in a Gaussian beam, using a Bromwhich formulation", J.Opt.Soc.Am.A 5,1427-1443 (1988).
 * Eq.(22) - complex conjugate.
 *
 * Davis beam is based on: L. W. Davis, "Theory of electromagnetic beams," Phys.Rev.A 19, 1177-1179 (1979).
 * Eqs.(15a),(15b) - complex conjugate; in (15a) "Q" changed to "Q^2" (typo).
 *
 * Barton beam is based on: J. P. Barton and D. R. Alexander, "Fifth-order corrected electromagnetic-field components
 * for a fundamental Gaussian-beam," J.Appl.Phys. 66,2800-2802 (1989).
 * Eqs.(25)-(28) - complex conjugate.
 *
 * Copyright (C) ADDA contributors
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
#include "interaction.h"
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
extern const opt_index opt_beam;

// used in CalculateE.c
double C0dipole,C0dipole_refl; // inherent cross sections of exciting dipole (in free space and addition due to surface)
// used in param.c
const char *beam_descr; // string for log file with beam parameters

// LOCAL VARIABLES
static double s,s2;            // beam confinement factor and its square
static double scale_x,scale_z; // multipliers for scaling coordinates
static doublecomplex ki,kt;    // abs of normal components of k_inc/k0, and ktran/k0
static doublecomplex ktVec[3]; // k_tran/k0
static double p0;              // amplitude of the incident dipole moment
#ifndef NO_FORTRAN
static int n0;                 // Bessel beam order
int vorticity;				   // Vorticity of vortex beams (n0 for Bessel beams)
static double alpha0;          // half-cone angle
static doublecomplex K,Kt,Kz;  // wave-vector components
static doublecomplex M[4];     // Components of matrix M
#endif
/* TO ADD NEW BEAM
 * Add here all internal variables (beam parameters), which you initialize in InitBeam() and use in GenerateB()
 * afterwards. If you need local, intermediate variables, put them into the beginning of the corresponding function.
 * Add descriptive comments, use 'static'.
 */
 
// EXTERNAL FUNCTIONS

#ifndef NO_FORTRAN 
void bjndd_(const int *n,const double *x,double *bj,double *dj,double *fj);
#endif

//======================================================================================================================

void InitBeam(void)
// initialize beam; produce description string
{
	double w0; // beam width
	const char *tmp_str; // temporary string
	/* TO ADD NEW BEAM
	 * Add here all intermediate variables, which are used only inside this function.
	 */
	
	// initialization of global option index for error messages
	opt=opt_beam;
	// beam initialization
	beam_asym=(beam_center_0[0]!=0 || beam_center_0[1]!=0 || beam_center_0[2]!=0);
	if (!beam_asym) vInit(beam_center); // not to calculate it for each orientation
	vorticity = 0;
	switch (beamtype) {
		case B_PLANE:
			if (IFROOT) beam_descr="plane wave";
			if (surface) {
				/* here we assume that prop_0 will not change further (e.g., by rotation of particle), 
				 * i.e. prop=prop_0 in GenerateBeam() below
				 */				
				if (prop_0[2]==0) PrintError("Ambiguous setting of beam propagating along the surface. Please specify "
					"the incident direction to have (arbitrary) small positive or negative z-component");
				if (msubInf && prop_0[2]>0) PrintError("Perfectly reflecting surface ('-surf ... inf') is incompatible "
					"with incident direction from below (including the default one)");
				// Here we set ki,kt,ktVec and propagation directions prIncRefl,prIncTran
				if (prop_0[2]>0) { // beam comes from the substrate (below)
					// here msub should always be defined
					inc_scale=1/creal(msub);
					ki=msub*prop_0[2];
					/* Special case for msub near 1 to remove discontinuities for near-grazing incidence. The details
					 * are discussed in CalcFieldSurf() in crosssec.c.
					 */
					if (cabs(msub-1)<ROUND_ERR && cabs(ki)<SQRT_RND_ERR) kt=ki;
					else kt=cSqrtCut(1 - msub*msub*(prop_0[0]*prop_0[0]+prop_0[1]*prop_0[1]));
					// determine propagation direction and full wavevector of wave transmitted into substrate
					ktVec[0]=msub*prop_0[0];
					ktVec[1]=msub*prop_0[1];
					ktVec[2]=kt;
				}
				else if (prop_0[2]<0) { // beam comes from above the substrate
					inc_scale=1;
					vRefl(prop_0,prIncRefl);
					ki=-prop_0[2]; // always real
					if (!msubInf) {
						// same special case as above
						if (cabs(msub-1)<ROUND_ERR && cabs(ki)<SQRT_RND_ERR) kt=ki;
						else kt=cSqrtCut(msub*msub - (prop_0[0]*prop_0[0]+prop_0[1]*prop_0[1]));
						// determine propagation direction of wave transmitted into substrate
						ktVec[0]=prop_0[0];
						ktVec[1]=prop_0[1];
						ktVec[2]=-kt;
					}
				}
				else LogError(ONE_POS,"Ambiguous setting of beam propagating along the surface. Please specify the"
					"incident direction to have (arbitrary) small positive or negative z-component");
				vRefl(prop_0,prIncRefl);
				if (!msubInf) {
					vReal(ktVec,prIncTran);
					vNormalize(prIncTran);
				}
			}
			return;
		case B_DIPOLE:
			if (surface) {
				if (beam_center_0[2]<=-hsub)
					PrintErrorHelp("External dipole should be placed strictly above the surface");
				inc_scale=1; // but scaling of Mueller matrix is weird anyway
			}
			/* definition of p0 is important for scaling of many scattering quantities (that are normalized to incident
			 * irradiance). Alternative definition is p0=1, but then the results will scale with unit of length
			 * (breaking scale invariance)
			 */
			p0=1/(WaveNum*WaveNum*WaveNum);
			if (IFROOT) beam_descr="point dipole";
			return;
		case B_LMINUS:
		case B_DAVIS3:
		case B_BARTON5:
			if (surface) PrintError("Currently, Gaussian incident beam is not supported for '-surf'");
			// initialize parameters
			w0=beam_pars[0];
			TestPositive(w0,"beam width");
			s=1/(WaveNum*w0);
			s2=s*s;
			scale_x=1/w0;
			scale_z=s*scale_x; // 1/(k*w0^2)
			// beam info
			if (IFROOT) {
				switch (beamtype) {
					case B_LMINUS:
						tmp_str="L- approximation";
						break;
					case B_DAVIS3:
						tmp_str="3rd order approximation, by Davis";
						break;
					case B_BARTON5:
						tmp_str="5th order approximation, by Barton";
						break;
					default: LogError(ONE_POS,"Incompatibility error in GenerateB");
				}
				beam_descr=dyn_sprintf("Gaussian beam (%s)\n"
				                       "\tWidth="GFORMDEF" (confinement factor s="GFORMDEF")\n",tmp_str,w0,s);
			}
			return;
#ifndef NO_FORTRAN
		case B_BES_CS:
		case B_BES_CSp:
		case B_BES_M:
		case B_BES_LE:
		case B_BES_LM:
		case B_BES_TEL:
		case B_BES_TML:
			if (surface) PrintError("Currently, Bessel incident beam is not supported for '-surf'");
			// initialize parameters
			ConvertToInteger(beam_pars[0],"beam order",&n0);
			TestRangeII(abs(n0),"beam order (might cause the incorrect calculation of Bessel function)",0,899);
			vorticity = n0;
			alpha0 = Deg2Rad(beam_pars[1]);
			K =fabs(WaveNum);
			Kt=fabs(WaveNum)*sin(alpha0);
			Kz=fabs(WaveNum)*cos(alpha0);
			switch (beamtype) { // definition of elements of matrix M ((Mex,Mey),(Mmx,Mmy))
				case B_BES_CS:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					M[0]=0.5;	M[1]=0;
					M[2]=0;		M[3]=0.5;
					break;
				case B_BES_CSp:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					M[0]=0.5;	M[1]=0;
					M[2]=0;		M[3]=-0.5;
					break;
				case B_BES_M:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					if (beam_Npars==6) {
						M[0]=beam_pars[2];	M[1]=beam_pars[3];
						M[2]=beam_pars[4];	M[3]=beam_pars[5];
					}
					else {
						M[0]=beam_pars[2]+I*beam_pars[6];
						M[1]=beam_pars[3]+I*beam_pars[7];
						M[2]=beam_pars[4]+I*beam_pars[8];
						M[3]=beam_pars[5]+I*beam_pars[9];
					}
					break;
				case B_BES_LE:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					M[0]=0;		M[1]=0;
					M[2]=0;		M[3]=1;
					break;
				case B_BES_LM:
					TestRangeII(beam_pars[1],"half-cone angle",0,90);
					M[0]=0;		M[1]=1;
					M[2]=0;		M[3]=0;
					break;
				case B_BES_TEL:
					TestRangeNN(beam_pars[1],"half-cone angle for TEL type",0,90);
					M[0]=-K/Kt;		M[1]=0;
					M[2]=0;			M[3]=Kz/Kt;
					break;
				case B_BES_TML:
					TestRangeNN(beam_pars[1],"half-cone angle for TML type",0,90);
					M[0]=0; 		M[1]=Kz/Kt;
					M[2]=K/Kt;		M[2]=0;
					break;
				default: LogError(ONE_POS,"Incompatibility error in GenerateB");
			}
			symR=symX=symY=symZ=false;
			// beam info
			if (IFROOT) {
				switch (beamtype) {
					case B_BES_CS:
						tmp_str="circularly symmetric energy density";
						break;
					case B_BES_CSp:
						tmp_str="circularly symmetric energy density, alternative type";
						break;
					case B_BES_M:
						tmp_str="generalized";
						break;
					case B_BES_LE:
						tmp_str="linear electric field";
						break;
					case B_BES_LM:
						tmp_str="linear magnetic field";
						break;
					case B_BES_TEL:
						tmp_str="linear component of the TE";
						break;
					case B_BES_TML:
						tmp_str="linear component of the TM";
						break;
					default: LogError(ONE_POS,"Incompatibility error in GenerateB");
				}
				beam_descr=dyn_sprintf("Bessel beam (%s)\n"
									   "\tOrder: ""%d""\n\t""Half-cone angle: "GFORMDEF"\n",tmp_str,n0,alpha0);
			}
			return;
#endif // !NO_FORTRAN
		case B_READ:
			// the safest is to assume cancellation of all symmetries
			symX=symY=symZ=symR=false;
			if (surface) inc_scale=1; // since we can't know it, we assume the default case
			if (IFROOT) {
				if (beam_Npars==1) beam_descr=dyn_sprintf("specified by file '%s'",beam_fnameY);
				else beam_descr=dyn_sprintf("specified by files '%s' and '%s'",beam_fnameY,beam_fnameX);
			}
			// this type is unaffected by beam_center
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. The case should
	 * 1) save all the input parameters from array 'beam_pars' to local variables (defined in the beginning of this
	 *    source files)
	 * 2) test all input parameters (for that you're encouraged to use functions from param.h since they would
	 *    automatically produce informative output in case of error).
	 * 3) the symmetry breaking due to prop or beam_center is taken care of in VariablesInterconnect() in param.c.
	 *    But if there are other reasons why beam would break any symmetry, corresponding variable should be set to
	 *    false here. Do not set any of them to true, as they can be set to false by other factors.
	 *    symX, symY, symZ - symmetries of reflection over planes YZ, XZ, XY respectively.
	 *    symR - symmetry of rotation for 90 degrees over the Z axis
	 * 4) initialize beam_descr - descriptive string, which will appear in log file.
	 * 5) Consider the case of surface (substrate near the particle). If the new beam type is incompatible with it, add
	 *    an explicit exception, like "if (surface) PrintErrorHelp(...);". Otherwise, you also need to define inc_scale.
	 * All other auxiliary variables, which are used in beam generation (GenerateB(), see below), should be defined in
	 * the beginning of this file. If you need temporary local variables (which are used only in this part of the code),
	 * define them in the beginning of this function.
	 */
}

//======================================================================================================================

void GenerateB (const enum incpol which,   // x - or y polarized incident light
                doublecomplex *restrict b) // the b vector for the incident field
// generates incident beam at every dipole
{
	size_t i,j;
	doublecomplex psi0,Q,Q2;
	doublecomplex v1[3],v2[3],v3[3],gt[6];
	double ro2,ro4;
	double x,y,z,x2_s,xy_s;
	doublecomplex t1,t2,t3,t4,t5,t6,t7,t8,ctemp;
	const double *ex; // coordinate axis of the beam reference frame
	double ey[3];
	double r1[3];
	/* complex wave amplitudes of secondary waves (with phase relative to particle center);
	 * The transmitted wave can be inhomogeneous wave (when msub is complex), then eIncTran (e) is normalized
	 * counter-intuitively. Before multiplying by tc, it satisfies (e,e)=1!=||e||^2. This normalization is consistent
	 * with used formulae for transmission coefficients. So this transmission coefficient is not (generally) equal to
	 * the ratio of amplitudes of the electric fields. In particular, when E=E0*e, ||E||!=|E0|*||e||, where
	 * ||e||^2=(e,e*)=|e_x|^2+|e_y|^2+|e_z|^2=1
	 */
	doublecomplex eIncRefl[3],eIncTran[3];
#ifndef NO_FORTRAN
	double ro;
	int n1,q; // for Bessel beams
	doublecomplex vort;  // vortex phase of Bessel beam rotated on 90 deg
	doublecomplex fn[5]; // for general functions f(n,ro,phi) of Bessel beams (fn-2, fn-1, fn, fn+1, fn+2 respectively)
	double phi,arg,td1[abs(n0)+3],td2[abs(n0)+3],jn1[abs(n0)+3]; // for Bessel beams
#endif
	const char *fname;
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
		vCopy(incPolY,ey);
	}

	switch (beamtype) {
		case B_PLANE: // plane is separate to be fast (for non-surface)
			if (surface) {
				/* With respect to normalization we use here the same assumption as in the free space - the origin is in
				 * the particle center, and amplitude of incoming plane wave is equal to 1. Then irradiance of the beam
				 * coming from below is c*Re(msub)/(8pi), different from that coming from above.
				 * Original incident (incoming) beam propagating from the vacuum (above) is Exp(i*k*r.a), while - from
				 * the substrate (below) is Exp(i*k*msub*r.a). We assume that the incoming beam is homogeneous in its
				 * original medium.
				 */
				doublecomplex rc,tc; // reflection and transmission coefficients
				double hbeam=hsub+beam_center[2]; // height of beam center above the surface 
				if (prop[2]>0) { // beam comes from the substrate (below)
					//  determine amplitude of the reflected and transmitted waves; here msub is always defined
					if (which==INCPOL_Y) { // s-polarized
						cvBuildRe(ex,eIncTran);
						tc=FresnelTS(ki,kt);
					}
					else { // p-polarized
						crCrossProd(ey,ktVec,eIncTran);
						tc=FresnelTP(ki,kt,1/msub);
					}
					// phase shift due to the beam center relative to the origin and surface
					cvMultScal_cmplx(tc*cexp(I*WaveNum*((kt-ki)*hbeam-crDotProd(ktVec,beam_center))),eIncTran,eIncTran);
					// main part
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						// b[i] = eIncTran*exp(ik*kt.r); subtraction of beam_center is avoided by the factor above
						cvMultScal_cmplx(cexp(I*WaveNum*crDotProd(ktVec,DipoleCoord+j)),eIncTran,b+j);
					}
				}
				else if (prop[2]<0) { // beam comes from above the substrate
					// determine amplitude of the reflected and transmitted waves
					if (which==INCPOL_Y) { // s-polarized
						cvBuildRe(ex,eIncRefl);
						if (msubInf) rc=-1;
						else {
							cvBuildRe(ex,eIncTran);
							rc=FresnelRS(ki,kt);
						}
					}
					else { // p-polarized
						vInvRefl_cr(ex,eIncRefl);
						if (msubInf) rc=1;
						else {
							crCrossProd(ey,ktVec,eIncTran);
							cvMultScal_cmplx(1/msub,eIncTran,eIncTran); // normalize eIncTran by ||ktVec||=msub
							rc=FresnelRP(ki,kt,msub);
						}
					}
					// phase shift due to the beam center relative to the surface
					cvMultScal_cmplx(rc*imExp(2*WaveNum*creal(ki)*hbeam),eIncRefl,eIncRefl); // assumes real ki
					// main part
					for (i=0;i<local_nvoid_Ndip;i++) {
						j=3*i;
						// together with above phase shift, beam_center affects only the common phase factor
						LinComb(DipoleCoord+j,beam_center,1,-1,r1); 
						// b[i] = ex*exp(ik*r.a) + eIncRefl*exp(ik*prIncRefl.r)
						cvMultScal_RVec(imExp(WaveNum*DotProd(r1,prop)),ex,b+j);
						cvLinComb1_cmplx(eIncRefl,b+j,imExp(WaveNum*DotProd(r1,prIncRefl)),b+j);
					}
				}
			}
			else for (i=0;i<local_nvoid_Ndip;i++) { // standard (non-surface) plane wave
				j=3*i;
				LinComb(DipoleCoord+j,beam_center,1,-1,r1); // this affects only the common phase factor
				ctemp=imExp(WaveNum*DotProd(r1,prop)); // ctemp=exp(ik*r.a)
				cvMultScal_RVec(ctemp,ex,b+j); // b[i]=ctemp*ex
			}
			return;
		case B_DIPOLE: {
			double dip_p[3]; // dipole moment, = p0*prop
			vMultScal(p0,prop,dip_p);
			for (i=0;i<local_nvoid_Ndip;i++) { // here we explicitly use that dip_p is real
				j=3*i;
				LinComb(DipoleCoord+j,beam_center,1,-1,r1);
				(*InterTerm_real)(r1,gt);
				cSymMatrVecReal(gt,dip_p,b+j);
				if (surface) { // add reflected field
					r1[2]=DipoleCoord[j+2]+beam_center[2]+2*hsub;
					(*ReflTerm_real)(r1,gt);
					cReflMatrVecReal(gt,dip_p,v1);
					cvAdd(v1,b+j,b+j);
				}
			}
			/* calculate dipole inherent cross sections (for decay rate enhancements)
			 * General formula is C0=4pi*k*Im(p0*.G(r0,r0).p0) and it is used for reflection; but for direct
			 * interaction it is expected to be singular, so we use an explicit formula for point dipole:
			 * C0=(8pi/3)*k^4*|p0|^2. Thus we discard the choice of "-int ...", but still the used value should be
			 * correct for most formulations, e.g. poi, fcd, fcd_st, igt_so. Moreover, it is also logical, since the
			 * exciting dipole is really point one, in contrast to the dipoles composing the particle.
			 */
			double temp=p0*WaveNum*WaveNum; // in principle, t1=1/k, but we keep a general formula
			C0dipole=2*FOUR_PI_OVER_THREE*temp*temp;
			if (surface) {
				r1[0]=r1[1]=0;
				r1[2]=2*(beam_center[2]+hsub);
				(*ReflTerm_real)(r1,gt);
				double tmp;
				/* the following expression uses that dip_p is real and a specific (anti-)symmetry of the gt
				 * a general expression is commented out below
				 */
				tmp=dip_p[0]*dip_p[0]*cimag(gt[0])+2*dip_p[0]*dip_p[1]*cimag(gt[1])+dip_p[1]*dip_p[1]*cimag(gt[3])
				   +dip_p[2]*dip_p[2]*cimag(gt[5]);
//				v1[0]=dip_p[0]; v1[1]=dip_p[1]; v1[2]=dip_p[2];
//				cReflMatrVec(gt,v1,v2);
//				tmp=cDotProd_Im(v2,v1);
				C0dipole_refl=FOUR_PI*WaveNum*tmp;
			}
			return;
		}
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
				Q=1/(2*z-I);
				psi0=-I*Q*cexp(I*Q*ro2);
				// ctemp=exp(ik*z0)*psi0, z0 - non-scaled coordinate (z/scale_z)
				ctemp=imExp(WaveNum*z/scale_z)*psi0;
				// the following logic (if-else-if...) is hard to replace by a simple switch
				if (beamtype==B_LMINUS) cvMultScal_RVec(ctemp,ex,b+j); // b[i]=ctemp*ex
				else {
					/* It is possible to rewrite the formulae below to avoid division by ro2, but we prefer
					 * dimensionless variables. The value for ro2=0 doesn't really matter (cancels afterwards).
					 * The current code should work OK even for very small ro2
					 */
					if (ro2==0) x2_s=0;
					else x2_s=x*x/ro2;
					Q2=Q*Q;
					ro4=ro2*ro2;
					// some combinations that are used more than once
					t4=s2*ro2*Q2; // t4=(s*ro*Q)^2
					t5=I*ro2*Q;   // t5=i*Q*ro^2
					t6=ro4*Q2;    // t6=ro^4*Q^2
					t7=x*s*Q;     // t7=x*s*Q
					if (beamtype==B_DAVIS3) {
						// t1=1+s^2(-4Q^2*x^2-iQ^3*ro^4)=1-t4(4x2_s+t5)
						t1 = 1 - t4*(4*x2_s+t5);
						// t2=0
						t2=0;
						// t3=-s(2Qx)+s^3(8Q^3*ro^2*x+2iQ^4*ro^4*x-4iQ^2x)=2t7[-1+iQ*s2*(-4t5+t6-2)]
						t3 = 2*t7*(-1 + I*Q*s2*(-4*t5+t6-2));
					}
					else if (beamtype==B_BARTON5) {
						if (ro2==0) xy_s=0; // see comment for x2_s above
						else xy_s=x*y/ro2;
						t8=8+2*t5; // t8=8+2i*Q*ro^2
						/* t1 = 1 + s^2(-ro^2*Q^2-i*ro^4*Q^3-2Q^2*x^2)
						 *    + s^4[2ro^4*Q^4+3iro^6*Q^5-0.5ro^8*Q^6+x^2(8ro^2*Q^4+2iro^4*Q^5)]
						 *    = 1 + t4*{-1-2x2_s-t5+t4*[2+3t5-0.5t6+x2_s*t8]}
						 */
						t1 = 1 + t4*(-1 - 2*x2_s - t5 + t4*(2+3*t5-0.5*t6+x2_s*t8));
						// t2=s^2(-2Q^2*xy)+s^4[xy(8ro^2*Q^4+2iro^4*Q^5)]=xy_s*t4(-2+t4*t8)
						t2=xy_s*t4*(-2+t4*t8);
						/* t3 = s(-2Qx) + s^3[(6ro^2*Q^3+2iro^4*Q^4)x] + s^5[(-20ro^4*Q^5-10iro^6*Q^6+ro^8*Q^7)x]
						 *    = t7{-2+t4[6+2t5+t4(-20-10t5+t6)]}
						 */
						t3 = t7*(-2 + t4*(6 + 2*t5 + t4*(-20-10*t5+t6)));
					}
					else LogError(ONE_POS,"Inconsistency in beam definition"); // to remove warnings
					// b[i]=ctemp(ex*t1+ey*t2+ez*t3)
					cvMultScal_RVec(t1,ex,v1);
					cvMultScal_RVec(t2,ey,v2);
					cvMultScal_RVec(t3,prop,v3);
					cvAdd2Self(v1,v2,v3);
					cvMultScal_cmplx(ctemp,v1,b+j);
				}
			}
			return;
#ifndef NO_FORTRAN
		case B_BES_CS:
		case B_BES_CSp:
		case B_BES_M:
		case B_BES_LE:
		case B_BES_LM:
		case B_BES_TEL:
		case B_BES_TML:
			if (which==INCPOL_X) vort = cpow(I,n0);
			else vort = 1.;
			for (i=0;i<local_nvoid_Ndip;i++) {
				j=3*i;
				LinComb(DipoleCoord+j,beam_center,1,-1,r1);
				x=DotProd(r1,ex);
				y=DotProd(r1,ey);
				z=DotProd(r1,prop);
				ro2=x*x+y*y;
				ro=sqrt(ro2);	// radial distance in a cylindrical coordinate system
				phi=atan2(y,x);	// angular coordinate in a cylindrical coordinate system
				// common factor
				ctemp=cexp(I*n0*phi)*cexp(I*Kz*z)/K/K*vort;
				arg=Kt*ro;
				if (arg<ROUND_ERR){
					if (n0 == 0) fn[2]=1.;
					else fn[2]=0;
					fn[0]=0; fn[1]=0; fn[3]=0; fn[4]=0;
				}
				else {
					n1=abs(n0)+2;
					if (n0 <= -3) {
						bjndd_(&n1, &arg, jn1, td1, td2);
						if (n1 % 2 == 0) q=1;	else q=-1;
						fn[0] =  q*jn1[-n0+2]*cexp(-2*I*phi);
						fn[1] = -q*jn1[-n0+1]*cexp(-I*phi);
						fn[2] =  q*jn1[-n0];
						fn[3] = -q*jn1[-n0-1]*cexp(I*phi);
						fn[4] =  q*jn1[-n0-2]*cexp(2*I*phi);
					}
					if (n0 >= 2) {
						bjndd_(&n1, &arg, jn1, td1, td2);
						fn[0] = jn1[n0-2]*cexp(-2*I*phi);
						fn[1] = jn1[n0-1]*cexp(-I*phi);
						fn[2] = jn1[n0];
						fn[3] = jn1[n0+1]*cexp(I*phi);
						fn[4] = jn1[n0+2]*cexp(2*I*phi);
					}
					if (n0 == -2) {
						bjndd_(&n1, &arg, jn1, td1, td2);
						fn[0] =  jn1[4]*cexp(-2*I*phi);
						fn[1] = -jn1[3]*cexp(-I*phi);
						fn[2] =  jn1[2];
						fn[3] = -jn1[1]*cexp(I*phi);
						fn[4] =  jn1[0]*cexp(2*I*phi);
					}
					if (n0 == -1) {
						bjndd_(&n1, &arg, jn1, td1, td2);
						fn[0] = -jn1[3]*cexp(-2*I*phi);
						fn[1] =  jn1[2]*cexp(-I*phi);
						fn[2] = -jn1[1];
						fn[3] =  jn1[0]*cexp(I*phi);
						fn[4] =  jn1[1]*cexp(2*I*phi);
					}
					if (n0 == 0) {
						bjndd_(&n1, &arg, jn1, td1, td2);
						fn[0] =  jn1[2]*cexp(-2*I*phi);
						fn[1] = -jn1[1]*cexp(-I*phi);
						fn[2] =  jn1[0];
						fn[3] =  jn1[1]*cexp(I*phi);
						fn[4] =  jn1[2]*cexp(2*I*phi);
					}
					if (n0 == 1) {
						bjndd_(&n1, &arg, jn1, td1, td2);
						fn[0] = -jn1[1]*cexp(-2*I*phi);
						fn[1] =  jn1[0]*cexp(-I*phi);
						fn[2] =  jn1[1];
						fn[3] =  jn1[2]*cexp(I*phi);
						fn[4] =  jn1[3]*cexp(2*I*phi);
					}
				}

				t1 = (((K*K+Kz*Kz)/2.*M[0] +  K*Kz*M[3])*fn[2] +
						  Kt*Kt/4.*((M[0]+I*M[1])*fn[0] + (M[0]-I*M[1])*fn[4]));	// Ex
				t2 = (((K*K+Kz*Kz)/2.*M[1] -  K*Kz*M[2])*fn[2] +
						I*Kt*Kt/4.*((M[0]+I*M[1])*fn[0] - (M[0]-I*M[1])*fn[4]));	// Ey
				t3 = ((I*Kz*(M[0]+I*M[1]) + K*(M[2]+I*M[3]))*fn[1] -
					  (I*Kz*(M[0]-I*M[1]) - K*(M[2]-I*M[3]))*fn[3])*Kt/2.;			// Ez

				cvMultScal_RVec(t1,ex,v1);
				cvMultScal_RVec(t2,ey,v2);
				cvMultScal_RVec(t3,prop,v3);
				cvAdd2Self(v1,v2,v3);
				cvMultScal_cmplx(ctemp,v1,b+j);
			}
			return;
#endif // !NO_FORTRAN		
		case B_READ:
			if (which==INCPOL_Y) fname=beam_fnameY;
			else fname=beam_fnameX; // which==INCPOL_X
			ReadField(fname,b);
			return;
	}
	LogError(ONE_POS,"Unknown type of incident beam (%d)",(int)beamtype);
	/* TO ADD NEW BEAM
	 * add a case above. Identifier ('B_...') should be defined inside 'enum beam' in const.h. This case should set
	 * complex vector 'b', describing the incident field in the particle reference frame. It is set inside the cycle for
	 * each dipole of the particle and is calculated using
	 * 1) 'DipoleCoord' - array of dipole coordinates;
	 * 2) 'prop' - propagation direction of the incident field;
	 * 3) 'ex' - direction of incident polarization;
	 * 4) 'ey' - complementary unity vector of polarization (orthogonal to both 'prop' and 'ex');
	 * 5) 'beam_center' - beam center in the particle reference frame (automatically calculated from 'beam_center_0'
	 *                    defined by '-beam_center' command line option).
	 * If the new beam type is compatible with '-surf', include here the corresponding code. For that you will need
	 * the variables, related to surface - see vars.c after "// related to a nearby surface".
	 * If you need temporary local variables (which are used only in this part of the code), either use 't1'-'t8' or
	 * define your own (with more informative names) in the beginning of this function.
	 */
}
