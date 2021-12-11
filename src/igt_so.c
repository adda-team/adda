/* Functions implementing approximate averaging of the Green's tensor over the voxel volume. 
 * Designed to be conveniently used separately from ADDA as well.
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
#include "igt_so.h" // corresponding headers
// system headers
#include <math.h>

//=====================================================================================================================

static inline double CalcFterm(const double Rvec[static restrict 3],const int mu,const int nu,const double Vr_div_R)
// component (mu,nu) of the F tensor, which is triple antiderivaive of -4pi*k^2*Gst(R)
{
	if (mu==nu) return atan(Vr_div_R/Rvec[mu]/Rvec[mu]);
	else return -atanh(Vr_div_R/Rvec[mu]/Rvec[nu]);
}

//=====================================================================================================================

static inline double CalcBterm(const double Rvec[static restrict 3],const int mu,const int nu,const double R,
	const double Vr,const double Vr_div_R)
// component (mu,nu) of the B tensor, which is triple antiderivaive of 8pi*G1(R)
{
	// TODO: replace log by atanh (as in paper), reuse values of artan and atanh from F (where possible)
	double B_term=0;
	if (mu==nu) {
		double R_tau2;
		for (int tau=0;tau<3;tau++) {
			if (tau==mu) B_term+=2*Vr/Rvec[tau]*log(Rvec[tau] + R);
			else {
				R_tau2=Rvec[tau]*Rvec[tau];
				B_term+=Vr/Rvec[tau]*log(Rvec[tau] + R) - R_tau2*atan(Vr_div_R/R_tau2);
			}
		}
	}
	else B_term=-0.5*(Vr*R/Rvec[mu]/Rvec[nu] + (Rvec[mu]*Rvec[mu]+Rvec[nu]*Rvec[nu])*atanh(Vr_div_R/Rvec[mu]/Rvec[nu]));

	return B_term;
}

//=====================================================================================================================

static inline doublecomplex_t CalcAterm(const int mu,const int nu,doublecomplex_t G,const double R,const double qmunu,
	const double invr3,const double halfk2)
// component (mu,nu) of tensor A=4pi*k^2*(G-Gst-G1)
{
	double Gst,G1;
	if (mu==nu) {
		Gst=-invr3*(1 - 3*qmunu);
		G1=halfk2/R*(1 + qmunu);
	}
	else {
		Gst=invr3*3*qmunu;
		G1=halfk2/R*qmunu;
	}
	return G - Gst - G1;
}

//=====================================================================================================================

static inline void InitUvars(const double Rvec[static restrict 3],const int mu,const int nu,const double ds_x,
	const double ds_y,const double ds_z,const double D,const double invr2,
	double *d_mu2,double *relativeR2,double *D2,double *relativeRrd)
// initialize some variables, which are further used in calculation of U-terms
{
	double d_mu, d_nu;
	// TODO: the following will not be needed if we use dsVec[3] everywhere instead of three variables
	// default is used below to remove warnings of uninitialized
#define SET_d(name) \
	{ \
	switch (name) {\
		case 0: d_##name=ds_x; break;\
		case 1: d_##name=ds_y; break;\
		default: d_##name=ds_z; break;\
	}\
	}
	SET_d(mu);
	SET_d(nu);

#undef SET_d

	// TODO: this seems to be equal to qmunu in other parts of the code => can be optimized
	*relativeR2=Rvec[mu]*Rvec[nu]*invr2;
	*D2=D*D;
	*d_mu2 = d_mu*d_mu;
	double d_nu2=d_nu*d_nu;
	*relativeRrd=invr2;
	if (mu==nu)  *relativeRrd *= 2*(*d_mu2)*Rvec[mu]*Rvec[mu];
	else  *relativeRrd *= Rvec[mu]*Rvec[nu]*((*d_mu2) + d_nu2);
}

//=====================================================================================================================

static inline doublecomplex_t CalcUofA(const double Rvec[static restrict 3],const int mu,const int nu,const double u,
	const double D,const double k2,doublecomplex_t ksi,doublecomplex_t ksi2,doublecomplex_t G,const double _3u_minus_1,
	const double _5u_minus_1,const double _7u_minus_1,const double invr2,const double factor1,const double factor2,
	doublecomplex_t factor3,const double ds_x,const double ds_y,const double ds_z)
/* calculates U(A[mu,nu]), where A=4pi*k^2*(G-Gst-G1), and
 * U(f) = Sum[i=1..3,d_i^2*D[f,r_i,2]] - linear combination of second derivatives
 */
{
	double relativeR2,D2,d_mu2,relativeRrd;
	InitUvars(Rvec,mu,nu,ds_x,ds_y,ds_z,D,invr2,&d_mu2,&relativeR2,&D2,&relativeRrd);

	// first part (G1)
	double first_part=3*relativeR2*_5u_minus_1*D2 - 6*relativeRrd;
	if (mu==nu) first_part+=_3u_minus_1*D2 + 2*d_mu2;
	first_part*=factor1;

	// second part (Gst)
	double second_part=5*relativeR2*_7u_minus_1*D2 - 10*relativeRrd;
	if (mu==nu) second_part+=-_5u_minus_1*D2 + 2*d_mu2;
	second_part*=factor2;

	// third part (G)
	/* TODO: when ksi is small, difference of G-Gst-G1 should be computed term-wise, potentially using Taylor expansion
	 * of exponents
	 */
	doublecomplex_t _3_3ksi_ksi2=3 - 3*ksi + ksi2;
	doublecomplex_t third_part=(relativeR2*_7u_minus_1*D2 - 2*relativeRrd)*(5*_3_3ksi_ksi2 + ksi2*(1 - ksi));
	if (mu==nu) third_part+=(-_5u_minus_1*D2 + 2*d_mu2)*_3_3ksi_ksi2 - _3u_minus_1*D2*ksi2*(1 - ksi);
	third_part*=factor3;
	third_part-=k2*G*u*D2;

	return third_part - second_part - first_part;
}

//=====================================================================================================================

static inline doublecomplex_t CalcUofG(const double Rvec[static restrict 3],const int mu,const int nu,const double u,
	const double D,const double k2,doublecomplex_t ksi,doublecomplex_t ksi2,doublecomplex_t G,const double _3u_minus_1,
	const double _5u_minus_1,const double _7u_minus_1,const double invr2,doublecomplex_t factor3,const double ds_x,
	const double ds_y,const double ds_z)
// calculates U(4pi*k^2*G[mu,nu]), U(f) = Sum[i=1..3,d_i^2*D[f,r_i,2]] - linear combination of second derivatives
{
	double relativeR2,D2,d_mu2,relativeRrd;
	InitUvars(Rvec,mu,nu,ds_x,ds_y,ds_z,D,invr2,&d_mu2,&relativeR2,&D2,&relativeRrd);

	doublecomplex_t _3_3ksi_ksi2=3 - 3*ksi + ksi2;
	doublecomplex_t conv2_G=(relativeR2*_7u_minus_1*D2 - 2*relativeRrd)*(5*_3_3ksi_ksi2 + ksi2*(1 - ksi));
	if (mu==nu) conv2_G+=(-_5u_minus_1*D2 + 2*d_mu2)*_3_3ksi_ksi2 - _3u_minus_1*D2*ksi2*(1 - ksi);

	return factor3*conv2_G - k2*G*u*D2;
}

//=====================================================================================================================

static inline void InitIGTvars(const double rvec[static restrict 3],double qmunu[static restrict 6],
	doublecomplex_t result[static restrict 6],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,double *rr,double *invr3,double *kr,double *kr2,double *D,double *u,double *invvol,double *ds_x2,
	double *ds_y2,double *ds_z2,double *halfk2,double *invr,doublecomplex_t *expval,doublecomplex_t *ksi)
// initialize common variables, inluding point value of G
{
	double qvec[3]; // scaled unit directional vector {rx,ry,rz}/r
	*invvol=1.0/ds_x/ds_y/ds_z;
	*halfk2=0.5*wave_num*wave_num;
	*rr=sqrt(rvec[0]*rvec[0]+rvec[1]*rvec[1]+rvec[2]*rvec[2]);
	*invr=1.0/(*rr);
	*invr3=(*invr)*(*invr)*(*invr);
	*kr=wave_num*(*rr);
	*kr2=(*kr)*(*kr);

	qvec[0]=rvec[0]*(*invr);
	qvec[1]=rvec[1]*(*invr);
	qvec[2]=rvec[2]*(*invr);

	qmunu[0]=qvec[0]*qvec[0];
	qmunu[1]=qvec[0]*qvec[1];
	qmunu[2]=qvec[0]*qvec[2];
	qmunu[3]=qvec[1]*qvec[1];
	qmunu[4]=qvec[1]*qvec[2];
	qmunu[5]=qvec[2]*qvec[2];

	const double t1=(3-(*kr2)), t2=-3*(*kr), t3=((*kr2)-1);

	*expval=(*invr3)*cexp(I*(*kr));

// TODO: make the following style consistent with the rest of the code. Such optimization seems to be too much.
#define INTERACT_DIAG(ind) { result[ind]=((t1*qmunu[ind]+t3) + I*((*kr)+t2*qmunu[ind]))*(*expval); }
#define INTERACT_NONDIAG(ind) { result[ind]=(t1+I*t2)*qmunu[ind]*(*expval); }

	INTERACT_DIAG(0);    // xx
	INTERACT_NONDIAG(1); // xy
	INTERACT_NONDIAG(2); // xz
	INTERACT_DIAG(3);    // yy
	INTERACT_NONDIAG(4); // yz
	INTERACT_DIAG(5);    // zz

#undef INTERACT_DIAG
#undef INTERACT_NONDIAG

	// TODO: make the following three variables internal to this function (not used outside)
	*ds_x2=ds_x*ds_x;
	*ds_y2=ds_y*ds_y;
	*ds_z2=ds_z*ds_z;
	// TODO: check if D is needed anywhere. If only D^2 is used, we can save one sqrt
	*D=sqrt((*ds_x2) + (*ds_y2) + (*ds_z2));
	// TODO: use qvec here instead of rvec ... /rr
	*u=((*ds_x2)*rvec[0]*rvec[0] + (*ds_y2)*rvec[1]*rvec[1]  + (*ds_z2)*rvec[2]*rvec[2])/((*D)*(*D)*(*rr)*(*rr));
	*ksi=I*(*kr);
}

//=====================================================================================================================

static inline void CalcTripleDiffFB(const double rvec[static restrict 3],const double ds_x,const double ds_y,
	const double ds_z,double F_term[static restrict 6],double B_term[static restrict 6])
// computes triple difference of B and F over the corners of the voxel
{
	double signum;
	double Rvec[3];
	int ind0,comp,mu,nu;
	int indX,indY,indZ;
	// TODO: rename Vr to VR and Vr_div_R analogously
	double R,Vr,Vr_div_R;

	for (ind0=0;ind0<6;ind0++) {
		F_term[ind0]=0;
		B_term[ind0]=0;
	}

	// iterate over 8 corners of the dipole (voxel)
	for (indX=-1;indX<2;indX+=2) for (indY=-1;indY<2;indY+=2) for (indZ=-1;indZ<2;indZ+=2) {
		Rvec[0]=rvec[0] + 0.5*indX*ds_x;
		Rvec[1]=rvec[1] + 0.5*indY*ds_y;
		Rvec[2]=rvec[2] + 0.5*indZ*ds_z;
		signum = 1.0*indX*indY*indZ;
		R=sqrt(Rvec[0]*Rvec[0]+Rvec[1]*Rvec[1]+Rvec[2]*Rvec[2]);
		Vr=Rvec[0]*Rvec[1]*Rvec[2];
		// TODO: the following need special care for R->0
		Vr_div_R=Vr/R;
		//iterate over all mu_nu components [3x3] but 6 independent components 2,3 = 3,2; 3,1 = 1,3; 1,2 = 2,1
		for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
			F_term[comp]+=signum*CalcFterm(Rvec,mu,nu,Vr_div_R);
			B_term[comp]+=signum*CalcBterm(Rvec,mu,nu,R,Vr,Vr_div_R);
		}
	}
}

//=====================================================================================================================

void CalcIGTso_near(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6])
// expression for very small distances. The same as medium one, but excluding U(A), which can be too large for small r
{
	double qmunu[6]; // outer-product of qvec (r/|r|) {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr,invr3,kr,kr2; // |r|, |r|^-1, |r|^-3, kr, (kr)^2
	doublecomplex_t expval; // exp(ikr)/|r|^3
	doublecomplex_t ksi; // ikr
	double D; // sqrt(ds_x^2 + ds_y^2 + ds_z^2)
	double u; // sum(ds_x^2*Rx^2 + ds_y^2*Ry^2  + ds_z^2*Rz^2)/(D^2*R^2)
	double invvol; // 1/V
	double ds_x2; // ds_x*ds_x
	double ds_y2; // ds_y*ds_y
	double ds_z2; // ds_z*ds_z
	double halfk2; // k^2/2

	InitIGTvars(rvec,qmunu,result,wave_num,ds_x,ds_y,ds_z,&rr,&invr3,&kr,&kr2,&D,&u,&invvol,&ds_x2,&ds_y2,&ds_z2,
		&halfk2,&invr,&expval,&ksi);

	// local vars
	int comp,mu,nu;
	double F_term[6];
	double B_term[6];
	doublecomplex_t A_term[6];

	CalcTripleDiffFB(rvec,ds_x,ds_y,ds_z,F_term,B_term);

	// iterate over all mu_nu components [3x3] but 6 independent components 2,3 = 3,2; 3,1 = 1,3; 1,2 = 2,1
	for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++)
		A_term[comp]=CalcAterm(mu,nu,result[comp],rr,qmunu[comp],invr3,halfk2);

	for (comp=0;comp<6;comp++) result[comp]=invvol*(halfk2*B_term[comp] - F_term[comp]) + A_term[comp];
}

//=====================================================================================================================

void CalcIGTso_far(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6])
{
	double qmunu[6]; // outer-product of qvec (r/|r|) {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr,invr3,kr,kr2; // |r|, |r|^-1, |r|^-3, kr, (kr)^2
	doublecomplex_t expval; // exp(ikr)/|r|^3
	doublecomplex_t ksi; // ikr
	double D; // sqrt(ds_x^2 + ds_y^2 + ds_z^2)
	double u; // sum(ds_x^2*Rx^2 + ds_y^2*Ry^2  + ds_z^2*Rz^2)/(D^2*R^2)
	double invvol; // 1/V
	double ds_x2; // ds_x*ds_x
	double ds_y2; // ds_y*ds_y
	double ds_z2; // ds_z*ds_z
	double halfk2; // k^2/2

	InitIGTvars(rvec,qmunu,result,wave_num,ds_x,ds_y,ds_z,&rr,&invr3,&kr,&kr2,&D,&u,&invvol,&ds_x2,&ds_y2,&ds_z2,
		&halfk2,&invr,&expval,&ksi);

	// local vars
	int comp,mu,nu;
	doublecomplex_t UofG[6];

	// iterate over all mu_nu components [3x3] but 6 independent components 2,3 = 3,2; 3,1 = 1,3; 1,2 = 2,1
	double _3u_minus_1=3*u - 1;
	double _5u_minus_1=5*u - 1;
	double _7u_minus_1=7*u - 1;
	double R_2=rr*rr;
	double invr2=1/R_2;
	double k2=wave_num*wave_num;
	doublecomplex_t factor3=expval*invr2;
	doublecomplex_t ksi2=-k2*R_2;

	for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++)
		UofG[comp] = CalcUofG(rvec,mu,nu,u,D,k2,ksi,ksi2,result[comp],_3u_minus_1,_5u_minus_1,_7u_minus_1,invr2,
			                  factor3,ds_x,ds_y,ds_z);

	for (comp=0;comp<6;comp++) result[comp]=result[comp] + UofG[comp]/24;
}

//=====================================================================================================================

void CalcIGTso_medium(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6])
// the most extensive formulae to be used for intermediate (medium) range
{
	double qmunu[6]; // outer-product of qvec (r/|r|) {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr,invr3,kr,kr2; // |r|, |r|^-1, |r|^-3, kr, (kr)^2
	doublecomplex_t expval; // exp(ikr)/|r|^3
	doublecomplex_t ksi; // ikr
	double D; // sqrt(ds_x^2 + ds_y^2 + ds_z^2)
	double u; // sum(ds_x^2*Rx^2 + ds_y^2*Ry^2  + ds_z^2*Rz^2)/(D^2*R^2)
	double invvol; // 1/V
	double ds_x2; // ds_x*ds_x
	double ds_y2; // ds_y*ds_y
	double ds_z2; // ds_z*ds_z
	double halfk2; // k^2/2

	InitIGTvars(rvec,qmunu,result,wave_num,ds_x,ds_y,ds_z,&rr,&invr3,&kr,&kr2,&D,&u,&invvol,&ds_x2,&ds_y2,&ds_z2,
		&halfk2,&invr,&expval,&ksi);

	//local vars
	int comp,mu,nu;
	double F_term[6];
	double B_term[6];
	doublecomplex_t A_term[6];
	doublecomplex_t UofA[6];

	CalcTripleDiffFB(rvec,ds_x,ds_y,ds_z,F_term,B_term);

	//iterate over all mu_nu components [3x3] but 6 independent components 2,3 = 3,2; 3,1 = 1,3; 1,2 = 2,1
	double _3u_minus_1=3*u - 1;
	double _5u_minus_1=5*u - 1;
	double _7u_minus_1=7*u - 1;
	double R_2=rr*rr;
	double invr2=1/R_2;
	double k2=wave_num*wave_num;
	double factor1=halfk2*invr3;
	double factor2=3*invr3*invr2;
	doublecomplex_t factor3=expval*invr2;
	doublecomplex_t ksi2=-k2*R_2;

	for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
		A_term[comp]=CalcAterm(mu,nu,result[comp],rr,qmunu[comp],invr3,halfk2);
		UofA[comp]=CalcUofA(rvec,mu,nu,u,D,k2,\
			ksi,ksi2,result[comp],\
			_3u_minus_1,_5u_minus_1,_7u_minus_1,\
			invr2,factor1,factor2,factor3,ds_x,ds_y,ds_z);
	}

	for (comp=0;comp<6;comp++) result[comp]=invvol*(halfk2*B_term[comp] - F_term[comp]) + A_term[comp] + UofA[comp]/24;
}

//=====================================================================================================================

void CalcIGTso(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6])
// main entrance function, see detailed description in the header
{
#define MAX(A,B) (((A) < (B)) ? (B) : (A))
	double d=MAX(ds_x,ds_y); d=MAX(d,ds_z);
#undef MAX
	double r2=rvec[0]*rvec[0] + rvec[1]*rvec[1] + rvec[2]*rvec[2];

	if (wave_num*wave_num*r2 > 25) CalcIGTso_far(rvec,wave_num,ds_x,ds_y,ds_z,result);
	else if (25*r2 > d*d) CalcIGTso_medium(rvec,wave_num,ds_x,ds_y,ds_z,result);
	else CalcIGTso_near(rvec,wave_num,ds_x,ds_y,ds_z,result);
}
