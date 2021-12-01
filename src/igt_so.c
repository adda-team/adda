// system headers
#include <float.h> // for DBL_EPSILON
#include <stdlib.h>
#include <stdbool.h> // for bool
#include <stddef.h> // for size_t
#include <math.h>
#include "igt_so.h"


double calculate_F_term(double Rvec[static 3],\
                        const int mu,\
                        const int nu,\
                        const double Vr_div_R)
{
    if (mu == nu) {
        return atan(Vr_div_R/Rvec[mu]/Rvec[mu]);
    } else {
        return -atanh(Vr_div_R/Rvec[mu]/Rvec[nu]);
    }
}

double calculate_B_term(double Rvec[static 3],\
                        const int mu,\
                        const int nu,\
                        const double R,\
                        const double Vr,\
                        const double Vr_div_R)
{
    double B_term = 0;
    if (mu == nu) {
        double R_tau2;
        for (int tau = 0; tau < 3; tau++) {
            if (tau == mu) {
                B_term += 2*Vr/Rvec[tau]*log(Rvec[tau] + R);
            } else {
                R_tau2 = Rvec[tau]*Rvec[tau];
                B_term += Vr/Rvec[tau]*log(Rvec[tau] + R) - R_tau2*atan(Vr_div_R/R_tau2);
            }
        }

    } else {
        B_term = -0.5*(Vr*R/Rvec[mu]/Rvec[nu] + (Rvec[mu]*Rvec[mu] + Rvec[nu]*Rvec[nu])*atanh(Vr_div_R/Rvec[mu]/Rvec[nu]));
    }
    return B_term;
}

doublecomplex_t calculate_A_term(const int mu,\
                                 const int nu,\
                                 doublecomplex_t G,\
                                 const double R,\
                                 const double qmunu,\
                                 const double invr3,\
                                 const double halfk2)
{
    double Gst, Gext;
    if (mu == nu) {
        Gst = -invr3*(1 - 3*qmunu);
        Gext = halfk2/R*(1 + qmunu);
    } else {
        Gst = invr3*3*qmunu;
        Gext =  halfk2/R*qmunu;
    }
    return G - Gst - Gext;
}

doublecomplex_t calculate_conv2_A_term(double Rvec[static 3],\
                                       const int mu,\
                                       const int nu,\
                                       const double u,\
                                       const double D,\
                                       const double k2,\
                                       doublecomplex_t ksi,\
                                       doublecomplex_t ksi2,\
                                       doublecomplex_t G,\
                                       const double _3u_minus_1,\
                                       const double _5u_minus_1,\
                                       const double _7u_minus_1,\
                                       const double invr2,\
                                       const double factor1,\
                                       const double factor2,\
                                       doublecomplex_t factor3,\
                                       const double ds_x,\
                                       const double ds_y,\
                                       const double ds_z)
{
    double d_mu, d_nu;

#define SET_d(name) \
 { \
    switch (name) {\
    case 0:d_##name=ds_x; break;\
    case 1:d_##name=ds_y; break;\
    case 2:d_##name=ds_z; break;\
    }\
}
    SET_d(mu);
    SET_d(nu);

#undef SET_d

    double relativeR2 = Rvec[mu]*Rvec[nu]*invr2;
    double D2 = D*D;
    double d_mu2 = d_mu*d_mu;
    double d_nu2 = d_nu*d_nu;
    double relativeRrd = invr2;
    if (mu == nu) {
        relativeRrd *= 2*d_mu2*Rvec[mu]*Rvec[mu];
    } else {
        relativeRrd *= Rvec[mu]*Rvec[nu]*(d_mu2 + d_nu2);
    }

    //first part
    //double factor1 = halfk2/R/R_2;
    double first_part = 3*relativeR2*_5u_minus_1*D2 - 6*relativeRrd;
    if (mu == nu) {
        first_part += _3u_minus_1*D2 + 2*d_mu2;
    }
    first_part *= factor1;


    //second part
    //double factor2 = 3/R_2/R_2/R;
    double second_part = 5*relativeR2*_7u_minus_1*D2 - 10*relativeRrd;
    if (mu == nu) {
        second_part += -_5u_minus_1*D2 + 2*d_mu2;
    }
    second_part *= factor2;

    doublecomplex_t _3_3ksi_ksi2 = 3 - 3*ksi + ksi2;
    doublecomplex_t third_part = (relativeR2*_7u_minus_1*D2 - 2*relativeRrd)*(5*_3_3ksi_ksi2 + ksi2*(1 - ksi));
    if (mu == nu) {
        third_part += (-_5u_minus_1*D2 + 2*d_mu2)*_3_3ksi_ksi2 - _3u_minus_1*D2*ksi2*(1 - ksi);
    }
    third_part *= factor3;
    third_part -= k2*G*u*D2;


    return third_part - second_part - first_part;
}

void calc_igt_so(double qvec[static 3],\
                 const double wave_num,\
                 const double ds_x,\
                 const double ds_y,\
                 const double ds_z,\
                 doublecomplex_t result[static restrict 6])
{
    double qvec1[3], qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
    double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
    doublecomplex_t expval; // exp(ikR)/|R|^3
    doublecomplex_t ksi; // ikR
    int ind0,comp,mu,nu;//,indmunu,mu1,nu1;
    int indX,indY,indZ;
    double D;//sqrt(ds_x^2 + ds_y^2 + ds_z^2)
    double u;//sum(ds_x^2*Rx^2 + ds_y^2*Ry^2  + ds_z^2*Rz^2)/(D^2*R^2)
    double Rvec[3];
    double signum;
    double invvol;
    double ds_x2;
    double ds_y2;
    double ds_z2;
    double halfk2;
    double F_term[6];
    double B_term[6];
    doublecomplex_t A_term[6];
    doublecomplex_t conv2_A_term[6];
    double invr;

    halfk2 = 0.5*wave_num*wave_num;
    rr=sqrt(qvec[0]*qvec[0]+qvec[1]*qvec[1]+qvec[2]*qvec[2]);
    invr = 1.0/rr;
    invr3=invr*invr*invr;
    kr=wave_num*rr;
    kr2=kr*kr;

    qvec1[0] = qvec[0]*invr;
    qvec1[1] = qvec[1]*invr;
    qvec1[2] = qvec[2]*invr;

    qmunu[0] = qvec1[0]*qvec1[0];
    qmunu[1] = qvec1[0]*qvec1[1];
    qmunu[2] = qvec1[0]*qvec1[2];
    qmunu[3] = qvec1[1]*qvec1[1];
    qmunu[4] = qvec1[1]*qvec1[2];
    qmunu[5] = qvec1[2]*qvec1[2];

    //InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

    const double t1=(3-kr2), t2=-3*kr, t3=(kr2-1);

    expval=invr3*cexp(I*kr);

#define INTERACT_DIAG(ind) { result[ind] = ((t1*qmunu[ind]+t3) + I*(kr+t2*qmunu[ind]))*expval; }
#define INTERACT_NONDIAG(ind) { result[ind] = (t1+I*t2)*qmunu[ind]*expval; }

    INTERACT_DIAG(0);    // xx
    INTERACT_NONDIAG(1); // xy
    INTERACT_NONDIAG(2); // xz
    INTERACT_DIAG(3);    // yy
    INTERACT_NONDIAG(4); // yz
    INTERACT_DIAG(5);    // zz

#undef INTERACT_DIAG
#undef INTERACT_NONDIAG

    ds_x2 = ds_x*ds_x;
    ds_y2 = ds_y*ds_y;
    ds_z2 = ds_z*ds_z;
    D = sqrt(ds_x2 + ds_y2 + ds_z2);
    u=(ds_x2*qvec[0]*qvec[0] + ds_y2*qvec[1]*qvec[1]  + ds_z2*qvec[2]*qvec[2])/(D*D*rr*rr);
    ksi=I*kr;

    ds_x2 = ds_x*ds_x;
    ds_y2 = ds_y*ds_y;
    ds_z2 = ds_z*ds_z;
    D = sqrt(ds_x2 + ds_y2 + ds_z2);
    u=(ds_x2*qvec[0]*qvec[0] + ds_y2*qvec[1]*qvec[1]  + ds_z2*qvec[2]*qvec[2])/(D*D*rr*rr);
    ksi=I*kr;

    for (ind0=0;ind0<6;ind0++) {
        F_term[ind0] = 0;
        B_term[ind0] = 0;
    }

    //iterate over corners of the ij dipole, 8 corners
    for (indX=-1;indX<2;indX+=2) for (indY=-1;indY<2;indY+=2) for (indZ=-1;indZ<2;indZ+=2) {
        Rvec[0] = qvec[0] + 0.5*indX*ds_x;
        Rvec[1] = qvec[1] + 0.5*indY*ds_y;
        Rvec[2] = qvec[2] + 0.5*indZ*ds_z;
        signum = 1.0*indX*indY*indZ;
        double R = sqrt(Rvec[0]*Rvec[0]+Rvec[1]*Rvec[1]+Rvec[2]*Rvec[2]);
        double Vr = Rvec[0]*Rvec[1]*Rvec[2];
        double Vr_div_R = Vr/R;
        //iterate over all mu_nu components [3x3] but 6 independent components 2,3 = 3,2; 3,1 = 1,3; 1,2 = 2,1
        for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
            F_term[comp] += signum*calculate_F_term(Rvec, mu, nu, Vr_div_R);
            B_term[comp] += signum*calculate_B_term(Rvec, mu, nu, R, Vr, Vr_div_R);
        }
    }

    //iterate over all mu_nu components [3x3] but 6 independent components 2,3 = 3,2; 3,1 = 1,3; 1,2 = 2,1
    double _3u_minus_1 = 3*u - 1;
    double _5u_minus_1 = 5*u - 1;
    double _7u_minus_1 = 7*u - 1;
    double R_2 = rr*rr;
    double invr2 = 1/R_2;
    double k2 = wave_num*wave_num;
    double factor1 = halfk2*invr3;
    double factor2 = 3*invr3*invr2;
    doublecomplex_t factor3 = expval*invr2;
    doublecomplex_t ksi2 = -k2*R_2;

    for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
        A_term[comp] = calculate_A_term(mu, nu, result[comp], rr, qmunu[comp], invr3, halfk2);
        conv2_A_term[comp] = calculate_conv2_A_term(qvec, mu, nu, u, D, k2,\
                                                            ksi, ksi2, result[comp],\
                                                            _3u_minus_1, _5u_minus_1, _7u_minus_1,\
                                                            invr2, factor1, factor2, factor3, ds_x, ds_y, ds_z);
    }

    invvol = 1.0/ds_x/ds_y/ds_z;

    //0.0416666666666667 is 1/24;
    for (comp=0;comp<6;comp++) {
        result[comp] = invvol*(halfk2*B_term[comp] - F_term[comp]) + A_term[comp] + 0.0416666666666667*conv2_A_term[comp];
    }
}