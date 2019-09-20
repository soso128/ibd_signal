#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <string>
#include <map>

using namespace std;

/****************************************************************************/
/*          This code contains useful functions to generate the             */
/*    DSNB spectrum as well as the IBD spectrum in SK (22.5kton.year)       */
/*    It is used to generate a library than can then be imported as a       */
/*    Python module but it can also be used as a standalone if a main       */
/*    method is added.                                                      */
/****************************************************************************/

/**Prototypes ***************************************************************/
double dsigma(double e_nu,double e_e,double costheta);
float enu(float ee,double costheta);
void genedir(double emin, double emax, double *ep, double *cth, double *phi);
double dsigmasv_max(double ep);
//void getspec_(double *emin, double *emax, double *ux, double *uy, double *uz, double *en, double *upx, double *upy, double *upz);

/*************** Constants ***************************************************************/
double c = 2.99792458e10;
double pi=3.14159265358979323846264338327950288;
double distconv=0.1973269631;
double gf=1.16637e-5;
double tantheta=0.2257/0.97419;
double mp=938.272;
double mn=939.5654;
double me=0.511;
double delta=mn-mp;

/*Rotation matrix in 3D given axis*************************************************/
void rot(double ux, double uy, double uz, double theta, double r[3][3]){
    double ct = cos(theta);
    double st = sin(theta);
    r[0][0] = ct + pow(ux, 2) * (1 - ct);
    r[1][1] = ct + pow(uy, 2) * (1 - ct);
    r[2][2] = ct + pow(uz, 2) * (1 - ct);
    r[0][1] = ux * uy * (1 - ct) - uz * st;
    r[1][0] = ux * uy * (1 - ct) + uz * st;
    r[0][2] = ux * uz * (1 - ct) + uy * st;
    r[2][0] = ux * uz * (1 - ct) - uy * st;
    r[1][2] = uy * uz * (1 - ct) - ux * st;
    r[2][1] = uy * uz * (1 - ct) + ux * st;
}

/*IBD cross-section (Beacom-Vogel)*************************************************/
double dsigma(double e_nu,double e_e,double costheta)
{
  double pi=3.14159265358979323846264338327950288;
  double distconv=0.1973269631;
  double gf=1.16637e-5,tantheta=0.2257/0.97419;
  double mp=938.272013,mn=939.56556,me=0.510998910;
  double f=1,g=1.2695,f2=3.706;
  double delta_inner=0.024;
 
  double cos2theta=1/(1+tantheta*tantheta);
  double sigma0=gf*gf*cos2theta/pi*(1+delta_inner)*distconv*distconv*1e-32;
  double delta=mn-mp,m=(mn+mp)/2;
  double me2=me*me,delta2=delta*delta;
  double fbar=2*(f+f2);
  double k=g*(fbar+2*g),s=f*f+3*g*g,a=f*f-g*g;

  double sig,e0,v0,ve,efrac1,efrac2,efrac3;

  e0=e_nu-delta;
  if (e0<=me) return(0);
  v0=sqrt(1-me2/(e0*e0));
  e_e=e0*(1-e_nu*(1-v0*costheta)/m)-(delta2-me2)/m;
  ve=sqrt(1-me2/((e_e)*(e_e)));

  efrac1=e_nu/m;
  efrac2=2*efrac1-delta/m;
  efrac3=efrac1-delta/m;

  sig=(1-2*efrac1)*(ve*s+a*costheta)+k*efrac2*(costheta-v0)+
    v0*efrac1*(-a+2*ve*s*costheta+3*a*costheta*costheta);
  sig=sig*e0*e0+me2*((v0*(k-s)+2*ve*s)*efrac3-
		     (k*efrac2+a*(1-2*efrac3)-s*efrac1)*costheta);
  sig-=2*delta2*efrac3*(ve*s+a*costheta);
  return(sig*sigma0/2);
}

/*IBD cross-section (Strumia-Vissani)*************************************************/
double dsigma_sv(float enu, double costheta){
    double hbarc2 = 0.389379365e-21;
    double alpha = 1./128;
    double g10 = -1.270;
    double me = 0.511;
    double MV2 = 0.71e6;
    double MA2 = 1.0e6;
    double mpi = 139.0;
    double xi = 3.706;
    double delta = 1.293;
    double M = 938.9;
    double mn = M + delta/2;
    double mp = M - delta/2;
    double Gf = 1.1663787e-11;
    double cthc = 0.9746;
    double epsilon = enu/mp;
    double delta2 = (pow(mn, 2) - pow(mp, 2) - pow(me, 2))/(2 * mp);
    double kappa = pow(1 + epsilon, 2) - pow(epsilon * costheta, 2);
    double ee = ((enu - delta2) * (1 + epsilon) + epsilon * costheta * sqrt(pow(enu - delta, 2) - pow(me, 2) * kappa))/kappa;
    double t = pow(mn, 2) - pow(mp, 2) - 2 * mp * (enu - ee);
    double f1 = (1 - (1 + xi) * t/pow(2 * M, 2))/((1 - t/pow(2 * M, 2)) * pow(1 - t/MV2, 2));
    double f2 = xi/((1 - t/pow(2 * M, 2)) * pow(1 - t/MV2, 2));
    double g1 = g10/pow(1 - t/MA2, 2);
    double g2 = 2 * pow(M, 2) * g1/(pow(mpi, 2) - t);
    double A = (t - pow(me, 2)) * (4 * pow(f1, 2) * (4 * pow(M, 2) + t + pow(me, 2)) + 4 * pow(g1, 2) * (-4 * pow(M, 2) + t + pow(me, 2)) + pow(f2, 2) * (pow(t/M, 2) + 4 * t + 4 * pow(me, 2)) + 4 * pow(me, 2) * t * pow(g2/M, 2) + 8 * f1 * f2 * (2 * t + pow(me, 2)) + 16 * pow(me, 2) * g1 * g2) - pow(delta, 2) * ((pow(2 * f1, 2) + t * pow(f2/M, 2)) * (pow(2 * M, 2) + t - pow(me, 2)) + 4 * pow(g1, 2) * (pow(2 * M, 2) - t + pow(me, 2)) + pow(2 * me * g2, 2) * (t - pow(me, 2))/pow(M, 2) + 8 * f1 * f2 * (2 * t - pow(me, 2)) + pow(4 * me, 2) * g1 * g2) - 32 * pow(me, 2) * M * delta * g1 * (f1 + f2);
    A /= 16;
    double B = t * g1 * (f1 + f2) + 0.25 * pow(me, 2) * delta * (pow(f2, 2) + f1 * f2 + 2 * g1 * g2)/M;
    double C = 0.25 * (pow(f1, 2) + pow(g1, 2)) - 1./16 * t * pow(f2/M, 2);

    double smu = 2 * mp * (enu + ee) - pow(me, 2);
    double smp2 = 2 * mp * enu;
    double m2 = A - smu * B + pow(smu, 2) * C;
    double pe = sqrt(pow(ee, 2) - pow(me, 2));
    double fact = pe * epsilon/(1 + epsilon * (1 - ee/pe * costheta));
    double rad = 1 + alpha/pi * (6.0 + 3./2 * log(mp/(2 * ee)) + 1.2 * pow(me/ee, 1.5));
    double dsigmadee = 2 * mp * hbarc2 * pow(Gf * cthc/smp2, 2)/(2 * pi) * m2 * fact * rad;
    return dsigmadee;
}

/*Max IBD cross-section (SV) for random number generation***********/
double dsigmasv_max(double ep){
    double th_step = 100;
    double smax = -1;
    for(int i = -th_step; i <= th_step; i += 1){
        double cth = i * 1./th_step;
        double snew = dsigma_sv(enu(ep, cth), cth);
        if (snew > smax) smax = snew;
    }
    return smax;
}

/*Draw random positron energy, theta, and phi for given energy bin ********/
void genedir(double emin, double emax, float *ep, double *cth, double *phi){
    *ep = random()/((double) RAND_MAX) * (emax - emin) + emin;
    double smax = dsigmasv_max(*ep);
    *cth = random()/((double) RAND_MAX) * 2 - 1;
    while(dsigma_sv(enu(*ep,*cth), *cth)/smax < random()/((double) RAND_MAX)){
        *cth = random()/((double) RAND_MAX) * 2 - 1;
    }
    *phi = random()/((double) RAND_MAX) * 2 * pi;
}

/*Enu from positron energy*************************************************/
float enu(float ee,double costheta)
{
  double mp=938.272013,mn=939.56556,me=0.510998910;
  double delta=mn-mp,d=(delta*delta-me*me)/(2*mp);
  double pe=sqrt(ee*ee-me*me);
  
  return((ee+delta+d)/(1-(ee-pe*costheta)/mp));
}

/*Get positron direction and neutrino energy given positron energy and nu direction *******/
extern "C"{
    void getspec_(float *emin, float *emax, float *ux, float *uy, float *uz, float *ep, float *en, float *upx, float *upy, float *upz){
        double cth, phi;
        genedir(*emin,*emax,ep,&cth,&phi);
        double r[3][3];
        // Rotate neutrino direction by theta
        // Rotation axis is anything orthogonal to nu direction
        double axis1[] = {-*uy/sqrt(pow(*ux, 2) + pow(*uy, 2)), *ux/sqrt(pow(*ux, 2) + pow(*uy, 2)), 0};
        if (fabs(*ux) < 1e-10 && fabs(*uy) < 1e-10){
            axis1[0] = *uz/sqrt(pow(*ux, 2) + pow(*uz, 2));
            axis1[1] = 0;
            axis1[2] = -*ux/sqrt(pow(*ux, 2) + pow(*uz, 2));
        }
        rot(axis1[0], axis1[1], axis1[2], acos(cth), r);
        float uux, uuy, uuz;
        uux = r[0][0] * (*ux) + r[0][1] * (*uy) + r[0][2] * (*uz);
        uuy = r[1][0] * (*ux) + r[1][1] * (*uy) + r[1][2] * (*uz);
        uuz = r[2][0] * (*ux) + r[2][1] * (*uy) + r[2][2] * (*uz);
        // Rotate again by phi around initial neutrino direction
        rot(*ux, *uy, *uz, phi, r);
        *upx = r[0][0] * (uux) + r[0][1] * (uuy) + r[0][2] * (uuz);
        *upy = r[1][0] * (uux) + r[1][1] * (uuy) + r[1][2] * (uuz);
        *upz = r[2][0] * (uux) + r[2][1] * (uuy) + r[2][2] * (uuz);
        // Neutrino energy
        *en = enu(*ep, cth);
    }
}
