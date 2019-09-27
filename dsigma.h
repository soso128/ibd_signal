#ifndef DSIGMA_H
#define DSIGMA_H
double dsigma(double e_nu,double e_e,double costheta);
float enu(float ee,double costheta);
void genedir(double emin, double emax, double *ep, double *cth, double *phi);
double dsigmasv_max(double ep);
float weight(float enu, float low, float bwidth, double *spectrum);
//void getspec_(double *emin, double *emax, double *ux, double *uy, double *uz, double *en, double *upx, double *upy, double *upz);
#endif
