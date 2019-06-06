#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <string>
#include <vector>
#include <map>
//#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>

using namespace std;
//namespace py = pybind11;

/****************************************************************************/
/*          This code contains useful functions to generate the             */
/*    DSNB spectrum as well as the IBD spectrum in SK (22.5kton.year)       */
/*    It is used to generate a library than can then be imported as a       */
/*    Python module but it can also be used as a standalone if a main       */
/*    method is added.                                                      */
/****************************************************************************/

/**Prototypes ***************************************************************/
double dsigma(double e_nu,double e_e,double costheta);
double enu(double ee,double costheta);
double integrated_imf(double mmin, double mmax, double xi1, double xi2);
double integrated_mimf(double mmin, double mmax, double xi1, double xi2);
double csfr(double z, double rho0, double alpha, double beta, double gamma, double z1, double z2);
double spectrum(double enu, map<string,double>);
double cosmology(double z);
double snrate_integrand(double z, double enu, map<string,double>, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2);
double snrate(double enu, map<string,double>, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2);
double skrate_integrand(double ctheta, double ee, map<string,double>, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2);
double skrate(double ee, map<string,double>, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2);
//extern "C" {void srnspectrum_(double *ee)};

/*************** Constants ***************************************************************/
double c = 2.99792458e10;
double secondsyear = 3600 * 24 * 365;
//double secondsyear = 3600 * 24 * 365.25;
double cmMpc = 3.08568025e24;
//double cmMpc = 3.26e6 * c * secondsyear;
double GeVtoErg = 0.0016;
double mcut = 0.5;
double min_m = 0.1;
double eta = -10;
double H0 = 70.90644 * (1e5/cmMpc);
double OmegaM = 0.3;
double OmegaL = 0.7;
//double Enutot = 3e53 * 1e3/GeVtoErg;
double mSun = 2e33;
double pi=3.14159265358979323846264338327950288;
double distconv=0.1973269631;
double gf=1.16637e-5;
double tantheta=0.2257/0.97419;
double mp=938.272013;
double mn=939.56556;
double me=0.510998910;
double delta=mn-mp;

/*IMF integrated ***************************************************************/
double integrated_imf(double mmin, double mmax, double xi1 = 2.35, double xi2 = 2.35){
    if (mmin > mcut){ 
        if (xi2 != 1)
            return (pow(mmax, 1 - xi2) - pow(mmin, 1 - xi2))/(1 - xi2);
        else
            return log(mmax/mmin);
    }
    if (mmax < mcut){ 
        if (xi1 != 1)
            return (pow(mmax, 1 - xi1) - pow(mmin, 1 - xi1))/(1 - xi1);
        else
            return log(mmax/mmin);
    }
    double up = xi2 != 1 ? (pow(mmax, 1 - xi2) - pow(mcut, 1 - xi2))/(1 - xi2) : log(mmax/mcut);
    double down = xi1 != 1 ? (pow(mcut, 1 - xi1) - pow(mmin, 1 - xi1))/(1 - xi1) : log(mcut/mmin);
    return  up + down;
}

/*M x IMF integrated ***************************************************************/
double integrated_mimf(double mmin, double mmax, double xi1 = 2.35, double xi2 = 2.35){
    if (mmin > mcut){ 
        if (xi2 != 2)
            return (pow(mmax, 2 - xi2) - pow(mmin, 2 - xi2))/(2 - xi2);
        else
            return log(mmax/mmin);
    }
    if (mmax < mcut){
        if (xi1 != 2)
            return (pow(mmax, 2 - xi1) - pow(mmin, 2 - xi1))/(2 - xi1);
        else
            return log(mmax/mmin);
    }
    double up = xi2 != 2 ? (pow(mmax, 2 - xi2) - pow(mcut, 2 - xi2))/(2 - xi2) : log(mmax/mcut);
    double down = xi1 != 2 ? (pow(mcut, 2 - xi1) - pow(mmin, 2 - xi1))/(2 - xi1) : log(mcut/mmin);
    return  up + down;
}

/*Cosmic star formation rate (from Horiuchi-Beacom-Dweck 2008 paper)************************/
double csfr(double z, double rho0 = 0.0178, double alpha = 3.4, double beta = -0.5, double gamma = -4.5, double z1 = 1, double z2 = 4){
    double CSFRtoErg = 1/(secondsyear * pow(cmMpc, 3));
    double B = pow(1 + z1, 1 - alpha/beta);
    double C = pow(1 + z1, (beta - alpha)/gamma) * pow(1 + z2, 1 - beta/gamma);
    return CSFRtoErg * rho0 * exp(log(pow(1 + z, alpha * eta) + pow((1 + z)/B, beta * eta) + pow((1 + z)/C, gamma * eta))/eta);
}

/*Pinched Fermi-Dirac spectrum*******************************************************/
double pinched(double enu, double emean, double alpha){
    return pow(1 + alpha, 1 + alpha)/(pow(emean, 2) * tgamma(alpha + 1)) * pow(enu/emean, alpha) * exp(-(1 + alpha) * enu/emean);
}

/*Neutrino spectrum ********************************************/
double spectrum(double enu, map<string,double> spec_params){
    // Blackbody
    if (spec_params["type"] == 0){
        double tnu = spec_params["tnu"];
        double Enutot = spec_params["lumi"];
        return Enutot/6 * 120/(7 * pow(pi, 4)) * pow(enu, 2)/pow(tnu, 4) * 1/(exp(enu/tnu) + 1);
    }
    //Pinched (no flavors, no BH)
    if (spec_params["type"] == 1){
        double emean = spec_params["emean"];
        double alpha = spec_params["alpha"];
        double Enutot = spec_params["lumi"];
        return Enutot * pinched(enu, emean, alpha);
    }
    //Pinched (flavors, no BH)
    if (spec_params["type"] == 2){
        double pbar = spec_params["pbar"];
        double Enutot = spec_params["lumi"];
        double emeanE = spec_params["emeanE"];
        double alphaE = spec_params["alphaE"];
        double emeanX = spec_params["emeanX"];
        double alphaX = spec_params["alphaX"];
        double lumX = spec_params["lumX"];
        return Enutot * (pbar * pinched(enu, emeanE, alphaE) + (1 - pbar) * lumX * pinched(enu, emeanX, alphaX));
    }
    //Pinched (flavors, BH)
    if (spec_params["type"] == 3){
        double pbar = spec_params["pbar"];
        double Enutot = spec_params["lumi"];
        double fBH = spec_params["fBH"];
        double emeanE_NS = spec_params["emeanE_NS"];
        double emeanE_BH = spec_params["emeanE_BH"];
        double alphaE = spec_params["alphaE"];
        double emeanX_NS = spec_params["emeanX_NS"];
        double emeanX_BH = spec_params["emeanX_BH"];
        double alphaX = spec_params["alphaX"];
        double lumE_BH = spec_params["lumE_BH"];
        double lumX_NS = spec_params["lumX_NS"];
        double lumX_BH = spec_params["lumX_BH"];
        //cout << pinched(enu, emeanE_NS, alphaE) << " " << pinched(enu, emeanE_BH, alphaE) << endl;
        return Enutot * (pbar * ((1 - fBH) * pinched(enu, emeanE_NS, alphaE) + fBH * lumE_BH * pinched(enu, emeanE_BH, alphaE)) + (1 - pbar) * ((1 - fBH) * lumX_NS * pinched(enu, emeanX_NS, alphaX) + fBH * lumX_BH * pinched(enu, emeanX_BH, alphaX)));
    }
    return 0;
}

/*Cosmological factor*************************************************/
double cosmology(double z){
    return 1/sqrt(OmegaM * pow(1 + z, 3) + OmegaL);
}

/*DSNB integrand (function of redshift)*************************************************/
double snrate_integrand(double z, double enu, map<string,double> spec_params, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2){
    return csfr(z, rho0, alpha, beta, gamma, z1, z2) * integrated_imf(8, 50, xi1, xi2)/integrated_mimf(0.1, 100, xi1, xi2) * spectrum(enu * (1 + z), spec_params) * cosmology(z);
}

/*DSNB spectrum (function of Enu)*************************************************/
double snrate(double enu, map<string,double> spec_params, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2){
    double z = 0;
    double dz = 0.2;
    double res = snrate_integrand(z, enu, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2)/2;
    z += dz/2;
    while(z < 5){
        res += 2 * snrate_integrand(z, enu, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2);
        res += snrate_integrand(z + dz/2, enu, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2);
        z += dz;
    }
    res += snrate_integrand(z, enu, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2)/2;
    res *= dz/3 * c/H0;
    return res;
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

/*Enu from positron energy*************************************************/
double enu(double ee,double costheta)
{
  double mp=938.272013,mn=939.56556,me=0.510998910;
  double delta=mn-mp,d=(delta*delta-me*me)/(2*mp);
  double pe=sqrt(ee*ee-me*me);
  
  return((ee+delta+d)/(1-(ee-pe*costheta)/mp));
}

/*Integrand SK DSNB rate (function of costheta)***************************************/
double skrate_integrand(double ctheta, double ee, map<string,double> spec_params, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2){
    double en = enu(ee, ctheta);
    return snrate(en, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2) * dsigma(en, 0, ctheta) * 0.047390723e42;
}

double random2(){
    return rand()/((double) RAND_MAX);
}

extern "C" {
    void srnspectrum_(double &ep, double &en, int *sd, char *finp){
        // Reformat file name properly because fortran sucks
        char finput[200];
        int cr = 0;
        while(!isspace(finp[cr])){
            finput[cr] = finp[cr];
            cr++;
        }
        finput[cr] = '\0';
        // Open spectral parameter file and read it
        FILE *fin = fopen(finput, "r");
        map<string,double> spec_params;
        map<string,double> snrate_params;
        string snrate_k[] = {"rho0", "alpha", "beta", "gamma", "z1", "z2", "xi1", "xi2"};
        vector<string> snrate_keys(snrate_k, snrate_k + 8);
        char key_char[500];
        double param;
        while(!feof(fin)){
            fscanf(fin, "%s %lf\n", &key_char, &param);
            const string key(key_char);
            if (find(snrate_keys.begin(), snrate_keys.end(), key) != snrate_keys.end())
                snrate_params.insert(pair<string,double>(key, param));
            else
                spec_params.insert(pair<string,double>(key, param));
        }
        fclose(fin);
        srand(*sd); // Set random seed using time
        //double rho0 = 0.0213;
        //double alpha = 3.6;
        //double beta = -0.3;
        //double gamma = -2.5;
        //double z1 = 1;
        //double z2 = 1;
        //double xi1 = 2.35;
        //double xi2 = 2.35;
        //spec_params.insert(pair<string,double>("lumi", 3e53 * 1e3/0.0016));
        //spec_params.insert(pair<string,double>("type", 0));
        //spec_params.insert(pair<string,double>("tnu", 6));
        double ctheta = 2 * (random2() - 0.5); // cos(theta) random
        // Find maximum flux
        double fmax = 0;
        for (int i = 0; i < 300; i++){
            double fl = skrate_integrand(ctheta, i * 0.1, spec_params, snrate_params["rho0"], snrate_params["alpha"], snrate_params["beta"], snrate_params["gamma"], snrate_params["z1"], snrate_params["z2"], snrate_params["xi1"], snrate_params["xi2"]);
            if (fl > fmax) fmax = fl;
        }
        // Pick random positron energy using rejection method
        while(true){
            ep = random2() * 30;
            double fl = skrate_integrand(ctheta, ep, spec_params, snrate_params["rho0"], snrate_params["alpha"], snrate_params["beta"], snrate_params["gamma"], snrate_params["z1"], snrate_params["z2"], snrate_params["xi1"], snrate_params["xi2"]);
            //cout << fl << " " << fmax << random2() << endl;
            //cout << ep << endl;
            if (random2() < fl/fmax) break;
        }
        en = enu(ep, ctheta);
    }
}

/*SK DSNB rate (function of positron energy)*********************************************/
double skrate(double ee, map<string,double> spec_params, double rho0, double alpha, double beta, double gamma, double z1, double z2, double xi1, double xi2){
    // This method is slower than the one used in scipy quad
    double z = -1;
    double dz = 0.01;
    double res = skrate_integrand(z, ee, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2)/2;
    z += dz/2;
    while(z < 1){
        res += 2 * skrate_integrand(z, ee, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2);
        res += skrate_integrand(z + dz/2, ee, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2);
        z += dz;
    }
    res += skrate_integrand(z, ee, spec_params, rho0, alpha, beta, gamma, z1, z2, xi1, xi2)/2;
    res *= dz/3;
    return res;
}

//int main(){
    //double ee;
    //srnspectrum_(&ee);
    //cout << ee << endl;
//}

//int main(){
    //double x = 0;
    //double tstart = time(NULL);
    //while (x <= 40){
        ////printf("%f %e\n", x, skrate(30, 6, 0.0213, 3.6, -0.1, -2.5, 1, 4, 2.35, 2.35));
        //double s = skrate(30, 6, 0.0213, 3.6, -0.1, -2.5, 1, 4, 2.35, 2.35);
        //x += 0.1;
    //}
    //cout << time(NULL) - tstart << endl;
//}

/*pybind11 bindings************************************************/

//PYBIND11_MODULE(snrate, m){
    //m.doc() = "pybind11 plugin for DSNB spectra";
    //m.def("dsigma", &dsigma, "Beacom-Vogel cross-section");
    //m.def("enu", &enu, "Neutrino energy");
    //m.def("cosmology", &cosmology, "Cosmological factor");
    //m.def("spectrum", &spectrum, "Neutrino spectrum");
    //m.def("integrated_imf", &integrated_imf, "Integrated IMF", py::arg("mmin"), py::arg("mmax"), py::arg("xi1"), py::arg("xi2"));
    //m.def("integrated_mimf", &integrated_mimf, "Integrated M x IMF", py::arg("mmin"), py::arg("mmax"), py::arg("xi1"), py::arg("xi2"));
    //m.def("csfr", &csfr, "Cosmic star formation rate", py::arg("z"), py::arg("rho0"), py::arg("alpha"), py::arg("beta"), py::arg("gamma"), py::arg("z1"), py::arg("z2"));
    //m.def("skrate", &skrate, "SK DSNB rate as a function of positron energy");
    //m.def("snrate", &snrate, "DSNB spectrum");
//}
