#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "mcinfo.h"
#include "TMath.h"
#include "dsigma.h"

class MCInfo;

using namespace std;

// Read spectrum from text file
// Get bin width, lowest nu energy
// Also return integral of spectrum (trapezoidal method)
double read_spectrum(char *fname, float &low, float &bwidth, int &nbins, double *spectrum){
    FILE *file = fopen(fname, "r");
    char sdum[1000];
    double x[2500];
    fscanf(file, "%[^\n]\n", sdum);
    int count = 0;
    double integral = 0;
    while(!feof(file)){
        fscanf(file, "%lf %lf\n", &(x[count]), &(spectrum[count]));
        if (count > 0) integral += 0.5 * (spectrum[count - 1] + spectrum[count]) * (x[count] - x[count - 1]);
        count++;
    }
    nbins = count;
    low = x[0];
    bwidth = x[1] - x[0];
    fclose(file);
    return integral;
}

// Normalize spectrum for positrons
int normspec(double *spectrum, int nbins, double norm){
    for(int i = 0; i < nbins; i++){
        spectrum[i] /= norm;
    }
    return 0;
}

int main(int argc, char **argv){
    if (argc < 5){
        cout << "Usage ./reweight <0 for antinu, 1 for positron, 2 for Enu interacted> <output file> <spectrum file> <input file> [livetime (days)]" << endl;
    }
    double spectrum[2500];
    float low, bwidth;
    int is_positron = atoi(argv[1]);
    double livetime = 2790.1;
    if (argc > 5)
        livetime = atof(argv[5]);
    cout << "Livetime is " << livetime << " days" << endl; 
    int nbins;
    double norm = read_spectrum(argv[3], low, bwidth, nbins, spectrum);
    if (is_positron == 1)
        normspec(spectrum, nbins, norm);

    MCInfo *mc = new MCInfo;

    TFile *f = new TFile(argv[4]);
    TTree *tr = (TTree*) f->Get("data");
    TFile *fout = new TFile(argv[2], "RECREATE");
    Float_t wt;
    Float_t enu;
    TBranch *br = (TBranch*) tr->Branch("weight", &wt, "weight/F");
    tr->SetBranchAddress("MC", &mc);
    TTree *newtree = tr->CloneTree(0);
    int nentries = tr->GetEntries();
    for(int i = 0; i < nentries; i++){
        if (i % 1000 == 0) cout << i << "/" << nentries << endl;
        tr->GetEntry(i);
        float pe2 = mc->pvc[1][0] * mc->pvc[1][0] + mc->pvc[1][1] * mc->pvc[1][1] + mc->pvc[1][2] * mc->pvc[1][2];
        // Call weight function
        if (is_positron == 1){
            float ep = sqrt(pe2 + 0.511 * 0.511);
            wt = weight_ep(ep, low, bwidth, nbins, spectrum, 1, 1);
        }
        else if (is_positron == 0){
            float enu = sqrt(mc->pvc[0][0] * mc->pvc[0][0] + mc->pvc[0][1] * mc->pvc[0][1] + mc->pvc[0][2] * mc->pvc[0][2]);
            float ctheta = (mc->pvc[1][0] * mc->pvc[0][0] + mc->pvc[1][1] * mc->pvc[0][1] + mc->pvc[1][2] * mc->pvc[0][2])/(enu * sqrt(pe2));
            wt = weight_enu(enu, ctheta, low, bwidth, nbins, spectrum, livetime);
        }
        else{
            float enu = sqrt(mc->pvc[0][0] * mc->pvc[0][0] + mc->pvc[0][1] * mc->pvc[0][1] + mc->pvc[0][2] * mc->pvc[0][2]);
            wt = weight_enu_interacted(enu, low, bwidth, nbins, spectrum);
        }
        newtree->Fill();
    }
    fout->Write();
    fout->Close();
}
