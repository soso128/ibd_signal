#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "mcinfo.h"
#include "TMath.h"
// Include function that computes weights
// in parent directory
#include "dsigma.h"

class MCInfo;

using namespace std;

// Read spectrum from text file
// Get bin width, lowest nu energy
int read_spectrum(char *fname, float &low, float &bwidth, double *spectrum){
    FILE *file = fopen(fname, "r");
    char sdum[1000];
    double x[1000];
    fscanf(file, "%[^\n]\n", sdum);
    int count = 0;
    while(!feof(file)){
        fscanf(file, "%lf %lf\n", &(x[count]), &(spectrum[count]));
        count++;
    }
    low = x[0];
    bwidth = x[1] - x[0];
    fclose(file);
    return 0;
}

int main(int argc, char **argv){
    if (argc < 3){
        cout << "Usage ./reweight <output file> <spectrum file> <input files>" << endl;
    }
    double spectrum[1000];
    float low, bwidth;
    read_spectrum(argv[2], low, bwidth, spectrum);

    MCInfo *mc = new MCInfo;

    TFile *f = new TFile(argv[3]);
    TTree *tr = (TTree*) f->Get("data");
    TFile *fout = new TFile(argv[1], "RECREATE");
    Float_t wt;
    Float_t enu;
    TBranch *br = (TBranch*) tr->Branch("weight", &wt, "weight/F");
    tr->SetBranchAddress("MC", &mc);
    TTree *newtree = tr->CloneTree(0);
    int nentries = tr->GetEntries();
    for(int i = 0; i < nentries; i++){
        tr->GetEntry(i);
        float enu = sqrt(mc->pvc[0][0] * mc->pvc[0][0] + mc->pvc[0][1] * mc->pvc[0][1] + mc->pvc[0][2] * mc->pvc[0][2]);
        // Call weight function
        wt = weight(enu, low, bwidth, spectrum);
        newtree->Fill();
    }
    fout->Write();
    fout->Close();
}
