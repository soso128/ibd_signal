#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include "TIterator.h"
#include "tqrealroot.h"
#include <algorithm>

using namespace std;
class TQReal;


void count(){
    ofstream of;
    of.open("run_stat_t2k_akutsu.txt");
    for (int k = 0; k < 9; k++){
        TString dirname(Form("/disk01/usr5/rakutsu/t2k/Neutron/work/dummy/out/inbtwRun%d_%d/", k, k+1));
        TSystemDirectory tdir(dirname, dirname);
        TList *files = tdir.GetListOfFiles();
        TSystemFile *file;
        TString fname;
        TIter next(files);
        TQReal *tq;
        tq = new TQReal;

        while((file = (TSystemFile*) next())){
            fname = file->GetName();
            if (file->IsDirectory() || !(fname.EndsWith(".root"))) continue;
            TFile *f = TFile::Open(fname);
            TTree *tr = (TTree*) f->Get("data");
            int nentries = tr->GetEntries();
            int nh = 0;
            TBranch *br = tr->GetBranch("TQREAL");
            br->SetAddress(&tq);
            int ntot = min(nentries, 100);
            for(int i = 0; i < ntot; i++){
                tr->GetEntry(i);
                nh += tq->nhits;
            }
            nh /= ntot;
            TObjArray *table = fname.Tokenize(".");
            TString run = ((TObjString*) table->At(1))->String();
            printf("%s %d %d\n", run.Data(), nentries, nh);
            of << run << " " << nentries << " " << nh << "\n";
            f->Close();
        }
    }
    of.close();
}

int main(){
    count();
}
