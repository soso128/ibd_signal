#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TBranch.h"

#include <iostream>
#include <fstream>
#include <string>
//#include <vector>
using namespace std;

#include "skheadC.h"
#include "skparmC.h"
#include "sktqC.h"
#include "skbadcC.h"

#include "tqrealroot.h"  // include HEADER information
#include "loweroot.h"
#include "mcinfo.h"
#include "stdio.h" // check whether file exist
#include "stdlib.h" 
#include <iomanip> // set io
#include "skroot.h"
#include <algorithm>

const Float_t mintime = 18000; //ns
const Float_t maxtime = 535000;
const Float_t deadtime = 900;
const Int_t NEVENTGAP = 1500; //~2000 events per MC file, 2 MC/1 dummy, so 1500 is enough
Int_t twind = 517000;
Int_t tbeam = 500000;
Int_t NHITUP = 65000;
Int_t NHITLOW = 50000;
Float_t CENTERT = 500000;

int   access(const   char   *filename,   int   amode);

extern "C" {
	void initialize_(void);
}
//
int main(int argc, char **argv)
{
	if (argc != 4)
	{
		cout << "Usage: " << argv[0] << " fname_out fname_in nrun" << endl;
		exit(1);
	}

	TString fname_out = argv[1];
	TString fname_in = argv[2];
	TString nrun_char = argv[3];
        int nfile = nrun_char.Atoi();
	cout << "Input file: " << fname_in << endl;
	cout << "Output file: " << fname_out << endl;
	cout << "Run number: " << nfile << endl;

        // Find run time bin and the corresponding T2K runs to use
        // Info in a text file
        FILE *fruns = fopen("/home/elhedri/SK2p2MeV/mc/generate/ibd_signal/run_info_t2k.txt", "r");
        int timebin = -1;
        map<int,int> runs;
        int run, bin, dum;
        int startbin = -1, oldbin = -1;
        int count = 0;
        // Fill run table and find time bin for current run
        while(!feof(fruns)){
            fscanf(fruns, "%d %d %d\n", &run, &bin, &dum);
            if (bin != oldbin){
                if (timebin != -1) break;
                startbin = count;
                oldbin = bin;
            }
            if (run == nfile){
                timebin = bin;
            }
            pair<int,int> p(run,dum);
            runs.insert(p);
            count++;
        }
        // Remove everything not in the right time bin
        // (for the first bin, nothing to remove)
        // So at the end we have a hash table that, for a given run, tells us
        // whether it is dummy or real bin
        if (startbin > 0)
            for(int q = 0; q < startbin; q++)
                runs.erase(runs.begin());
        fclose(fruns);
        cout << "Time bin is " << timebin << endl;

        // Load T2K data for the run period we want
	TChain *t2kch = new TChain("data");
	//dummy trigger chain
        for(map<int,int>::iterator iter = runs.begin(); iter != runs.end(); iter++){
            int irun = iter->first;
            char inname[500];
            for (int k=0; k < 9; k++){
                sprintf(inname, "/disk01/usr5/rakutsu/t2k/Neutron/work/dummy/out/inbtwRun%d_%d/t2k.0%d*.root", k, k+1, irun);
                t2kch->Add(inname);
            }
        }
        int nt2k = t2kch->GetEntries();
        cout << nt2k << " t2k entries available for run " << nfile << " in bin " << timebin << endl;
	//cout << "Processing run " << n << endl;
        if (t2kch->GetEntries() == 0){
            cout << "no T2K data for time period " << timebin << endl;
            exit(42);
        }

	TQReal *TQI = new TQReal;
	TQReal *TQA = new TQReal;
	Header *HEADER = new Header;

	t2kch->SetBranchAddress("TQREAL", &TQI);
	t2kch->SetBranchAddress("TQAREAL", &TQA);
	t2kch->SetBranchAddress("HEADER", &HEADER);

	// Input file
	//
	Float_t toffset;
	Int_t entry; // 2 MC event ~ 1 dummy trigger event
	cout << "Combining MC..." << endl;
	//open input sample file (with reconstructions/after skdetsim, etc)
	//
	// Input file
	TFile *fileIn = new TFile(fname_in.Data());
	TTree *treeIn = (TTree*) fileIn->Get("data");
	cout << "Processing " << fname_in.Data() << endl;

	// Get HEADER and LOWE, MU branches. LOWE & MU are for copying to output file
	Header *lomu_head = new Header();
	LoweInfo *lomu_lowe = new LoweInfo();
	TQReal *lomu_tqi = new TQReal; // TQ (in-gate)
	TQReal *lomu_tqa = new TQReal; // Anti-TQ (in-gate)
	MuInfo *lomu_mu = new MuInfo();
	MCInfo *lomu_mc = new MCInfo();
	TBranch *br_hd = treeIn->GetBranch("HEADER");
	TBranch *br_le = treeIn->GetBranch("LOWE");
	TBranch *br_mu = treeIn->GetBranch("MU");
	TBranch *br_mc = treeIn->GetBranch("MC");
	TBranch *br_tqi = treeIn->GetBranch("TQREAL");
	TBranch *br_tqa = treeIn->GetBranch("TQAREAL");
	br_hd->SetAddress(&lomu_head);
	br_le->SetAddress(&lomu_lowe);
	br_mu->SetAddress(&lomu_mu);
	br_mc->SetAddress(&lomu_mc);
	br_tqi->SetAddress(&lomu_tqi);
	br_tqa->SetAddress(&lomu_tqa);

	//open output file
	TFile *fileOut = new TFile(fname_out.Data(), "RECREATE");
	cout << "Output to " << fname_out.Data() << endl;
	TTree *treeOut = new TTree("data", "SK-IV 2p2-dummy data");
	Int_t bufsize = 2*1024*1024;
	Header *HEAD = new Header;
	TQReal *TQII = new TQReal; // TQ (in-gate)
	TQReal *TQAA = new TQReal; // Anti-TQ (in-gate)
	MCInfo *MC = new MCInfo;
        Int_t is_signal[100000];
	treeOut->Branch("issignal",  is_signal, "issignal[100000]/I");
	treeOut->Branch("HEADER",  "Header", &lomu_head, bufsize, 0);
	treeOut->Branch("TQREAL",  "TQReal", &TQII,  bufsize, 0);
	treeOut->Branch("TQAREAL", "TQReal", &TQAA,  bufsize, 0);
	treeOut->Branch("LOWE", "LoweInfo", &lomu_lowe,  bufsize, 0);
	treeOut->Branch("MC", "MCInfo", &lomu_mc,  bufsize, 0);

	// Initialize for root 
	initialize_();
        
        //SKTQZ
        Float_t prevt[MAXPM];
        Int_t ihtiflz[30*MAXPM], icabiz[30*MAXPM], itiskz[30*MAXPM], iqiskz[30*MAXPM], index[30*MAXPM];
        Float_t tiskz[30*MAXPM], qiskz[30*MAXPM];
        //SKTQAZ	
        Int_t ihtflz[30*MAXPMA], icabaz[30*MAXPMA], itaskz[30*MAXPMA], iqaskz[30*MAXPMA] ; 
        Int_t ihacab[30*MAXPMA], indexa[30*MAXPMA];
        Float_t taskz[30*MAXPMA], qaskz[30*MAXPMA]; 
        Int_t nqiskz_save = 0, nhitaz_save = 0;

	Int_t nrun, nsub, nev;
	cout << "Reading event loop: " << endl;
	//loop over events
	Int_t nentries = treeIn->GetEntries();
	cout << "Get 2p2 candidate entries " <<  treeIn->GetEntries() << endl;
	for (Int_t I2P2=0; I2P2<nentries; I2P2 ++) {
            cout << "Event " << I2P2 << endl;
		treeIn->GetEntry(I2P2);
                // Pick a T2K event at random
                int is_dummy, t2krun;
                do{
                    int entry = (int) (random()/((double) RAND_MAX) * nt2k);
                    //read t2k
                    t2kch->GetEntry(entry);
                    t2krun = HEADER->nrunsk;
                    is_dummy = runs[t2krun];
                } while (is_dummy && (TQI->nhits > NHITUP || TQI->nhits < NHITLOW));
		//dummy trigger timing goes from -500us < t < 500us 
		//get the time offset at the beginnig of time window and 
		//the center of the time window 
                if (is_dummy < 0){
                    cout << "Error: there should not have been t2k data for run " << t2krun << endl;
                    is_dummy = 1; // temporary fix since we use only Akutsu-san's data for now
                }
                // If dummy event, we can take either the first or the second half of the event
                // For a real beam event we always have to take the first 500 microseconds
                int noffset = is_dummy ? (int) (random()/((double) RAND_MAX) * 2) : 0;
		if (noffset == 0)
			toffset = TQI->T[0];
		else
			toffset = TQI->T[TQI->nhits-1] - twind ;

                // Save info from SKDetSim IBD event
		for (Int_t i = 0; i < lomu_tqi->nhits; i++)
		{
			ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,0); //left shift 16 bit
			icabiz[nqiskz_save] = lomu_tqi->cables[i]&0x0000FFFF;
			tiskz[nqiskz_save] = lomu_tqi->T[i];
			qiskz[nqiskz_save] = lomu_tqi->Q[i];
			nqiskz_save++;
		}

                //copy from dummy trigger TQ data into temporary TQ banks.
                for (Int_t i = 0; i < TQI->nhits; i ++)
                {
                        //twind is how much of the dummy trigger banks we copy.
                        //offset is how much the first hit that we will copy is offset from 0.
                        float tmax = is_dummy ? toffset + twind : tbeam;
                        if (TQI->T[i] >= toffset && TQI->T[i] <= tmax)
                        {
                                tiskz[nqiskz_save] = TQI->T[i]-toffset +mintime;
                                qiskz[nqiskz_save] = TQI->Q[i];
                                icabiz[nqiskz_save] = TQI->cables[i]&0x0000FFFF;
                                //had problems if I just directly copy ihtiflz, so just set every 
                                //hit to the same default value.
                                ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,0);
                                nqiskz_save++;
                        }
                        // For T2K beam event, we need extra events (time window > 500 microseconds)
                        // We take them from the event afterwards (or from the event before if current event is the last one of the batch)
                }
                if( !is_dummy ){
                    int extra_event = entry < nt2k - 1 ? entry + 1 : entry - 1;
                    t2kch->GetEntry(extra_event);
                    for (Int_t i = 0; i < TQI->nhits; i ++){
                        if (TQI->T[i] <= TQI->T[0] + twind - tbeam + toffset)
                        {
                                tiskz[nqiskz_save] = TQI->T[i]-TQI->T[0]  - toffset + tbeam + mintime;
                                qiskz[nqiskz_save] = TQI->Q[i];
                                icabiz[nqiskz_save] = TQI->cables[i]&0x0000FFFF;
                                //had problems if I just directly copy ihtiflz, so just set every 
                                //hit to the same default value.
                                ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,0);
                                nqiskz_save++;
                        }
                        else break;
                    }
                }
                //there are some neutron hits saved which are out 
                //of order from the rest of the dummy trigger hits. I sort to fix it.
                Float_t tiskz2[nqiskz_save], qiskz2[nqiskz_save], icabiz2[nqiskz_save]; 
                //Float_t taskz2[nhitaz_save], qaskz2[nhitaz_save],icabaz2[nhitaz_save]; 
                for(int i=0;i<nqiskz_save;i++){
                        tiskz2[i] = tiskz[i];
                        qiskz2[i] = qiskz[i];
                        icabiz2[i] = icabiz[i];
                }
                TMath::Sort(nqiskz_save, tiskz2, index, kFALSE); //increasing order
                for(int i=0;i<nqiskz_save;i++){
                        tiskz[i] = tiskz2[index[i]];
                        qiskz[i] = qiskz2[index[i]];
                        icabiz[i] = icabiz2[index[i]];
                        if (index[i] < lomu_tqi->nhits)
                            is_signal[i] = 1;
                        else
                            is_signal[i] = 0;
                }
                //If a hit is okay, save it in sktqz_ timing bank
                // Set TQI 
                TQII->Clear();
                Int_t newnqiskz = 0;
                for (Int_t i = 0; i < MAXPM; i++)
                {
                        prevt[i] = -9999999;
                }
                for (Int_t i = 0; i<nqiskz_save;i++)
                {
                        //Linyan
                        if (tiskz[i] - prevt[icabiz[i]] < deadtime) continue;//cutting for pmt dead time.
                        if (prevt[icabiz[i]] < tiskz[i]) prevt[icabiz[i]] = tiskz[i];

                        if (tiskz[i] > maxtime) continue;

                        TQII->T.push_back(tiskz[i]);
                        TQII->Q.push_back(qiskz[i]);
                        TQII->cables.push_back(icabiz[i] + (ihtiflz[i] << 16));
                        newnqiskz++;
                }
                TQII->nhits = newnqiskz;

                // LOWE & MU are just copied from lomu_lowe & lomu_mu
                treeOut->Fill();	
                nqiskz_save = 0, nhitaz_save = 0;
	}
	// Close root
	// Save output file
	fileIn->Close();
	fileOut->Write();
	fileOut->Close();

	// delete
	delete lomu_head;
	delete TQII;
	delete TQAA;
	delete lomu_lowe;
	delete lomu_mc;

	return 0;
}
