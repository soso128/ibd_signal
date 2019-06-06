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
#include "stdio.h" // check whether file exist
#include <iomanip> // set io
#include "skroot.h"

const Float_t mintime = 18000; //ns
const Float_t maxtime = 535000;
const Float_t deadtime = 900;
const Int_t NEVENTGAP = 1500; //~2000 events per MC file, 2 MC/1 dummy, so 1500 is enough
Int_t twind = 517000;
Int_t NHITUP = 55000;
Int_t NHITLOW = 50000;
Float_t CENTERT = 500000;

int   access(const   char   *filename,   int   amode);

extern "C" {
	void initialize_(void);
}
//
int main(int argc, char **argv)
{
	if (argc != 5)
	{
		cout << "Usage: " << argv[0] << " fname_out fname_in run_number" << endl;
		exit(1);
	}

	TString fname_out = argv[1];
	TString fname_in = argv[2];
	TString n = argv[3];
	Int_t nfile = n.Atoi();
	cout << "Input file: " << fname_in << endl;
	cout << "Output file: " << fname_out << endl;
	cout << "File number: " << n << endl;
	Int_t filenumber = n.Atoi();

	//start combine
	cout << "Start of comb MC " << endl;
	TChain *t2kch = new TChain("data");
	//dummy trigger chain
        char inname[500];
        sprintf(inname, "/disk02/usr6/elhedri/t2k/data/t2k.0%d.root", nfile);
	t2kch->Add(inname);
	NHITUP=60000;
	NHITLOW=55000;
	cout << "Processing run " << n << endl;
//	t2kch->Add("/disk01/lowe5/zhangyang/2p2/t2k/*.root");
	//t2kch->Add("/disk/lowe5/zhangyang/2p2/t2k/t2k.069700.root");
	//t2kch->Add("/disk/lowe5/zhangyang/2p2/t2k_old_data/*.root");

	TQReal *TQI = new TQReal;
	TQReal *TQA = new TQReal;
	Header *HEADER = new Header;

	t2kch->SetBranchAddress("TQREAL", &TQI);
	t2kch->SetBranchAddress("TQAREAL", &TQA);
	t2kch->SetBranchAddress("HEADER", &HEADER);
	// Input file
	//
	//nread is for MC

	Float_t toffset;
	Int_t nread=0;
	Int_t entry=nfile*NEVENTGAP; // 2 MC event ~ 1 dummy trigger event
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
        int *is_signal = new int[100000];
	treeOut->Branch("issignal",  &is_signal, "issignal[100000]/I");
	treeOut->Branch("HEADER",  "Header", &lomu_head, bufsize, 0);
	//        treeOut->Branch("TQREAL",  "TQReal", &lomu_tqi,  bufsize, 0);
	//        treeOut->Branch("TQAREAL", "TQReal", &lomu_tqa,  bufsize, 0);
	treeOut->Branch("TQREAL",  "TQReal", &TQII,  bufsize, 0);
	treeOut->Branch("TQAREAL", "TQReal", &TQAA,  bufsize, 0);
	treeOut->Branch("LOWE", "LoweInfo", &lomu_lowe,  bufsize, 0);
	treeOut->Branch("MC", "MCInfo", &lomu_mc,  bufsize, 0);
	//open previously created event list of input sample events

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
        int nsep = 2;
	cout << "Get 2p2 candidate entries " <<  treeIn->GetEntries() << endl;
	for (Int_t I2P2=0; I2P2<nentries; I2P2 ++) {
            cout << "Event " << I2P2 << endl;
//	for (Int_t I2P2=0; I2P2<100; I2P2 ++) {
		//                if (I2P2%2001 == 0)  continue;  //exclude overlapping events
//		cout<<I2P2<<endl;
		treeIn->GetEntry(I2P2);
		//1 dummy trigger events with  2 MC event, hence nread%2 conditional
		//nread is for MC
		if(nread%nsep == 0) 
		{
			//read dummy trigger
			t2kch->GetEntry(entry);
			entry++;
			//a goog dummy events should have nhits in 50000~55000
			while (TQI->nhits > NHITUP || TQI->nhits < NHITLOW)
			{
				t2kch->GetEntry(entry);
				entry++;
			}
		}
//		cout<<I2P2<<" "<<entry<<" "<<nread<<endl;
		if (nread%1000 == 0) cout << "Read " << nread << " events. " << endl;
		//dummy trigger timing goes from -500us < t < 500us 
		//get the time offset at the beginnig of time window and 
		//the center of the time window 
		if (nread%nsep ==0)
		{
			//Linyan
//			toffset = TQI->T[1000]; //shift to 1000th hit
			toffset = TQI->T[0]; //shift to 1000th hit
		}
		else if (nread % nsep == 1)
		{
			//Linyan
//			toffset = TQI->T[TQI->nhits-1] - CENTERT; // set toffset to center.
			toffset = TQI->T[TQI->nhits-1] - twind ; // set toffset to center.
			//cout << "*************TOFFSET********************************** " << toffset << endl;
		}
		//Linyan
		//nread++;
		//Temporary TQ information to save output TQ 
		//Int_t nqiskz, nhitaz, nqisk_raw, nhitaz_raw;
		//nqiskz = sktqz_.nqiskz;
		//nhitaz = sktqaz_.nhitaz;

		for (Int_t i = 0; i < lomu_tqi->nhits; i++)
		{
			//ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,16); //left shift 16 bit
			//cable hit flag, default 
			ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,0); //left shift 16 bit
			//lower 16 bit for calbe, and higher 16 bit for flag
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
                        if (TQI->T[i] >= toffset && TQI->T[i] <= toffset + twind)
                        {
                                tiskz[nqiskz_save] = TQI->T[i]-toffset +mintime;
                                qiskz[nqiskz_save] = TQI->Q[i];
                                icabiz[nqiskz_save] = TQI->cables[i]&0x0000FFFF;
                                //had problems if I just directly copy ihtiflz, so just set every 
                                //hit to the same default value.
                                //ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,16);
                                ihtiflz[nqiskz_save] = 2*(Int_t)TMath::Power(2,0);
                                nqiskz_save++;
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
                //TMath::Sort(nhitaz_save, taskz2, indexa, kFALSE);
                for(int i=0;i<nqiskz_save;i++){
                        tiskz[i] = tiskz2[index[i]];
                        qiskz[i] = qiskz2[index[i]];
                        icabiz[i] = icabiz2[index[i]];
                        if (index[i] < lomu_tqi->nhits)
                            is_signal[i] = 1;
                        else
                            is_signal[i] = 0;
                        //cout << "tiskz " << tiskz[i] << endl;
                }
                //cout << "***********************normal*********************************" << endl;
                //If a hit is okay, save it in sktqz_ timing bank
                // Set TQI 
                TQII->Clear();
                Int_t newnqiskz = 0;
                //cout << "***********************newnqiskz*************** " << newnqiskz  << endl;
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
		nread++;
                //if (nread==10000) break;
	}
	// Close root
	// Save output file
	fileIn->Close();
	fileOut->Write();
	fileOut->Close();
	cout << "Combining finished. Read " << nread << " events. " << endl;

	// delete
	delete lomu_head;
	delete TQII;
	delete TQAA;
	delete lomu_lowe;
	delete lomu_mc;

	return 0;
}
