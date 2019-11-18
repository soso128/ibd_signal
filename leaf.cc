//
// leaf.cc : Low Energy Additional Fitting
// made by Hiroshi Ito
//        in 22 July 2019
//
// (intput) root file
//    HEAD, TQREAL, TQAREAL, LOWE, MC
// (output) root file
//    HEAD, TQREAL, TQAREAL, LOWE, MC, Cherenkov angle, Q50, etc...
//

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TVector3.h"

#include "TVirtualFitter.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "TreeManager.h"
#include "skroot.h"
#include "tqrealroot.h"
#include "mcinfo.h"
#include "loweroot.h"

#include "skbadcC.h"
#include "sktqC.h"
#include <iostream>
//#include <skbadcC.h>

#include "skparmC.h"
#include "geotnkC.h"
#include "geopmtC.h"

#include "skheadC.h"
//#include "skparmC.h"
//#include "geopmtC.h"
//#include "skbadcC.h"

#include <skwaterlenC.h>
#include <DataDefinition.h>
#include <Reformatter.h>
#include <vcworkC.h>
#include <zbsC.h>
#include <zbsmc.h>
#include "apmueC.h"

#include "typepdf.hh"
#include "./geo.hh"

SuperManager* Smgr;
TreeManager* mgr;
Header*          HEAD;
TQReal*          TQREALINFO;
TQReal*          TQAREALINFO;
LoweInfo*        LOWE;
AtmpdInfo*       ATMPD;
UpmuInfo*        UPMU;
MuInfo*          MU;
SLEInfo*         SLE;
MCInfo*          MC;

int NMIS;
int MISCH[1000];
int like_all;
int like_fail;
bool verbose;
TH1D *hcheren;

double cherenkov(int isMC, bool verbose,int &ierr);
double GetLinearLikelihood(TH1D *hcheren);
int getNhits(Float_t *v, int start_index, Float_t width, int nhits);
Double_t lfrsqrt(Double_t a2, Double_t b2, Double_t c2);
void MinimumWidthWindow(Float_t *t, int ALLHITS, float width, int &maxnhits, int &maxindex);
void MaxHitsInWindow(Float_t *t, int ALLHITS, float width, int &maxnhits, int &maxindex);
Bool_t checkMiss (const Int_t cab);
float dwallf( float *pos ); 
float effwallf(int idod, float *pos, float *dir, float &poswx, float &poswy, float &poswz);
extern "C" void muechk_( float* , int* );
void lfcharge(float twindow, float *vertex, float &charge);

#define nstep 100
// MAXPM has alreadly defined at tqrealroot.h
const double angledown=0.,angleup=90.;

extern "C" {
  void  init_geom_(void);
  void  skbadopt_(int *);
  void  skbadch_(int *, int *, int *);
}


int main(int argc, char* argv[])
{
    // check arguments
    if (argc < 3) {
	cout << endl;
	cout << "    usage: " << argv[0]
	     << " input file " 
	     << " output file " << endl;
	cout << endl;
	return 1;
    }

    //
    // open skroot file (read only)
    //
    int id = 10;
    string fname_in = argv[1];
    string fname_out = argv[2];

    // obtain tree
    Smgr = SuperManager::GetManager();
    Smgr ->CreateTreeManager(id,"",fname_out.c_str(),0);
    mgr = Smgr->GetTreeManager(id);
    mgr -> SetInputFile(fname_in.c_str());
    mgr -> Initialize();

    int ierr;
    int nmue; 
    float q50;
    float angle, pilike;
    float ovaq; 
    float poswall[3]={0};
    float dwall, effwall;
    float mc_npv[3]={0};
    char name[1000];

    TTree*tree= mgr->GetTree();
    TTree* theOTree = mgr->GetOTree();

    TBranch* headerbranch = 0;
    TBranch* tqrealbranch = 0;
    TBranch* tqarealbranch = 0;
    TBranch* lowebranch = 0;
    TBranch* mubranch = 0;
    TBranch* slebranch = 0;
    TBranch* mcbranch = 0;

    HEAD         = new Header;
    LOWE         = new LoweInfo;
    MU           = new MuInfo;
    SLE          = new SLEInfo;
    MC           = new MCInfo;
    TQREALINFO   = new TQReal;
    TQAREALINFO  = new TQReal;

    tree->SetBranchAddress("HEADER", &HEAD , &headerbranch );
    tree->SetBranchAddress("TQREAL", &TQREALINFO , &tqrealbranch );
    tree->SetBranchAddress("TQAREAL", &TQAREALINFO , &tqarealbranch );
    tree->SetBranchAddress("LOWE" , &LOWE , &lowebranch );
    tree->SetBranchAddress("MU" , &MU , &mubranch );
    tree->SetBranchAddress("SLE" , &SLE  , &slebranch  );
    tree->SetBranchAddress("MC" , &MC  , &mcbranch  );

    theOTree->SetBranchAddress("HEADER", &HEAD , &headerbranch );
    theOTree->SetBranchAddress("TQREAL", &TQREALINFO , &tqrealbranch );
    theOTree->SetBranchAddress("TQAREAL", &TQAREALINFO , &tqarealbranch );
    theOTree->SetBranchAddress("LOWE" , &LOWE , &lowebranch );
    theOTree->SetBranchAddress("MU" , &MU , &mubranch );
    theOTree->SetBranchAddress("MC" , &MC  , &mcbranch  );
    theOTree->Branch("angle",&angle,"angle/F");
    theOTree->Branch("pilike",&pilike,"pilike/F");
    theOTree->Branch("q50",&q50,"q50/F");
    theOTree->Branch("ovaq",&ovaq,"ovaq/F");
    theOTree->Branch("nmue", &nmue, "nmue/I");
    theOTree->Branch("dwall",&dwall,"dwall/F");
    theOTree->Branch("effwall",&effwall,"effwall/F");
    theOTree->Branch("poswall[3]",&poswall,"poswall[3]/F");
    theOTree->Branch("mc_npv[4]",&mc_npv,"mc_npv[4]/F");// (x,y,z,t)/(cm, ns)

    like_all = 0;
    like_fail = 0;
    hcheren = new TH1D("hcheren", "Opening Cherenkov Angle from 3-PMT combinations", nstep, angledown, angleup);

    //
    // bad ch file import (Takeuchi-code import by H.Ito, August, 2019)
    //
    tree->GetEntry(0);
    int nrun      = HEAD->nrunsk;//atoi(argv[3]);
    int nsub      = 1;
    int imask     = 23;
    int log_level = 1;
    int istat;
    skheadg_.sk_geometry = 4; // SK-IV
    combad_.log_level_skbadch = log_level; // outputs from skbadch

    //--------------------------------------
    // call Fortran routines
    //--------------------------------------
    skbadopt_(&imask);
    skbadch_(&nrun, &nsub, &istat);
    //--------------------------------------

    if(false){
      cout << "nrun/nsub/imask = "
	   << nrun  << " "
	   << nsub  << " "
	   << imask << " "
	   << endl;
      cout << "nbad/imaskbad/imaskbadopt/log_level_skbadch = "
	   << combad_.nbad  << " "
	   << combad_.imaskbad << " "
	   << combad_.imaskbadopt << " "
	   << combad_.log_level_skbadch << " "
	   << endl;
      return 0;
    }

    const int NMIS=combad_.nbad;
    int MISCH[NMIS];
    cout<<"nbadch "<<NMIS<<": ";
    for(int j=0;j<NMIS;j++){
      MISCH[j]=combad_.isqbad[j];
      cout<<MISCH[j]<<" ";
    }
    cout<<endl;
      
    Int_t pn,max1,max2,min1,min2;

    // main loop
    int ntotal = tree->GetEntries();
    for (int i = 0; i < ntotal; i++) {
      tree->GetEntry(i);

      // *** Cherenkov angle calc
      //   option 1) input format: root 0, zbs 1.
      //          2) verbose: [false] quite, [true] comment
      angle = cherenkov(0,false,ierr);

      // *** Iida's pilike
      pn = hcheren->GetMaximumBin();
      min1 = pn - 3;      max1 = pn + 3;
      min2 = pn - 10;     max2 = pn + 10;
      pilike = ( hcheren->Integral( min1 + 1, max1 - 1 ) )
              /( hcheren->Integral( min2 + 1, max2 - 1 ) ); 
      //pilike = hcheren->Integral(min1+1,max1-1) /
      //  (hcheren->Integral(min2+1,max2-1) - hcheren->Integral(min1+1,max1-1));
      //cout << "  pilike: " << pilike << endl;

      // *** Q50 calc
      lfcharge(50.,LOWE->bsvertex,q50);

      // *** bonsai ovaQ
      ovaq = pow(LOWE->bsgood[1],2) - pow(LOWE->bsdirks,2);
	
      // *** dwall & effwall
      poswall[0]=0; poswall[1]=0; poswall[2]=0;
      dwall = dwallf( LOWE->bsvertex ); 
      effwall = effwallf(1, LOWE->bsvertex, LOWE->bsdir,poswall[0],poswall[1],poswall[2]);
      
      // *** nmue 
      int muerr = 1;
      muechk_( LOWE->bsvertex, &muerr );
      nmue = apmue_.apnmue;
      //cout << "  nmue: " << nmue << endl;


      // *** mc info
      for(int j=0;j<3;j++){
	MC->energy[j] = sqrt(pow(MC->pvc[j][0],2) 
			     + pow(MC->pvc[j][1],2)
			     + pow(MC->pvc[j][2],2) );
      }

      // *** mc neutron capture (vertex, time)
      mc_npv[0]= *reinterpret_cast<float*>(&MC->mcinfo[0]);// pos_x
      mc_npv[1]= *reinterpret_cast<float*>(&MC->mcinfo[1]);// pos_y
      mc_npv[2]= *reinterpret_cast<float*>(&MC->mcinfo[2]);// pos_z
      mc_npv[3]= *reinterpret_cast<float*>(&MC->mcinfo[3]);// time  

      // root file filling
      theOTree->Fill();
      mgr->Clear();

      if(ntotal>100 &&
	 i%(ntotal/100)==0 &&
	 i>0) cout<<"."<<flush;
    }
    cout<<endl;

    Smgr->DeleteTreeManager(id);
    SuperManager::DestroyManager();
    return 0;
}




//****************************************
// Cherenkov Angle calculation
// (option1) input format 0:root, 1:zbs
// (option2) debug true: on fault: off 
// (option3) error info 0:normal, -1:fail
//**************************************** 
double cherenkov(int isMC, bool verbose, int &ierr){

  ierr=0;
  bool limithits=true;  // use angle
  //limithits=false;    // use angle_15ns
  Int_t goodcount=0, ipmt;
  int icabbit=0xffff;
  if (isMC==0){  // root 
    for (int i=0;i<TQREALINFO->nhits;i++){
      int realcable = TQREALINFO->cables[i] & icabbit;
      int icflag = TQREALINFO->cables[i]>>16;
      if (realcable<= MAXPM && !checkMiss(realcable)){
	if (icflag &1){
	  goodcount+=1;
	}
      }
    }
  }
  /*  else if(isMC==1){ // zbs
    for (int i=0;i<sktqz_.nqiskz;i++){
      if (sktqz_.ihtiflz[i]&0x02 &&
	  sktqz_.icabiz[i] <= MAXPM &&
	  !combad_.ibad[sktqz_.icabiz[i]-1] &&
	  !checkMiss(sktqz_.icabiz[i]) &&
	  sktqz_.ihtiflz[i]&0x01){
	goodcount+=1;
      }
    }
    }*/
  

  const Int_t MAXN10 = 200;
  const int ALLHITS = goodcount;
  const Float_t C_WATER = 21.5833;

  Int_t   cabiz[ALLHITS], cabiz2[ALLHITS];
  Float_t tiskz[ALLHITS], tiskz2[ALLHITS];
  Float_t qiskz[ALLHITS], qiskz2[ALLHITS];
  Int_t   index[ALLHITS], nindex[MAXN10];
  Double_t dhit[3], dtmp;
  Float_t hitv_x[ALLHITS];
  Float_t hitv_y[ALLHITS];
  Float_t hitv_z[ALLHITS];
  Double_t rsqrd, opang, rnosq;
  Int_t icount=0;

  float pvx[3];
  for(int j=0;j<3;j++)
    pvx[j]=LOWE->bsvertex[j];

  /*
  if (isMC==1){
    for (int i=0; i<sktqz_.nqiskz; i++)
      {
	if (sktqz_.ihtiflz[i]&0x02 && sktqz_.icabiz[i] <= MAXPM &&
	    !combad_.ibad[sktqz_.icabiz[i]-1] &&
	    !checkMiss(sktqz.icabiz[i]) && 
	    sktqz_.ihtiflz[i]&0x01) 
	  {
	    cabiz2[icount] = sktqz_.icabiz[i];
	    tiskz2[icount] = sktqz_.tiskz[i];
	    qiskz2[icount] = sktqz_.qiskz[i];
	    icount+=1;
	  }
      }
  }
  else*/
  if(isMC==0){
    if (verbose)
    cout << "About to load up TQREAL with " << TQREALINFO->nhits << endl;
    for (int i=0; i<TQREALINFO->nhits; i++){
      int realcable = TQREALINFO->cables[i] & icabbit;
      int icflag = TQREALINFO->cables[i]>>16;
      if (realcable <= MAXPM && !checkMiss(realcable)){
	if (icflag&1){
	  if (verbose)
	  cout<< "realcable "<<realcable<<" i "<<i <<std::endl;
	  cabiz2[icount] = realcable;
	  tiskz2[icount] = TQREALINFO->T[i];
	  qiskz2[icount] = TQREALINFO->Q[i];
	  icount+=1;
	}
      }
    }
  }

  // TOF subtraction for all hits in copied arrays
  for (int i=0; i<ALLHITS; i++) {
    Float_t tof;
    tof = TMath::Sqrt((pvx[0] - xyz[cabiz2[i]-1][0]) * (pvx[0] - xyz[cabiz2[i]-1][0])
		      +(pvx[1] - xyz[cabiz2[i]-1][1]) * (pvx[1] - xyz[cabiz2[i]-1][1])
		      +(pvx[2] - xyz[cabiz2[i]-1][2]) * (pvx[2] - xyz[cabiz2[i]-1][2]))
                 / C_WATER;
    tiskz2[i] -= tof;
  }
  // Sort hits by TOF-corrected time
  TMath::Sort(ALLHITS, tiskz2, index, kFALSE); // In increasing order
  for (int i=0; i<ALLHITS; i++){
    cabiz[i] = cabiz2[ index[i] ];
    tiskz[i] = tiskz2[ index[i] ];
    qiskz[i] = qiskz2[ index[i] ];
  }

  int N15      = 0;
  int N15index = 0;
  
  if (limithits)   {
    N15 = 183;
    MinimumWidthWindow(tiskz, ALLHITS, 15., N15, N15index);
  }
  else {
    MaxHitsInWindow(tiskz, ALLHITS, 15., N15, N15index);
  }

  int maxncmb = N15*N15*N15 / 6;
  vector<double> abc2(N15*N15, 0.0);
  vector<double> maxcmb(maxncmb, 0.0);

  //making unit direction vector form vertex to selected tubes
  if (verbose)
    printf("N15 is %d at index %d.\n", N15, N15index);

  //zero histogram          
  hcheren->Reset();

  for (int ii=N15index; ii < N15index+N15; ii++){
    for (int jj=0; jj<3; jj++){
      dhit[jj] = xyz[cabiz[ii]-1][jj] - pvx[jj];
    }
    dtmp= 1.0/TMath::Sqrt(dhit[0]*dhit[0] +
			  dhit[1]*dhit[1] + 
			  dhit[2]*dhit[2]);
    int intmp= ii - N15index;
    hitv_x[intmp]= dhit[0]*dtmp;
    hitv_y[intmp]= dhit[1]*dtmp;
    hitv_z[intmp]= dhit[2]*dtmp;
  }
  int tmpindex_a=0, tmpindex_b=0, tmpindex_c=0;
  //Calculate lengths of difference vectors between directions
  for (int aa = 0; aa < N15-1; aa++) {
    for (int bb=aa+1; bb<N15; bb++) {
      tmpindex_a=aa*N15+bb;
      abc2[tmpindex_a] = 
	 (hitv_x[aa]-hitv_x[bb])*(hitv_x[aa]-hitv_x[bb]) 
	+(hitv_y[aa]-hitv_y[bb])*(hitv_y[aa]-hitv_y[bb])
	+(hitv_z[aa]-hitv_z[bb])*(hitv_z[aa]-hitv_z[bb]);
    }
  }
  //Do all direction triangles and fill histogram                                                                        
  for (int aa=0;aa<N15-2;aa++){
    for (int bb=aa+1;bb<N15-1;bb++){
      tmpindex_a=aa*N15+bb;
      dtmp= abc2[tmpindex_a];
      for (int cc=bb+1;cc<N15;cc++){
	tmpindex_b=aa*N15+cc;
	tmpindex_c=bb*N15+cc;
	rsqrd = lfrsqrt(dtmp, abc2[tmpindex_b], abc2[tmpindex_c]);
	rnosq = TMath::Sqrt(rsqrd);
	opang = TMath::ASin(rnosq)*180.0/TMath::Pi();
	hcheren->Fill(opang);
      }
    }
  }
  char xtitle[64];

  //Finding peak of max summation of neighboring 7 bins
  TString che_angle;
  Double_t unpack[nstep], b2deg=0.0, lsum=0.0;
  int nwindow=7, midpoint=3, look=100;
  int ilowbnd=41, ihibnd=67;
  int nrange, llook, ij=0, locale=0;
  Double_t lheight=0, peakangle=0.0;

  float plotcontent[nstep],angle,angle_likelihood, 
    angle_15ns,angle_likelihood_15ns;

  nrange = nstep-nwindow+1;
  b2deg  = 90.0/100.0;

  int anglebinnum= nstep;
  for(int hh=0;hh<nstep;hh++){
    unpack[hh]=hcheren->GetBinContent(hh+1);
    plotcontent[hh] = hcheren->GetBinContent(hh);
  }
  for(int ii=0;ii<nrange;ii++){
    lsum=0.0;
    for (int jj=0;jj<nwindow;jj++){
      ij=ii+jj;
      lsum+=unpack[ij];
    }
    if (lheight<lsum){
      lheight=lsum;
      locale=ii+1+midpoint;}
  }
  peakangle=locale*b2deg;
  if (verbose)
    printf("Cherenkov angle is %f\n",  peakangle);

  double linear_like = GetLinearLikelihood(hcheren);

  if (verbose)
    printf("likelihood value is %f.\n", linear_like);

  if (limithits) {
    angle =  peakangle;
    angle_likelihood = linear_like;
  }
  else {
    angle_15ns = peakangle;
    angle_likelihood_15ns = linear_like;
  }
  return angle;
}


double GetLinearLikelihood(TH1D *hcheren){
  
  float pli0=0.0, pli1=0.0;
  TF1 *pdfFcn = new TF1("typpdf", typPdf, 10, 90, 2);
  // TF1 option: name, formula
  //             xmin, xmax, npar
  // typPdf = p0 + p1 x
  pdfFcn->SetParameters(-0.00277, 0.000277);
  
  double nentry = hcheren->GetEntries();
  TH1D *hnorche= (TH1D*)hcheren->Clone();
    
  hnorche->Scale(1.0/nentry);
  TVirtualFitter::SetDefaultFitter("Minuit2");
  TFitResultPtr r = hnorche->Fit(pdfFcn, "RWLSQ");
  // fit option ... 
  //    = "R"  Use the Range specified in the function range
  //    = "WL" Use Loglikelihood method and bin contents are not integer,
  //    = "S"  The result of the fit is returned in the TFitResultPtr
  //      (see below Access to the Fit Result)
  //     = "Q"  Quiet mode (minimum printing)

  like_all++;
  if (!r->IsValid())
    like_fail++;

  pli0 = pdfFcn->GetParameter("p0");
  pli1 = pdfFcn->GetParameter("p1");
  
  Double_t  logm, yi, caltmp, sumlikeli=0.0, binlike=0.0;
  int nvpar, nparx;
  int numbin = hnorche->GetNbinsX();
  int fitlow = hnorche->FindBin(10.0);
  Double_t hwidth = hnorche->GetBinWidth(1);
  
  //calculate likelihood
  for (int jj=fitlow;jj<numbin;jj++){
    yi=hnorche->GetBinContent(jj);
    caltmp=(pli1*jj*hwidth)+pli0;
    if (caltmp>0){
      binlike=yi*TMath::Log10(caltmp);}
    else {binlike=0.0;}
    sumlikeli+=binlike;
  }
  logm= sumlikeli;
  delete pdfFcn;
  return logm;
}


int getNhits(Float_t *v, int start_index, Float_t width, int nhits){
  int i = start_index;
  while (1) {
    i++;
    if((i > nhits-1 ) || (TMath::Abs((v[i]-v[start_index])) > width)) break;
  }
  return TMath::Abs(i - start_index);
}


Double_t lfrsqrt(Double_t a2, Double_t b2, Double_t c2)
{
  Double_t rsquare=0.0;
  Double_t tmpr2=0.0;
  tmpr2= a2*(2.*b2-a2) + b2*(2.*c2-b2) + c2*(2.*a2-c2);
  if (tmpr2>0){
    rsquare= a2*b2*c2/tmpr2;
  }else {
    return 2;
  }
  return rsquare;
}


void MaxHitsInWindow(Float_t *t, int ALLHITS, float width, int &maxnhits, int &maxindex)
{
  maxnhits = 0;
  maxindex = 0;

  //Finding max N15 and its index, limiting potential number of hits
  for (int in = 0; in < ALLHITS; in++) {
    // Calculate hits in 15 ns window                          
    int N = getNhits(t, in, width, ALLHITS);
    if ( N > maxnhits) {
      maxnhits = N;
      maxindex = in;
    }
  }                         
  if (verbose)
    cout << "MaxHitsInWindow, maxnhits=" << maxnhits << ", maxindex=" << maxindex << endl;
}


void MinimumWidthWindow(Float_t *t, int ALLHITS, float width, int &maxnhits, int &maxindex)
{
  double dt   = 100000.;
  for (int i = 0; i < ALLHITS-maxnhits; i++) {
    if (t[i+maxnhits] - t[i] < dt) {
      dt = t[i+maxnhits] - t[i];
      maxindex = i;
    }
  }
  if (dt > width) {
    if (verbose)
      cout << "dt too wide: ";
    MaxHitsInWindow(t, ALLHITS, width, maxnhits, maxindex);
  }
  else {
    if (verbose)
      cout << "MinimumWidthWindow, maxnhits=" << maxnhits << ", maxindex=" << maxindex << endl;
  }
}



//
// check miss and bad chennel PMT
//
Bool_t checkMiss (const Int_t cab)
{
  for (Int_t i=0; i<NMIS; i++) {
    if ( MISCH[i] == cab ) return kTRUE;
  }
  return kFALSE;
}


//
// dwall 
//
float dwallf( float *pos )
{
  double rr, r;
  double dwall, z1, z2;

  rr = pos[ 0 ] * pos[ 0 ] + pos[ 1 ] * pos[ 1 ];
  r  = sqrt( rr );

  dwall = RINTK - r;
  z1 = ZPINTK - pos[ 2 ];
  z2 = pos[ 2 ] - ZMINTK;

  if ( z1 < dwall ) { 
    dwall = z1; 
  }
  if ( z2 < dwall ) { 
    dwall = z2; 
  }

  return dwall;
} 


//
// effective wall length
//
float effwallf(int idod, float *pos, float *dir,
	       float &poswx, float &poswy, float &poswz)
{
  float  a,b,c,t,x,y,z;

  if ( idod == 1 ){  //from ID wall
    // calc cross point at the barrel wall
    a = dir[0]*dir[0] + dir[1]*dir[1];
    b = dir[0]*pos[0] + dir[1]*pos[1];
    c = pos[0]*pos[0] + pos[1]*pos[1] - (RINTK)*(RINTK);
    t = (-b - sqrt(b*b - a*c)) / a;
    x = dir[0]*t + pos[0];
    y = dir[1]*t + pos[1];
    z = dir[2]*t + pos[2];

    if ( z > ZPINTK) {   //if crossed at top
      t = (ZPINTK - pos[2]) / dir[2];
      x = dir[0]*t + pos[0];
      y = dir[1]*t + pos[1];
      z = ZPINTK;
    }
    if ( z < ZMINTK) {   //if crossed at bottom
      t = (ZMINTK - pos[2]) / dir[2];
      x = dir[0]*t + pos[0];
      y = dir[1]*t + pos[1];
      z = ZMINTK;
    }

  }  else if ( idod == 0 ) {  // from OD wall
    //calc. cross point at barrel wall
    a = dir[0]*dir[0] + dir[1]*dir[1];
    b = dir[0]*pos[0] + dir[1]*pos[1];
    c = pos[0]*pos[0] + pos[1]*pos[1] - (RTKTK)*(RTKTK);
    t = (-b - sqrt(b*b - a*c)) / a;
    x = dir[0]*t + pos[0];
    y = dir[1]*t + pos[1];
    z = dir[2]*t + pos[2];

    if (z > ZPTKTK) { // if crossed at top
      t = (ZPTKTK - pos[2]) / dir[2];
      x = dir[0]*t + pos[0];
      y = dir[1]*t + pos[1];
      z = ZPTKTK;
    }

    if (z < ZMTKTK) { // if crossed at bottom
      t = (ZMTKTK - pos[2]) / dir[2];
      x = dir[0]*t + pos[0];
      y = dir[1]*t + pos[1];
      z = ZMTKTK;
    }
  }

  float effwallfs = sqrt((pos[0]-x)*(pos[0]-x)+(pos[1]-y)*(pos[1]-y)+(pos[2]-z)*(pos[2]-z));
  poswx=x;
  poswy=y;
  poswz=z;

  return effwallfs;
}



//
// charge integral determine
//
void lfcharge(float twindow, float *vertex, float &charge)
{
  int j=0;
  const Float_t C_WATER = 21.5833;
  Float_t tof = 0;
  int icflag = 0;
  int icount = 0;
  int icabbit=0xffff;
  const int nhitmax = 4000;
  int cabiz[nhitmax]={0};
  float tiskz[nhitmax]={0};
  float qiskz[nhitmax]={0};
  
// ---- get time ordered tof subtracted (vertex)
  for (int i=0; i<TQREALINFO->nhits; i++){
    j=TQREALINFO->cables[i] & icabbit;
    icflag = TQREALINFO->cables[i]>>16;
    if (icflag&1 && !checkMiss(j)){
      cabiz[icount] = j;//TQREALINFO->cables[i];
      tiskz[icount] = TQREALINFO->T[i];
      qiskz[icount] = TQREALINFO->Q[i];
      tof = TMath::Sqrt(  pow(vertex[0] - xyz[cabiz[icount]-1][0], 2)
		        + pow(vertex[1] - xyz[cabiz[icount]-1][1], 2)
		        + pow(vertex[2] - xyz[cabiz[icount]-1][2], 2)) / C_WATER;
      tiskz[icount] -= tof;
      icount++;
    }
  }
  // sort has already done!
  if(false){
    float dmy[4]={0};
    for(int k=0;k<icount-1;k++){
      for(int m=k;m<icount-1;m++){
	if(tiskz[k] > tiskz[m+1]){
	  dmy[0] = tiskz[k]; tiskz[k] = tiskz[m+1]; tiskz[m+1] = dmy[0];
	  dmy[1] = cabiz[k]; cabiz[k] = cabiz[m+1]; cabiz[m+1] = dmy[0];
	  dmy[2] = qiskz[k]; qiskz[k] = qiskz[m+1]; qiskz[m+1] = dmy[0];
	}}}  
  }


// ---- find the center of the distribution:
  int nhit=0, nhit_max=0, j0=0;
  if (icount > 0){
    for(int i = 0; i < icount; i++){
      nhit  = 0;
      j = i+1;
      while (j < icount  &&  tiskz[j]-tiskz[i] <= twindow){
	nhit++;
	j++;
      }
      if(nhit_max < nhit){
	nhit_max = nhit;
	j0=i;
      }
    }
  }
  else{
    nhit = 0;
  }

  charge=0.;
  j=j0+1;
  while (j < icount  &&  tiskz[j]-tiskz[j0] <= twindow){
    charge+=qiskz[j];
    j++;
  }
}


