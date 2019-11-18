#include "skmc_manager.h"

SKMCManager* SKMCManager::theManager = 0;

SKMCManager::SKMCManager(){};
SKMCManager::~SKMCManager(){};

void SKMCManager::Initialize()
{

  // Create New TTree
  skmcTree = new TTree("skmc", "Information from detector simulation");

  this->Clear();

  // Create New Branches
  skmcTree->Branch("evID", &evID, "evID/I");
  skmcTree->Branch("numCap", &numCap, "numCap/I");
  skmcTree->Branch("particleCap", &particleCap);
  skmcTree->Branch("materialCap", &materialCap);
  skmcTree->Branch("neutronEnergyCap", &neutronEnergyCap);
  skmcTree->Branch("XpositionCap", &XpositionCap);
  skmcTree->Branch("YpositionCap", &YpositionCap);
  skmcTree->Branch("ZpositionCap", &ZpositionCap);
  skmcTree->Branch("timeCap", &timeCap);
  skmcTree->Branch("numGamCap", &numGamCap);

  skmcTree->Branch("iparentGam", &iparentGam);
  skmcTree->Branch("gammaEnergyGam", &gammaEnergyGam);
  skmcTree->Branch("XdirectionGam", &XdirectionGam);
  skmcTree->Branch("YdirectionGam", &YdirectionGam);
  skmcTree->Branch("ZdirectionGam", &ZdirectionGam);

}

void SKMCManager::SetTree( int* nevsk,
			   int* num_cap,
			   int ip_cap[10000],
			   int med_cap[10000],
			   float en_cap[10000],
			   float pos_cap[10000][3],
			   float t_cap[10000],
			   int num_gam[10000],
			   int ip_gam[1000][10000],
			   float en_gam[1000][10000],
			   float dir_gam[1000][10000][3]
			    )
{
  evID = *nevsk;
  numCap = *num_cap;

  for(int i=0; i<numCap;i++){
    particleCap.push_back(ip_cap[i]);
    materialCap.push_back(med_cap[i]);
    neutronEnergyCap.push_back(en_cap[i]);
    XpositionCap.push_back(pos_cap[i][0]);
    YpositionCap.push_back(pos_cap[i][1]);
    ZpositionCap.push_back(pos_cap[i][2]);
    timeCap.push_back(t_cap[i]);
    numGamCap.push_back(num_gam[i]);
    /*
    cout << XpositionCap.at(i) << endl;
    cout << YpositionCap.at(i) << endl;
    cout << ZpositionCap.at(i) << endl;
    cout << numGamCap.at(i) << endl;
    */
    for(int j=0; j < numGamCap.at(i); j++){
      //cout << i << " " << j << " " << ip_gam[j][i] << " " << en_gam[j][i] << endl;
      //cout << dir_gam[j][i][0] << " " << dir_gam[j][i][1] << " " << dir_gam[j][i][2] << endl;
      iparentGam.push_back(ip_gam[j][i]-1);
      gammaEnergyGam.push_back(en_gam[j][i]);
      XdirectionGam.push_back(dir_gam[j][i][0]);
      YdirectionGam.push_back(dir_gam[j][i][1]);
      ZdirectionGam.push_back(dir_gam[j][i][2]);
    }

  }
}

void SKMCManager::FillTree()
{
  skmcTree->Fill();
  this->Clear();
}

void SKMCManager::WriteTree()
{
  skmcTree->Write();
}

void SKMCManager::Clear()
{
  evID=0;
  numCap = 0;
  vector<int>().swap(particleCap);
  vector<int>().swap(materialCap);
  vector<float>().swap(neutronEnergyCap);
  vector<float>().swap(XpositionCap);
  vector<float>().swap(YpositionCap);
  vector<float>().swap(ZpositionCap);
  vector<float>().swap(timeCap);
  vector<int>().swap(numGamCap);

  vector<int>().swap(iparentGam);
  vector<float>().swap(gammaEnergyGam);
  vector<float>().swap(XdirectionGam);
  vector<float>().swap(YdirectionGam);
  vector<float>().swap(ZdirectionGam);
}
