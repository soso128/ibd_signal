
#include <TTree.h>
#include <TFile.h>

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>

using namespace std;

class SKMCManager {

 protected:
  static SKMCManager* theManager;

  virtual ~SKMCManager();
  SKMCManager();

  TTree *skmcTree;

  // Definision of Variables
  Int_t evID;
  Int_t numCap;

  vector<int> particleCap;
  vector<int> materialCap;
  vector<float> neutronEnergyCap;
  vector<float> XpositionCap;
  vector<float> YpositionCap;
  vector<float> ZpositionCap;
  vector<float> timeCap;
  vector<int> numGamCap;

  vector<int> iparentGam;
  vector<float> gammaEnergyGam;
  vector<float> XdirectionGam;
  vector<float> YdirectionGam;
  vector<float> ZdirectionGam;

 public:

  static SKMCManager* GetManager() {
    if ( SKMCManager::theManager == NULL ) {       
      SKMCManager::theManager = new SKMCManager;
      return SKMCManager::theManager;
    }
    else { 
      return SKMCManager::theManager;
    }
  }

  void Initialize();
  void SetTree( int* nevsk,
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
		);
  void FillTree();
  void WriteTree();
  void Clear();

};
