
#include "skmc_manager.h"

extern "C" void skmc_initialize_(){

  SKMCManager* mgr = SKMCManager::GetManager();
  mgr->Initialize();

}

extern "C" void skmc_set_tree_(int *nevsk,
			       int *num_cap,
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

  SKMCManager* mgr = SKMCManager::GetManager();
  mgr->SetTree( nevsk,
		num_cap,
		ip_cap,
		med_cap,
		en_cap,
		pos_cap,
		t_cap,
		num_gam,
		ip_gam,
		en_gam,
		dir_gam
		);

  mgr->FillTree();

}

extern "C" void skmc_write_()
{

  SKMCManager* mgr = SKMCManager::GetManager();
  mgr->WriteTree();

}
