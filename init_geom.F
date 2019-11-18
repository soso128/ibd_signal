      subroutine init_geom()
      implicit none

#include "skhead.h"
#include "skparm.h"
#include "geopmt.h"
      
c***  set SK_GEOMETRY      
      SK_GEOMETRY = 4
      
c***  init
      call geoset

      return 
      end
