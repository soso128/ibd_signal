      subroutine initialize()
      implicit none

#include "skhead.h"
#include "skparm.h"
#include "geopmt.h"
#include "sktq.h"
      
c***  set SK_GEOMETRY      
      SK_GEOMETRY = 4
      
c***  init
      SK_FILE_FORMAT = 1
      call kzinit
      call skoptn('31,30,26,25')
      call geoset

      return 
      end
