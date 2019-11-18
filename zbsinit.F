      subroutine zbsinit
      implicit none

#include "skhead.h"
#include "skparm.h"
#include "geopmt.h"
#include "sktq.h"

      integer idim
      parameter (idim=4000000)
      real zbs(idim)
      common/kzbs/zbs

C      CALL KZINIT
      SK_FILE_FORMAT = 0
      call kzlimit(idim)
      return
      end
