       subroutine set_rflistC(lun,fnum1, fnum2)
#include "skhead.h"
#include "skparm.h"
#include "sktq.h"
#include "geopmt.h"
#include "skroot_lowe.h"
#include "skday.h"
#include "skwt.h"
#include "skbadc.h"

#include "vcwork.h"
#include "vcvrtx.h"

           integer lun,fnum1, fnum2
           character*1000 filename
           character*8 fmt, fmt2
           character*4 nam1
           character*2 nam2
           fmt = '(I4.4)'
           fmt2 = '(I2.2)'
           write(nam1, fmt) fnum1
           write(nam2, fmt2) fnum2
           filename="/disk02/usr6/yashida/relic/mc/atmnu/sk4.apr16_lowfit_20191005/zbs/lowfit.atmnu."//nam1//"_"//nam2//".zbs"
           call set_rflist(lun, filename,'LOCAL',' ','RED',' ',' ',
     &     'recl=5670 status=old',' ',' ' )
           return
       end

