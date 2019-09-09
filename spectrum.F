**********************************************************************
*     ----------------------------------------
      SUBROUTINE spectrum(iflag,idate,itime,pos,seed,emin, emax)
*     ----------------------------------------
*     
*     ( purpose )
*       VECTOR GENERATION ( SRN NEUTRINO )
*
*     ( input ) 
*       IFLAG : see spcsunsl.F
*       POS(3)   : VERTEX POSITION
*       idate    : run date
*       itime    : run time
*       pos      : vertex position
*       seed     : random seed for spectrum generation
*       emin, emax: energy bin in which to generate events
*
*     ( output ) 
*       COMMON VCWORK,VCVRTX
*
*     ( creation date and author )
*       1995.11.16 ; First version by Y.Koshio
*       1998.11.05 ; modified by Y.Takeuchi
*                         without HEP (only B8) version
*
**********************************************************************
#include "vcwork.h"
#include "vcvrtx.h"

      real pos(3,10),postmp(3)
      real*8 ep, enu, eneutron, mneutron, mproton
      integer*4 seed
      integer idate(3),itime(3),neutron,ip(3)
      real dirtmp(3), dir(3,3), pmass,pmomabs(3),pmom(3,10)
      real dwall,tvtx(3)
      character*200 specparamfile

      call kzeclr

      nvct=3
      dwall=200.0
      mproton=938.272
      mneutron=939.5654

      ip(1)= -12
      ip(2)= -11
      ip(3)= 2112
      call rnpos(postmp,dwall)
      do j = 1,3
          pos(1,j)=postmp(1)
          pos(2,j)=postmp(2)
          pos(3,j)=postmp(3)
          tvtx(j)=0.0
      enddo
      call rndir(dirtmp)
      dir(1,1)=dirtmp(1)
      dir(2,1)=dirtmp(2)
      dir(3,1)=dirtmp(3)
c     Use Strumia-Vissani cross-section to get positron direction
c     and deduce the one for the neutron
      ep = emin 
      pmomabs(2)=sqrt(ep*ep-0.511*0.511)
      eneutron = enu + mproton - ep
      pmomabs(3)=sqrt(eneutron * eneutron - mneutron * mneutron)
      do j = 1,3
          pmom(j,1) = dir(j,1) * pmomabs(1);
          pmom(j,2) = dir(j,2) * pmomabs(2)
          pmom(j,3) = pmom(j,1) - pmom(j,2)
      enddo

      write(*,*) '**************************'
      write(*,*) 'EV, Vertex::',i,pos(1,1),pos(2,1),pos(3,1)
      write(*,*) '1st positron::',pmomabs(1),dir(1,1),dir(2,1),dir(3,1)
      write(*,*) '1st neutron::',pmomabs(2),dir(1,2),dir(2,2),dir(3,2)
      write(30,*) pos(1,1), pos(2,1), pos(3,1)


      call VCMPMK(NVCT,POS,PMOM,TVTX,IP)
C
      CALL VCMKVC
      CALL VCMKVX
C
      RETURN
      END
