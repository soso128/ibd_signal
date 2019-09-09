************************************************************************
*     -------------------------------------
      SUBROUTINE SKROOT_SET_MCINFO(lun)
*     -------------------------------------
*
*     (Purpose)
*        Set MC branch for skroot file
*
*     (Input)
*        lun : Logical Unit Number
*
*     (Output)
*       
*     (Return)
*
*     (Creation Date and Author)
*       2009.01.15 ; First version by Y.Takeuchi
*       2018.03.03 added LOWMC bank information  by Y.Takeuchi
*          call this routine before "skroot_set_tree()"
*          to save nrunsk/nsubsk/ndaysk/ntimsk from LOWMC
************************************************************************
      IMPLICIT NONE

#include "skhead.h"
#include "skparm.h"
#include "vcwork.h"
#include "vcvrtx.h"

      integer lun
      integer ivmcp,ivtap,ivabl,ivrbw,ivrpm,ivrys,lux
      integer iseed_ranlux,k1_ranlux,k2_ranlux
      integer iseed_rm48,k1_rm48,k2_rm48,iseed1,iseed2,lux2,iseed,k1,k2
      real    energy(100)
      real    darkds,qresds,tresds,twinds,trigds,gateds,beftds,deadds,sigwds
      integer mcrun, mcninfo
cc      real    mcinfo(255)
      integer   mcinfo(255) ! fixed by y.t (2018/01/26)

      integer ipbank, ndata, idata(100)
      real    rdata(100)

#ifdef NEUTRON
      integer i,j
      equivalence (idata(1), rdata(1))
#include "ncapgam.h"
#endif

c*** read MCHEAD
      CALL KZBLOC('MCHEAD',IPBANK)
      IF (IPBANK.ne.0) then
         CALL KZGET1('MCHEAD',0,NDATA,IDATA)
         IF (NDATA.gt.0) then
            IVMCP=IDATA(1)
            IVTAP=IDATA(2)
            IVABL=IDATA(3)
            IVRBW=IDATA(4)
            IVRPM=IDATA(5)
            IVRYS=IDATA(6)
            LUX=IDATA(7)
            ISEED_RANLUX=IDATA(8)
            K1_RANLUX=IDATA(9)
            K2_RANLUX=IDATA(10)
            ISEED_RM48=IDATA(11)
            K1_RM48=IDATA(12)
            K2_RM48=IDATA(13)
         endif
      endif

c*** read MCVECT & MCVERTEX
      call VCRDVCCM()

c*** read MCPARMCONV
      CALL KZBLOC('MCPARMCONV',IPBANK)
      IF (IPBANK.ne.0) then
         CALL KZGET1('MCPARMCONV',0,NDATA,RDATA)
          IF (NDATA.gt.0) then
            DARKDS = RDATA(1)
            TCNVSK = RDATA(2)
            QCNVSK = RDATA(3)
            TRESDS = RDATA(4)
            QRESDS = RDATA(5)
            DTHRSK = RDATA(6)
            TWINDS = RDATA(7)
            TRIGDS = RDATA(8)
            GATEDS = RDATA(9)
            BEFTDS = RDATA(10)
            DEADDS = RDATA(11)
            SIGWDS = RDATA(12)
         endif
      endif

#ifdef NEUTRON
*** Temporary add by Y.Koshio (for neutron capture)
      mcninfo = 0
#endif

c*** read LOWMC
      CALL KZBLOC('LOWMC',IPBANK)
      IF (IPBANK.ne.0) then
         CALL KZGET1('LOWMC',1,NDATA,IDATA)
         IF (NDATA.gt.0) then
            mcinfo(1) = IDATA(1) ! nrunsk 
            mcinfo(2) = IDATA(2) ! nsubsk
            mcinfo(3) = IDATA(3) ! ID(1)*10000+ID(2)*100+ID(3)
            mcinfo(4) = IDATA(4) ! IT(1)*10000+IT(2)*100+IT(3)
            mcinfo(5) = IDATA(5) ! tweight
            if (mcninfo.lt.5) mcninfo = 5

            nrunsk = mcinfo(1)
            nsubsk = mcinfo(2)
            ndaysk(1) = idata(3)/10000 
            ndaysk(2) = (idata(3)-ndaysk(1)*10000)/100
            ndaysk(3) = idata(3)-ndaysk(1)*10000-ndaysk(2)*100
            ntimsk(1) = idata(4)/10000 
            ntimsk(2) = (idata(4)-ntimsk(1)*10000)/100
            ntimsk(3) = idata(4)-ntimsk(1)*10000-ntimsk(2)*100
         endif
      endif
     
#ifdef NEUTRON
*** Temporary add by Y.Koshio (for neutron capture)
      CALL KZBLOC('NCAPGAM',IPBANK)
      IF (IPBANK.ne.0) then
         CALL KZGET2('NCAPGAM',0,0,ndata,idata)
         num_cap = idata(1)
         num_gam = idata(2)
         if(num_gam .ne. 1) then
            print*,'Num of gamma should be 1, but ',num_gam
            print*,'Call Koshio'
            stop
         endif
         print*,'Num of neutron capture ',num_cap
         mcninfo = mcninfo + num_cap*4

         do i = 1, num_cap
            do j = 1, num_gam
               CALL KZGET2('NCAPGAM',i,j,ndata,rdata)

c               ip_cap(i) = idata(1)
c               en_cap(i) = rdata(2)
c               pos_cap(1,i) = rdata(3)
c               pos_cap(2,i) = rdata(4)
c               pos_cap(3,i) = rdata(5)
c               t_cap(i) = rdata(6)
c               en_gam(i,j) = rdata(7)
c               dir_gam(1,i,j) = rdata(8)
c               dir_gam(2,i,j) = rdata(9)
c               dir_gam(3,i,j) = rdata(10)

               mcinfo((i-1)*4+1) = idata(3)
               mcinfo((i-1)*4+2) = idata(4)
               mcinfo((i-1)*4+3) = idata(5)
               mcinfo((i-1)*4+4) = idata(6)
            enddo
         enddo
      endif


c      print*, num_cap, num_gam

#endif


c*** set for skroot
      call skroot_set_mc(lun,ivmcp,ivtap,ivabl,ivrbw,ivrpm,ivrys,lux,
     $     iseed_ranlux,k1_ranlux,k2_ranlux,
     $     iseed_rm48,k1_rm48,k2_rm48,iseed1,iseed2,lux2,iseed,k1,k2,
     $     nvtxvc,pvtxvc,iflvvc,iparvc,timvvc,
     $     nvc,ipvc,energy,pvc,iorgvc,ivtivc,icrnvc,iflgvc,ivtfvc,
     $     darkds,tcnvsk,qcnvsk,tresds,qresds,dthrsk,twinds,trigds,
     $     gateds,beftds,deadds,sigwds,
     $     mcrun,mcninfo,mcinfo)
      
      RETURN
      END
