***********************************************************************
*     ------------------------------------------------------
*     Sample program for generating vector of solar neutrino
*     ------------------------------------------------------
*
*     (Creation Date and Author)
*       1995.01.19 ; First version by K.Inoue
*            11.16 ;      modified by Y.Koshio
*       1997.06.21 ;      modified by Y.Takeuchi for real-time job
*       1998.04.06 ;      modified by Y.Takeuchi add RDNOUT() and RDMIN()
*       2004.03.16 ;      modified by Y.Takeuchi 
*           for rlu(), add SK_GEOMETRY = 2
*           http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/v115/top.html
*       2007.10.25 ;      modified by Byeoungs Yang for SK-3
*           change SK_GEOMETRY into 3
*           change FLAG of SLVECTB8 into 4
*	2017.03.15 ; 	  modified by P.Weatherly to use 
*	    integer irnsk4 = 61525 initial run SK-IV
* 	    integer lrnsk4 = 75765 last run in SK-IV for update (2657days)
*       2017.05.25 ; 	  modified by y.takeuchi
*           added "set_rflist", etc.     
***********************************************************************
      implicit none

#include "skhead.h"

      INTEGER IERR,nrun,nsub,IDATE(3),ITIME(3),istat,i,ievent
      INTEGER iwrite,lfbadrun
      INTEGER ibad_subrun(61525:77958,3000)  !2657days
      integer idate_prev(3), itime_prev(3)
      REAL tweight,sdir(3),POS(3,10)
      character*200 timefile, fname,fnamespec
      integer iargc,istart,iend
      character*80 arg, fseed
      logical*4 ex
      integer i1, cnt
      integer*4 seed
      integer irnsk4, lrnsk4   !Added by P. Weatherly
      integer num,neutron
      parameter(num = 10)

      integer lenchr

Cpi   SET the initial run for SK-IV and the Last Run for the most recent update
      irnsk4 = 61525	! SK-IV Start Run
      lrnsk4 = 77958	! SK-IV Last  Run for update days 2400
      timefile=
     $'/home/sklowe/realtime_sk4_rep/solar_apr19/timevent/timevent.r061525.r077958'

      SK_GEOMETRY = 4




c*** init. multiple badrun 
     
      do nrun = irnsk4, lrnsk4
         do nsub = 1, 3000
            ibad_subrun(nrun, nsub) = 0
         enddo
      enddo

c*** open time event file
      print*,'Open timefile:',timefile
      open(10, file=timefile,
     $     status='old',form='unformatted')

      print*,'timefile opened'


 10   continue
**** read run number, subrun number, date(yy,mm,dd), time(hh,mm,ss), run time(sec), solar direction(direction cosine toward(?) sun), istat
         call slredtimev(10,nrun,nsub,idate,itime,tweight,sdir,istat)
      if(istat.ne.0) then
         if(nrun.ge.lrnsk4) goto 100
         print *, 'reading error! stop'
         stop
      endif
 
         if( nrun .lt. 60000 )then
            nrun = nrun + 65536         ! for SK-IV analysis
         endif

         if(nrun.lt.irnsk4) goto 10
         if(nrun.gt.lrnsk4) goto 100


c*** check duplicate date
      if (idate(1).eq.idate_prev(1) .and.
     $     idate(2).eq.idate_prev(2) .and.
     $     idate(3).eq.idate_prev(3) .and.
     $     itime(1).eq.itime_prev(1) .and.
     $     itime(2).eq.itime_prev(2) .and.
     $     itime(3).eq.itime_prev(3)) then

         if (nsub.gt.0 .and. nsub.lt.3000) then
            ibad_subrun(nrun, nsub) = 1
         else
            print *, 'data error!', nrun, nsub,idate,itime
            stop
         endif
      endif

c*** save previous date
      idate_prev(1) = idate(1)
      idate_prev(2) = idate(2)
      idate_prev(3) = idate(3)
      itime_prev(1) = itime(1)
      itime_prev(2) = itime(2)
      itime_prev(3) = itime(3)


c***  badsub run cut 
      if( lfbadrun(nrun, nsub).eq.1) then 
c            print *, "Bad run!"
            goto 10
         endif
         
c** multiple badsub run cut
         if(nrun.ge.irnsk4 .and. nrun.le.lrnsk4) then
            if (nsub.gt.0 .and. nsub.lt.3000) then
               if (ibad_subrun(nrun, nsub).ne.0) then
c                  print *, "Bad subrun!"
                  goto 10
                  endif
            endif
         endif

         if(istat.ne.0) goto 100	! 0 : OK, -1 : EOF
         print *, nrun, nsub, idate, itime, tweight
      goto 10

 100  continue
      close(10)

      call skclosef(21)
      print*,'#subruns: ievent = ',ievent
      print*,'#events/subrun: iwrite = ',iwrite

      STOP
      END
