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
      REAL tweight,sdir(3),POS(3,10), emin, emax
      character*200 timefile, fname
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

      if(iargc().lt.3)then
         write(6,*)
     &    'Usage: vectgen run# random_seed_file ZBS_file emin emax [posx] [posy] [posz]'
         stop
      endif

      call getarg(1,arg)
      read(arg,*) istart
      call getarg(2,fseed) ! in/out seed file
      call getarg(3,fname) ! output ZBS file name

c***  Get emin and emax 
      call getarg(4,arg)
      read(arg,*) emin
      call getarg(5,arg)
      read(arg,*) emax

c***  Get position if given
      if(iargc().gt.5) then
          DO i1= 1, 3
              call getarg(i1 + 5, arg)
              read(arg, *) pos(i1, 1)
          enddo
      else
          pos(1, 1) = -10000
      endif



c*** open random seed file, if exist
      inquire(file=fseed, exist=ex)
      if (ex) then
c***  set luxury level 3
         open(11, file=fseed,status='old',form='formatted')
         read(11,*) i1
         close(11)
         call rluxgo(3,i1,0,0)
         cnt = 0

      else
         write(6,*) 'Please set Random seed'
         stop

      endif

      CALL KZINIT
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call set_rflist(21,fname,'LOCAL',' ','WRT',' ',' ',
     &     'recl=5670 status=unknown',' ',' ')
      print *, 'open as LOCAL: ', fname(1:lenchr(fname))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Open vector output ZBS file
      print*,'Trying SKOPENF(21,1,Z,IERR)'
      CALL SKOPENF(21,1,'Z',IERR)

      print*,'Start Read '

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
         if(nrun.lt.istart) goto 10
         if(nrun.gt.istart) goto 100
         ievent = ievent + 1
         if(mod(ievent,10000).eq.0) print*,ievent,' is finish'

c*** make vector ("num" times per minute)
         DO I=1,num
cccccccccccccccccccccccccccccccccccccccccccccccc
            seed=cnt + i1
            cnt = cnt + 1
            CALL spectrum(9,IDATE,ITIME,pos, seed,emin, emax,1)
*     BAHCALL & PINSONNEAULT 2004 (total tmpflux, b8=5.79)
*     Winter et al PRC 73,025503 (2006) (tmpb8 spectrum)
*     BAHCALL @NEUTRIO98 (hep spectrum)
cccccccccccccccccccccccccccccccccccccccccccccccc
            call slmcmklow(nrun,nsub,idate,itime,tweight)
            iwrite = iwrite+1
            CALL KZWRIT(21)
            CALL KZECLR
         enddo
      goto 10

 100  continue
      close(10)

      call skclosef(21)
      print*,'#subruns: ievent = ',ievent
      print*,'#events/subrun: iwrite = ',iwrite

      STOP
      END




#include "spectrum.F"
