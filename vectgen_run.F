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
      integer ranpos

      integer lenchr

c     if(iargc().lt.6)then
c        write(6,*)
c    &    'Usage: vectgen run# random_seed_file ZBS_file emin emax nevents [posx] [posy] [posz]'
c        stop
c     endif

      call getarg(1,arg)
      read(arg,*) nrun
      call getarg(2,fseed) ! in/out seed file
      call getarg(3,fname) ! output ZBS file name

c***  Get emin and emax 
      call getarg(4,arg)
      read(arg,*) emin
      call getarg(5,arg)
      read(arg,*) emax
      call getarg(6,arg)
      read(arg,*) num

c***  Get position if given
      if(iargc().gt.7) then
          DO i1= 1, 3
              call getarg(i1 + 6, arg)
              read(arg, *) pos(i1, 1)
          enddo
          ranpos = 0
      else
          ranpos=1
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call set_rflist(21,fname,'LOCAL',' ','WRT',' ',' ',
     &     'recl=5670 status=unknown',' ',' ')
      print *, 'open as LOCAL: ', fname(1:lenchr(fname))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Open vector output ZBS file
      print*,'Trying SKOPENF(21,1,Z,IERR)'
      CALL SKOPENF(21,1,'Z',IERR)

c*** make vector
         DO I=1,num
cccccccccccccccccccccccccccccccccccccccccccccccc
            seed=cnt + i1
            cnt = cnt + 1
            CALL spectrum(9,IDATE,ITIME,pos, seed,emin, emax,ranpos)
cccccccccccccccccccccccccccccccccccccccccccccccc
            call slmcmklow(nrun,1,idate,itime,1)
            iwrite = iwrite+1
            CALL KZWRIT(21)
            CALL KZECLR
         enddo
 100  continue
      close(10)

      call skclosef(21)
      print*,'#subruns: ievent = ',ievent
      print*,'#events/subrun: iwrite = ',iwrite

      STOP
      END




#include "spectrum.F"
