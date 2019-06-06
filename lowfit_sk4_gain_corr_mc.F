c
c     lowfit_sk4_gain_corr.F         13-APR-2017     Y.Takeuchi
c
c     for SK-IV process apply lowe fitters with gain correction
c
c   output: (arg1)   skroot file 
c   input:  (arg2-)  skroot reformat files 
c
c
      implicit none

#include "skhead.h"
#include "skparm.h"
#include "sktq.h"
#include "skbadc.h"
#include "geopmt.h"
#include "skday.h"
#include "skwt.h"
#include "skroot_lowe.h"

      integer iret, lfflag, ipbank, i, lenchr, nread, nevent,ierr
      real    watert, wallsk
      integer icall, isave_event, lastrun, log_level
      real dark, darklf, rrr
      integer nsub_sav, nrun_sav, sdrun
      
      integer lun
      parameter (lun=10)

      integer iargc, narg
      character*200 fname_in
      character*200 fname_out
      character*200 fname_base

c*** number of hit limit for clusfit
      integer NHITCUT
      parameter (NHITCUT =   800)

c*** NTUPLE parameters
      character*200 fname_nt, mom
      real effwallf, forwallf

      INTEGER  NH, iquest
      REAL     H
      PARAMETER (NH=10000000)
      COMMON /PAWC/ H(NH)
      COMMON /QUEST/ IQUEST(100)

      integer istat,icycle, ntnum
      parameter(ntnum=62)
      real data(ntnum)
      character*8 nttag(ntnum)
      data nttag /'nrun', 'ene0', 
     $     'x0','y0','z0','dx0','dy0','dz0','energy0',
     $     'wall0','effwal0',
     $     'clx','cly','clz','clgdn','cldx','cldy','cldz','cln50',
     $     'clwall','cleffwal','clneff70','clenergy',
     $     'bsx','bsy','bsz','bsgdn','bsdx','bsdy','bsdz',
     $     'bswall','bseffwal','bsneff70','bsenergy',
     $     'clcossun','bscossun','sdz','bsforw',
     $     'bspatlik', 'clpatlik', 'bsneff80', 'bsneff90',
     $     'bsnef100','bsnef110','bsnef120','bsnef130',
     $     'bsnef140','bsnef150','bsnef160','bsneff',
     $     'bsdx2','bsdy2','bsdz2','bsene2',
     $     'bsovaq','bsclik','bseffwal','rr','elapday',
     $     'bsdirks', 'poswal', 'rrr'/
cccccccccccccccccccccccccccccccccccccccccccccccccc

c*** mc info for SKROOT files
#include "vcwork.h"
#include "vcvrtx.h"
      integer ivmcp,ivtap,ivabl,ivrbw,ivrpm,ivrys,lux
      integer iseed_ranlux,k1_ranlux,k2_ranlux
      integer iseed_rm48,k1_rm48,k2_rm48,iseed1,iseed2,lux2,iseed,k1,k2
      real    energy(100)
      real    darkds,qresds,tresds,twinds,trigds,gateds,beftds,deadds,sigwds
      integer mcrun, mcninfo
      integer mcinfo(255)

      integer elapday
      real rr
      
c*** check arguments
      narg = iargc()
      if (narg .lt. 2) then
         print *, 'Usage: mcfit bad_run_number input out_base'
         call exit(1)
      endif
c      call getarg(1,fname_in)
c      read(fname_in,*) run_num

      call getarg(2,fname_in)
      call getarg(3,fname_base)
      print *,fname_in, fname_base
      fname_nt   = fname_base(1:lenchr(fname_base)) // '.mcfit.nt'
      fname_out  = fname_base(1:lenchr(fname_base)) // '.mcfit.root'

c*** get output file name
      call skroot_open(lun,fname_out(1:lnblnk(fname_out)))

c*** prepare input files
      call skroot_set_input_file(lun,fname_in(1:lnblnk(fname_in)))
c     call set_rflist(10,fname_in(1:lnblnk(fname_in)),'LOCAL',' ','RED',' ',' ',
c    &     'recl=5670 status=old',' ',' ')
c     call skopenf(10,1,'z',ierr)
c     goto 9999
c     call skoptn('31,30')

c*** initilalize skroot files
      call skroot_init(lun)

c*** set SK_GEOMETRY
      SK_GEOMETRY = 4

c*** init.
      call kzinit
      call geoset
      call skbadopt(23)         ! for lowe events 
      call SKOPTN("31,30,26,25")

c***  init water transparency table
      call skrunday 
      call skwt_gain_corr 

c***  set run number 
c      call skbadch(nrun, 1, iret)
c      call lfwater(relapse(nrun), watert)
c      print *, 'nrun, watert, nbad0 = ', nrun, watert, nbad0
      
ccccccccccccccccccccccccccccccccccccccccccccccccc
c   1) Lots of output
c   2) only prints filenames
c   3) only prints when not found
c   4) do not print
c   others) only prints when found 

c*** suppress outputs from skbadch
      log_level_skbadch = 4
c*** suppress outputs from skdark
cc      log_level_skdark = 4
c*** suppress outputs from lfallfit_sk4_gain_corr
cc      log_level = 4
ccccccccccccccccccccccccccccccccccccccccccccccccc


c*** initialize NTUPLE
      call hlimit(-nh)
      print *, fname_nt
      IQUEST(10) = 65000
      call HROPEN(96, 'ntuple', fname_nt, 'N', 8000, istat)
      call HBOOKN(1, 'ntuple', ntnum, '//ntuple', 10000, nttag)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c*** main loop 
 10   CONTINUE 
         isave_event = 0
         CALL SKREAD(lun,*1002,*1001,*10,*10)
c         call skgetv
c         print *, mcninfo,mcinfo(1),mcinfo(2),mcinfo(3),mcinfo(4),mcinfo(5)
         
         nread = nread + 1
         if(mod(nread,10000).eq.0) print *,nread,nrunsk,nevsk,nqisk

cc         if(nrunsk.eq.sdrun) goto 10

c***  init. bonsai 
         if (icall.eq.0) then
            CALL cfbsinit(MAXPM,xyzpm)
            icall = 1
         endif

c*** load all branches
         call skroot_get_entry(lun)

c*** update dark rate, watert, badch
         if(nsubsk.ne.nsub_sav .or. nrunsk.ne.nrun_sav) then
            if(nrunsk.ne.nrun_sav) then
c***** clear variables
               call lfclear_all()
               dark = darklf(nrunsk)
c***** to skip  strange badrun 
c               if (  dark.lt.2000.0 .or. dark .gt. 20000.0  ) then
c                  sdrun = nrunsk
c                  write(*,*) '*************************************'
c                  write(*,*) 'Run number     :',nrunsk
c                  write(*,*) 'darklf =  ', dark, 'is very strage.'
c                  write(*,*) 'skip this run', nrunsk
c                  Write(*,*) '*************************************'
c                  GOTO 10
c               else
                  print *, 'nrunsk, dark = ', nrunsk, dark
c               endif

               call lfwater(relapse(nrunsk), watert)
               print *, 'nrunsk, watert = ', nrunsk, watert
            endif
            call skbadch(nrunsk,nsubsk,ierr)
cc            print *, 'nrunsk, nsubsk, nbad, nbad0 = ', nrunsk, nsubsk, nbad0, nbad
            nrun_sav = nrunsk
            nsub_sav = nsubsk
         endif

c*** apply lowfit
         if (nqisk.lt.NHITCUT) then
            nevent = nevent + 1
            call lfallfit_sk4_gain_corr_mc(watert, NHITCUT, 0, log_level, lfflag)
         endif

c         print *,'==============='
c         print *, nrunsk, nevsk, posmc
c         print *, dirmc
c         print *, pabsmc
c         print *, darkmc
c         print *, energymc

cccccccccccccccccccccccccccccccccccccccccc
c         print *,'===MC=== ',
c     $        nrunsk, posmc, dirmc(1,1),dirmc(2,1),dirmc(3,1),
c     $        pabsmc(1),energymc(1)
cccccccccccccccccccccccccccccccccccccccccc

         rrr = sqrt(poswal(1)**2 + poswal(2)**2)

c***  store LOWE branch
         call skroot_set_lowe(lun,bsvertex, bsresult, bsdir, 
     $        bsgood, bsdirks, bseffhit, bsenergy, bsn50, bscossun,
     $        clvertex, clresult, cldir, clgoodness,     
     $        cldirks, cleffhit, clenergy, cln50, clcossun,
     $        latmnum, latmh, lmx24, ltimediff, lnsratio, lsdir, 
     $        spaevnum, spaloglike, sparesq, spadt, spadll, spadlt,
     $        spamuyn, spamugdn,
     $        posmc, dirmc, pabsmc, energymc, darkmc, islekeep,
     $        bspatlik, clpatlik, lwatert, lninfo, linfo)
               
c*** store header & TQ info.
         call skroot_set_tree(lun)
c*** output root file 
         call skroot_fill_tree(lun) ! write to disk

         elapday= relapse(nrunsk)
         rr = (bsvertex(1)**2+bsvertex(2)**2)/10000.

c***  fill ntuple
         data(1) = nrunsk
         data(2) = energymc(2)   ! neutrino energy for b8mc
         data(3) = posmc(1)
         data(4) = posmc(2)
         data(5) = posmc(3)
         data(6) = dirmc(1,1)
         data(7) = dirmc(2,1)
         data(8) = dirmc(3,1)
         data(9) = energymc(1)   ! electron generated energy
         data(10) = wallsk(posmc)
         data(11) = effwallf(1, posmc, data(6), poswal)

         data(12) = clvertex(1)
         data(13) = clvertex(2)
         data(14) = clvertex(3)
         data(15) = clgoodness
         data(16) = cldir(1)
         data(17) = cldir(2)
         data(18) = cldir(3)
         data(19) = cln50
         data(20) = wallsk(clvertex)
         data(21) = cleffwal
         data(22) = cleffhit(1)   ! 70m
         data(23) = clenergy 

         data(24) = bsvertex(1)
         data(25) = bsvertex(2)
         data(26) = bsvertex(3)
         data(27) = bsgood(2)
         data(28) = bsdir(1)
         data(29) = bsdir(2)
         data(30) = bsdir(3)
         data(31) = wallsk(bsvertex)
         if (bsdir(1) .lt. 1.1) then
            data(32) = bseffwal
         else
            data(32) = -999.
         endif
         data(33) = bseffhit(1) ! 70m
         data(34) = bsenergy 

         data(35) = clcossun
         data(36) = bscossun
         data(37) = lsdir(3)   ! solar direction Z
         data(38) = forwallf(bsvertex, bsdir)
         data(39) = bspatlik
         data(40) = clpatlik
         data(41) = bseffhit(2)  ! 80m
         data(42) = bseffhit(3)  ! 90m
         data(43) = bseffhit(4)  ! 100m
         data(44) = bseffhit(5)  ! 110m
         data(45) = bseffhit(6)  ! 120m
         data(46) = bseffhit(7)  ! 130m
         data(47) = bseffhit(8)  ! 140m
         data(48) = bseffhit(9)  ! 150m
         data(49) = bseffhit(10) ! 160m
         data(50) = bseffh       ! at watert
         data(51) = bsdir_lfdir2(1)
         data(52) = bsdir_lfdir2(2)
         data(53) = bsdir_lfdir2(3)
         data(54) = bsenergy_lfdir2 
         data(55) = bsovaq
         data(56) = bsclik
         data(57) = bseffwal
         data(58) = rr
         data(59) = elapday
         data(60) = bsdirks
         data(61) = poswal(3)
         data(62) = rrr
         call HFN(1, data)

      goto 10

 1003 continue
      write(0,*)'read error2'
      stop
 1002 continue
      write(0,*)'read error'
      goto 10
 1001 continue

 9999 continue
c*** close output file
      call skroot_close(lun)
      call skroot_end

c*** for ntuple
      call HROUT(1, icycle, ' ')
      call HREND('ntuple')

c*** for bonsai 
      call cfbsexit
      print *, 'done. nrun/nsub/nread/nevent=', 
     $     nrunsk, nsubsk, nread, nevent

      stop
      end
      
