
      implicit none

      integer i,ivec(25)
      real rlu
      character*80 fseed
      logical*4 ex

      if(iargc().ne.1)then
         write(6,*) 'Usage: make_random random_seed_file'
         stop
      endif
      call getarg(1,fseed)


      inquire(file='seed.dat', exist=ex)
      if (ex) then
c***  set luxury level 3
         call rluxgo(3, 0, 0, 0)
         open(11, file='seed.dat',status='old',form='formatted')
         read(11,*) ivec
         CALL RLUXIN(IVEC)
         close(11)
      else
c***  set luxury level 3 and initial seed
         call rluxgo(3, 057852121, 0, 0)


      endif

      i=int(rlu()*100000000)
      open(20,file=fseed,status='new',form='formatted')
      write(20,*) i
      close(20)

c*** open random seed file, if exist
      open(11, file='seed.dat',status='unknown',form='formatted')
      CALL RLUXUT(IVEC)
      write(11,*) ivec
      close(11)

      STOP
      END
