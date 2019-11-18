
#SKOFL_ROOT = /home/sklowe/skofl/r17806 // still have problem used default version skofl_11c with special slmkb8.F
#SKOFL_ROOT=/usr/local/sklib_g77/skofl_14b
#SKOFL_ROOT = /usr/local/sklib_gcc4.8.2/skofl-trunk

SKOFL_ROOT = /home/elhedri/skofl/
#SKOFL_ROOT = /usr/local/sklib_gcc4.8.5/skofl-trunk
#SKOFL_ROOT = /home/sklowe/skofl/r25242/

include $(SKOFL_ROOT)/config.gmk
#FORTRANDEFINES += -DRUNNUM_REP=20000

LOCAL_INC	= 

LOCAL_LIBS	=   zbsinit.o set_rflistC.o fort_fopen.o \
		   $(APLIB) $(APLIB) \
	   	   -lsollib_4.0 -lsklowe_7.0 -lsollib_4.0 -lwtlib_5.1 -lbonsai_3.3 -lstmu -lska \
		-L/home/iida/relic/lib -lrelic\
		-L$(NEUT_ROOT)/lib/Linux_pc\
		-lfccomb -lprotonid -lpolfit -lexpq -lringlib -lseplib -lprtlib -lmuelib -lstmu -lffit -laplowe \
	        -lodlib  -laplib -lmslib -lmsfit -ltf -lThreeProb -lnumrec -lneutflux -lska -lfiTQun \
		-lmsfit  -lprtlib  \
                -lstmu -lapdrlib  \
                -lConnectionTableReader 

LIBSROOT=-L$(ROOTSYS)/lib/ -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

INCROOT=-I$(ROOTSYS)/include/

FORTRANINCLUDES = $(SITE_INCLUDES) -I. -I$(SKOFL_FORTRAN_INCDIR) -I$(SKOFL_FORTRAN_INCDIR)/lowe -I$(A_FORTRAN_INCDIR) 

#
#  Objects
#

OBJS   =  zbsinit.o set_rflistC.o fort_fopen.o dsigma vectgen vectgen_run make_random incorporate lem reweight fortran_interface_skmc.o skmc_manager.o lowfit_sk4_mc.o lowfit_sk4_mc leaf

all: $(OBJS)

lowfit_sk4_mc: lowfit_sk4_mc.o skmc_manager.o fortran_interface_skmc.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o lowfit_sk4_mc lowfit_sk4_mc.o skmc_manager.o fortran_interface_skmc.o $(INCROOT) -I$(SKOFL_INCDIR) -L$(SKOFL_LIBDIR) $(LIBSROOT) $(LDLIBS)

lem: lem.cpp init_geom.o
	g++ -g -O3 -o lem lem.cpp init_geom.o $(INCROOT) -I$(SKOFL_INCDIR) -L$(SKOFL_LIBDIR) $(LIBSROOT) $(LDLIBS)

dsigma: dsigma.h dsigma.cpp
	g++ -c -g dsigma.cpp

reweight: reweight.cpp dsigma.o
	g++ -g -O3 -o reweight reweight.cpp dsigma.o $(INCROOT) -I$(SKOFL_INCDIR) -L$(SKOFL_LIBDIR) $(LIBSROOT) $(LDLIBS)

vectgen: vectgen.o dsigma.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o vectgen vectgen.o dsigma.o $(LDLIBS) 

vectgen_run: vectgen_run.o dsigma.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o vectgen_run vectgen_run.o dsigma.o $(LDLIBS) 

incorporate: zbsinit.o set_rflistC.o fort_fopen.o initialize.o incorporate.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o incorporate incorporate.o initialize.o $(LDLIBS) 

clean: 
	$(RM) *.o *~ core fort.* $(OBJS)

install.exec: 


