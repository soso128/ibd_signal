
#SKOFL_ROOT = /home/sklowe/skofl/r17806 // still have problem used default version skofl_11c with special slmkb8.F
#SKOFL_ROOT=/usr/local/sklib_g77/skofl_14b
#SKOFL_ROOT = /usr/local/sklib_gcc4.8.2/skofl-trunk

SKOFL_ROOT = /home/elhedri/skofl/
#SKOFL_ROOT = /usr/local/sklib_gcc4.8.5/skofl-trunk
#SKOFL_ROOT = /home/sklowe/skofl/r25242/

include $(SKOFL_ROOT)/config.gmk
#FORTRANDEFINES += -DRUNNUM_REP=20000

LOCAL_INC	= 

LOCAL_LIBS	= -lsollib_4.0 -lsklowe_7.0 -lsollib_4.0 -lwtlib_5.1 -lbonsai_3.3 -lstmu -lska

#
#  Objects
#

OBJS   =  dsigma vectgen make_random lowfit_sk4_mc incorporate

all: $(OBJS)
	$(RM) *.o *~

dsigma: dsigma.h dsigma.cpp
	g++ -c dsigma.cpp

vectgen: vectgen.o dsigma.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o vectgen vectgen.o dsigma.o $(LDLIBS) 

incorporate: initialize.o incorporate.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o incorporate incorporate.o initialize.o $(LDLIBS) 

clean: 
	$(RM) *.o *~ core fort.* $(OBJS)

install.exec: 


