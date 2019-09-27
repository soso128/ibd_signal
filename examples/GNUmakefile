SKOFL_ROOT = /home/elhedri/skofl/
include $(SKOFL_ROOT)/config.gmk
FLAGS=-O3 -DNDEBUG -pipe -m64 -Wshadow -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated-declarations -pthread -MMD -MP
FFLAGS=-O3 -DNDEBUG -finit-local-zero -fno-backslash -fno-automatic -fno-second-underscore -ffixed-line-length-132 -fPIC

LOCAL_INCS=-I. -I.. -I$(ROOTSYS)/include/

OBJS=read_reweight.o dsigma read_reweight

LIBSROOT=-L$(ROOTSYS)/lib/ -lCint -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lGui

LOCAL_LIBS= -L. -L/usr/lib64/ $(LIBSROOT) -lm -ldl -rdynamic
#LIBS= -L. $(LIBSROOT) -lNet -lm -ldl -rdynamic

all: $(OBJS)

dsigma: ../dsigma.h ../dsigma.cpp
	g++ -c ../dsigma.cpp

read_reweight.o: ../dsigma.h read_reweight.cpp
	g++ -c -O3 -g -Wall -I$(SKOFL_INCDIR) $(LOCAL_INCS) read_reweight.cpp

read_reweight:  dsigma.o read_reweight.o
	LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR) $(CXX) $(CXXFLAGS) -o read_reweight dsigma.o read_reweight.o $(LDLIBS) $(LIBSROOT)
