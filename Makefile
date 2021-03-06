###############################################
# 	Makefile for LCLS csPad detector version of myana
#	Now morphed into Cheetah
#	Anton Barty, CFEL, December 2010
#	last updated by Jonas Sellberg on April 26, 2013
#	
#	Note: before this makefile will work it is necessary to compile 
#	the LCLS libraries.  This is done as follows:
#	> cd release/pdsdata
#	> make x86_64-linux-opt
#	(this creates the LCLS libraries necessary for compile to complete.
#	if not done, you'll get an error such as 'can not find libacqdata.so'.
#	Change the target CPU/OS combination if working on a different system, eg: OS-X)
#	
#   if the preprocessor flag CORRELATION_ENABLED is not set, cheetah does not need 
#   the giraffe library, but cannot perform correlation calculations, either
###############################################


#--------------------------------------------------------------
#define variables and directories
TARGET 			= cheetah
#ARCH 			= x86_64-linux
ARCH 			= x86_64-linux-opt

MYANADIR		= myana
CSPADDIR		= cspad-gjw
#LCLSDIR 		= release
LCLSDIR 		= /reg/neh/home/sellberg/source/cheetah/release
HDF5DIR 		= /reg/neh/home/barty/hdf5
TIFFDIR			= /reg/neh/home/feldkamp/tiff
FFTWDIR         	= /reg/neh/home/feldkamp/fftw
GIRAFFEDIR      	= /reg/neh/home/sellberg/source/giraffe
ROOTSYS			= /reg/g/pcds/package/root
#OBJFILES		= main.o XtcRun.o

#INCLUDEDIRS     = -Irelease \

INCLUDEDIRS     = -I$(LCLSDIR) \
                  -I$(HDF5DIR)/include \
                  -I$(TIFFDIR)/include \
                  -I$(FFTWDIR)/include \
                  -I$(GIRAFFEDIR)

#LIBDIRS         = -Lrelease/build/pdsdata/lib/$(ARCH)/ \

LIBDIRS         = -L$(LCLSDIR)/build/pdsdata/lib/$(ARCH)/ \
                  -L$(HDF5DIR)/lib \
                  -L$(TIFFDIR)/lib \
                  -L$(FFTWDIR)/lib \
                  -L$(GIRAFFEDIR)

LIBRARIES       = -lacqdata \
                  -lbld \
                  -lxtcdata \
                  -lopal1kdata \
                  -lcamdata \
                  -lpnccddata \
                  -lcontroldata \
                  -lipimbdata \
                  -lprincetondata \
                  -levrdata \
                  -lencoderdata \
                  -lgsc16aidata \
                  -llusidata \
                  -lcspaddata \
                  -lindexdata \
                  -lanadata \
                  -lcspad2x2data \
                  -lfexampdata \
                  -lphasicsdata \
                  -ltimepixdata \
                  -loceanopticsdata \
                  -lflidata \
                  -lhdf5 \
                  -ltiff \
                  -lpthread \
                  -lfftw3 \
                  -lm \
                  -lgiraffe_static 

CPP             = g++ -c -g
LD              = g++
CPP_LD_FLAGS    = -O4 -Wall
#CFLAGS          = $(INCLUDEDIRS) #no correlation
CFLAGS          = $(INCLUDEDIRS) -DCORRELATION_ENABLED #with correlation and giraffe dependence

LD_FLAGS        = -Wl,-rpath=$(LCLSDIR)/build/pdsdata/lib/$(ARCH)/:$(HDF5DIR)/lib
CFLAGS_ROOT     = $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS_ROOT    = $(shell $(ROOTSYS)/bin/root-config --libs) -Wl,-rpath=$(ROOTSYS)/lib:release/build/pdsdata/lib/$(ARCH)/


all: $(TARGET)


#--------------------------------------------------------------
# Standard myana libraries
$(MYANADIR)/main.o: $(MYANADIR)/main.cc $(MYANADIR)/myana.hh $(MYANADIR)/main.hh
	$(CPP) $(CFLAGS) -o $(MYANADIR)/main.o $<

$(MYANADIR)/XtcRun.o: $(MYANADIR)/XtcRun.cc $(MYANADIR)/XtcRun.hh $(MYANADIR)/main.hh
	$(CPP) $(CFLAGS) -o $(MYANADIR)/XtcRun.o $<
$(MYANADIR)/SplitEventQ.o: $(MYANADIR)/SplitEventQ.cc $(MYANADIR)/SplitEventQ.hh $(MYANADIR)/main.hh
	$(CPP) $(CFLAGS) -o $(MYANADIR)/SplitEventQ.o $<


# csPAD libraries
$(CSPADDIR)/myana_cspad-gjw.o: $(CSPADDIR)/myana_cspad-gjw.cc 
	$(CPP) $(CFLAGS) -o $(CSPADDIR)/myana_cspad-gjw.o $<

$(CSPADDIR)/CspadTemp.o: $(CSPADDIR)/CspadTemp.cc $(CSPADDIR)/CspadTemp.hh
	$(CPP) $(CFLAGS) -o $(CSPADDIR)/CspadTemp.o $<

$(CSPADDIR)/CspadGeometry.o: $(CSPADDIR)/CspadGeometry.cc $(CSPADDIR)/CspadGeometry.hh
	$(CPP) $(CFLAGS) -o $(CSPADDIR)/CspadGeometry.o $<

$(CSPADDIR)/CspadCorrector.o: $(CSPADDIR)/CspadCorrector.cc $(CSPADDIR)/CspadCorrector.hh
	$(CPP) $(CFLAGS) -o $(CSPADDIR)/CspadCorrector.o $<

$(CSPADDIR)/myana_cspad-gjw: $(MYANADIR)/main.o $(MYANADIR)/XtcRun.o \
  $(CSPADDIR)/myana_cspad-gjw.o \
  $(CSPADDIR)/CspadCorrector.o \
  $(CSPADDIR)/CspadGeometry.o \
  $(CSPADDIR)/CspadTemp.o \
  $(CSPADDIR)/myana_cspad-gjw.o
	@echo ""
	@echo "---Linking myana_cspad-gjw---"
	$(LD) $(CPP_LD_FLAGS) $(LD_FLAGS) -o $@ $^ $(LIBDIRS) $(LIBRARIES)
	@echo ""


#--------------------------------------------------------------
# cheetah objects
cheetah.o: cheetah.cpp \
  attenuation.h \
  setup.h \
  worker.h
	$(CPP) $(CFLAGS) $<

worker.o: worker.cpp worker.h \
  background.h \
  commonmode.h \
  correlation.h \
  hitfinder.h \
  worker.h
	$(CPP) $(CFLAGS) $<

setup.o: setup.cpp setup.h \
  attenuation.h \
  data2d.h \
  setup.h \
  worker.h
	$(CPP) $(CFLAGS) $<

data2d.o: data2d.cpp data2d.h 
	$(CPP) $(CFLAGS) $<

commonmode.o: commonmode.cpp commonmode.h \
  setup.h \
  worker.h \
  peakdetect.h \
  point.h
	$(CPP) $(CFLAGS) $<

background.o: background.cpp background.h \
  setup.h \
  worker.h
	$(CPP) $(CFLAGS) $<

hitfinder.o: hitfinder.cpp hitfinder.h \
  setup.h \
  worker.h
	$(CPP) $(CFLAGS) $<

attenuation.o: attenuation.cpp attenuation.h \
  setup.h \
  worker.h
	$(CPP) $(CFLAGS) $<

correlation.o: correlation.cpp correlation.h \
  setup.h \
  worker.h
	$(CPP) $(CFLAGS) $<

peakdetect.o: peakdetect.cpp peakdetect.h \
  pointvector.h \
  point.h
	$(CPP) $(CFLAGS) $<

pointvector.o: pointvector.cpp pointvector.h \
  point.h
	$(CPP) $(CFLAGS) $<

point.o: point.cpp point.h 
	$(CPP) $(CFLAGS) $<


#--------------------------------------------------------------
#compile the different parts of the cheetah and link them together
#the value of ‘$@’ is the target
#the value of ‘$^’ is a list of all the prerequisites of the rule
cheetah: cheetah.o \
  setup.o \
  worker.o \
  data2d.o \
  commonmode.o \
  background.o \
  hitfinder.o \
  attenuation.o \
  correlation.o \
  peakdetect.o \
  pointvector.o \
  point.o \
  $(MYANADIR)/main.o
	@echo ""
	@echo "---Linking the Cheetah---"
    #the order is important: first the objects, then the libraries they need
    #If A uses B, then A must be listed before B in the link command.
	$(LD) $(CPP_LD_FLAGS) $(LD_FLAGS) -o $@ $^ $(LIBDIRS) $(LIBRARIES)
	@echo ""


clean:
	rm -f *.o *.gch myana/*.o $(TARGET)

remake: clean all

.PHONY: all clean remake

# test data
test: cspad_cryst
	./cspad_cryst -f ~gjwillms/cfel-cspad/e55-r0461-s00-c00.xtc -n 2

gdb: cspad_cryst
	gdb ./cspad_cryst -f ~gjwillms/cfel-cspad/e55-r0461-s00-c00.xtc -n 2

valgrind: cspad_cryst
	valgrind ./cspad_cryst -f ~gjwillms/cfel-cspad/e55-r0461-s00-c00.xtc -n 2


