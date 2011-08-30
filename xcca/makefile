###############################################
# Makefile for standalone XCCA correlation software
###############################################

#--------------------------------------------------------------
#define variables and directories
TARGET 			= xcca

HDF5DIR         = /reg/neh/home/feldkamp/hdf5
TIFFDIR         = /reg/neh/home/feldkamp/tiff
FFTWDIR         = /reg/neh/home/feldkamp/fftw
#HDF5DIR         = /usr/local/hdf5
#TIFFDIR         = /usr/local
#FFTWDIR         = /usr/local
MATLIBDIR       = ./matlib_support_files
CHEETAHDIR      = ..

INCLUDEDIRS     = -I$(CHEETAHDIR) \
                  -I$(HDF5DIR)/include \
                  -I$(TIFFDIR)/include \
                  -I$(FFTWDIR)/include \
                  -I$(MATLIBDIR)

LIBDIRS         = -L$(HDF5DIR)/lib \
                  -L$(TIFFDIR)/lib \
                  -L$(FFTWDIR)/lib \
                  -L$(MATLIBDIR)

LIBRARIES       = -lhdf5 \
                  -ltiff \
                  -lpthread \
                  -lfftw3 \
                  -lm \
                  -lmat_support

CPP             = g++ -c -g
CPP_LD_FLAGS    = -O4 -Wall
CFLAGS          = $(INCLUDEDIRS)
LD              = g++
LD_FLAGS        = -Wl,-rpath=$(HDF5DIR)/lib

all: $(TARGET)

#-----------------------------------------------------compile objects
main.o: main.cpp \
  analyze.h \
  ../crosscorrelator.h \
  ../arrayclasses.h \
  ../arraydataIO.h
	$(CPP) $(CFLAGS) $<

analyze.o: analyze.cpp analyze.h \
  ../crosscorrelator.h \
  ../arrayclasses.h \
  ../arraydataIO.h
	$(CPP) $(CFLAGS) $<

crosscorrelator.o: ../crosscorrelator.cpp ../crosscorrelator.h \
  ../arrayclasses.h \
  ../arraydataIO.h
	$(CPP) $(CFLAGS) $<

arrayclasses.o: ../arrayclasses.cpp ../arrayclasses.h
	$(CPP) $(CFLAGS) $<

arraydataIO.o: ../arraydataIO.cpp ../arraydataIO.h \
  ../arrayclasses.h
	$(CPP) $(CFLAGS) $<

fouriertransformer.o: ../fouriertransformer.cpp ../fouriertransformer.h \
  ../arrayclasses.h
	$(CPP) $(CFLAGS) $<



#--------------------------------------------------------------
#link the different parts
#the value of ‘$@’ is the target
#the value of ‘$^’ is a list of all the prerequisites of the rule
$(TARGET): main.o \
  analyze.o \
  crosscorrelator.o \
  arrayclasses.o \
  arraydataIO.o \
  fouriertransformer.o
	@echo ""
	@echo "---Linking " $(TARGET) " ---"
    #the order is important: first the objects, then the libraries they need
    #If A uses B, then A must be listed before B in the link command.
	$(LD) $(CPP_LD_FLAGS) $(LD_FLAGS) -o $@ $^ $(LIBDIRS) $(LIBRARIES)
	@echo ""


clean:
	rm -f *.o *.gch $(TARGET)

remake: clean all

.PHONY: all clean remake

