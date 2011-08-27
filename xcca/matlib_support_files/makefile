###############################################
# Makefile for matlib support library
###############################################

#--------------------------------------------------------------
#define variables and directories
LIBNAME         = libmat_support.a

INCLUDEDIRS     = -I$(CHEETAHDIR) \
                  -I$(HDF5DIR)/include \
                  -I$(TIFFDIR)/include \
                  -I$(FFTWDIR)/include \
                  -I$(MATLIBDIR)

LIBDIRS         = -L$(HDF5DIR)/lib \
                  -L$(TIFFDIR)/lib \
                  -L$(FFTWDIR)/lib \
                  -L$(MATLIBDIR)/lib


CPP             = g++ -c -g
CFLAGS          = $(INCLUDEDIRS) -O4 -Wall
LD              = g++
LD_FLAGS        = -Wl

all: $(LIBNAME)


#-----------------------------------------------------compile objects
edf.o: edf.cpp edf.h
	$(CPP) $(CFLAGS) $<

edfkeywords.o: edfkeywords.cpp edfkeywords.h
	$(CPP) $(CFLAGS) $<

timer.o: timer.cpp timer.h
	$(CPP) $(CFLAGS) $<

util.o: util.cpp util.h
	$(CPP) $(CFLAGS) $<


#--------------------------------------------------------------
#link the different parts
#the value of ‘$@’ is the target
#the value of ‘$^’ is a list of all the prerequisites of the rule
$(LIBNAME): edf.o \
  edfkeywords.o \
  timer.o \
  util.o
	@echo ""
	@echo "---Building library " $(LIBNAME) " ---"
	ar ru $(LIBNAME) $^
	ranlib $(LIBNAME) 
	@echo ""


clean:
	rm -f *.o *.gch $(LIBNAME)

remake: clean all

.PHONY: all clean remake

