BUILD = ../build
VPATH = .:RCS:$(BUILD)

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h  



LD = gfortran
CF90 = gfortran
CF77 = gfortran

OPT_FFLAGS = -implicitnone -nozero -warn unused -c -O3 -real_size 64 -assume byterecl -assume nobuffered_io -convert big_endian
OPT_FFLAGS = -c -fdefault-real-8

PAR = 
DEBUG_FFLAGS = -fpe0 -traceback -warn unused -CB -u -g -c -real_size 64 -assume byterecl -assume nobuffered_io -convert big_endian
FFLAGS =  $(DEBUG_FFLAGS)
FFLAGS =  $(OPT_FFLAGS)

# run ulimit -c unlimited

F77FLG =  

F90FLG =  

LINKFLAGS =


CPPARCH = 
CPPMODEL  =  


CPPFLAGS = -traditional -P $(CPPARCH) $(CPPMODEL)

#LIBDIR = /local/pgi/linux86/6.0/lib
#LIBDIR = /local/intel/compiler90/ia32/fc/lib   -L /h/geve/lib/linux
#LIBDIR = /h/geve/lib/linux  -L /usr/lib/gcc-lib/i386-redhat-linux/3.2.3 
#LIBDIR = -L /h/geve/lib/linux  -L /usr/lib/gcc/x86_64-redhat-linux/3.4.3/32 -L/h/a152128/lib/linux 

LIBS =  -lg2c  -lblas -llapack -lfftw3 -lecl -lhash -lpthread
LIBS =  -lblas -llapack -lfftw3 


BINDIR = $(HOME)/bin

CPP = /usr/bin/cpp


# Rules for running cpp and updating files in build directory
.H.h:
	rm -f ./$(BUILD)/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.h


.F90.o:
	@rm -f ./$(BUILD)/$*.f90
	@cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.f90
	cd ./$(BUILD) ; $(CF90) $(FFLAGS) $(F90FLG) -o $*.o $*.f90  

.F.o:
	rm -f ./$(BUILD)/$*.f
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS) > ./$(BUILD)/$*.f
	cd ./$(BUILD) ; $(CF77) $(FFLAGS) $(F77FLG) -o $*.o $*.f  



include source.files
TARGET = process.lin


INC2 =$(INC1:.H=.h)
FILES =$(F90FILES) $(F77FILES) $(MODULES)
FFILES =$(F90FILES:.F90=.f90) $(F77FILES:.F=.f) $(MODULES:.F90=.f90)
OBJECTS = $(F90FILES:.F90=.o) $(F77FILES:.F=.o) 
OMOD = $(MODULES:.F90=.o) $(MODULES77:.F=.o)



all: $(TARGET) install


$(TARGET): $(INC2) $(OMOD) $(OBJECTS) 
	cd ./$(BUILD) ; $(LD) $(LINKFLAGS) -o $(TARGET) $(OMOD) $(OBJECTS) $(LIBDIR) $(LIBS) 

install:
	mv $(BUILD)/$(TARGET) $(BINDIR)


clean:
	rm -f $(TARGET)
	cd ./$(BUILD) ; rm -rf *.o *.mod *.f90 rii_files

new: source depend

source:
	mksource.sh > source.files

depend:
	mkdepend_linux.pl | sort -u > depends.file

tags: $(FILES)
	f90tags.sh

include depends.file
