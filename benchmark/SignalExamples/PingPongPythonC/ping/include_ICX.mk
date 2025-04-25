CC  = mpiicx
CXX = mpiicpx #If static EMPIRE_API is used mpifort has to be used for linking
FC  = mpiifx
LINKER = $(CXX)

CFLAGS    = -g -O3
CXXFLAGS  = $(CFLAGS)
FCFLAGS   = 
LFLAGS    = 
DEFINES   =
INCLUDES = -I$(EMPIRE_API_INC_ON_MACHINE)
LIBS      = $(EMPIRE_API_LIBSO_ON_MACHINE)


 


