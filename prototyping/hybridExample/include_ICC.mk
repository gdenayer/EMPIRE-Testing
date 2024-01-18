CC  = mpiicc
CXX = mpiicpc
FC  = ifort
LINKER = $(CXX)

CFLAGS    = -g -xHost -O3 -openmp -std=c99
CXXFLAGS  = $(CFLAGS)
FCFLAGS   = 
LFLAGS    = -openmp
DEFINES   =
INCLUDES  = 
LIBS      = 


 


