CC  = mpiicc
CXX = mpiicpc
FC  = ifort
LINKER = $(CXX)

CFLAGS    = -g -xHost -O3 -qopenmp -std=c99
CXXFLAGS  = $(CFLAGS)
FCFLAGS   =
LFLAGS    = -qopenmp
DEFINES   =
INCLUDES  =
LIBS      =
