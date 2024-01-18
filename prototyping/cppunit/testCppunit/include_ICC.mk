CC  = icc
CXX = icpc
FC  = ifort
LINKER = $(CXX)

CFLAGS    = -g -O3
CXXFLAGS  = $(CFLAGS)
FCFLAGS   = 
LFLAGS    =
DEFINES   =

INCLUDES =
INCLUDES += -I ../cppunit/include

LIBS =
LIBS += -L../cppunit/lib -lcppunit