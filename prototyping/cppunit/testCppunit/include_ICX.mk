CC  = icx
CXX = icpx
FC  = ifx
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
