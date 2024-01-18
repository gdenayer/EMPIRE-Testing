CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

CFLAGS    = -g -O3 -Wall -std=c99 -Wno-write-strings
CXXFLAGS  = -g -O3 -Wall -Wno-write-strings
FCFLAGS   = 
LFLAGS    =
DEFINES   =

INCLUDES =
INCLUDES += -I ../cppunit/include

LIBS =
LIBS += -L../cppunit/lib -lcppunit