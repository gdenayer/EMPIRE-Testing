CC  = mpiicc
CXX = mpiicpc
FC  = ifort
LINKER = $(CXX)

CFLAGS   = -g -xHost -O3 -std=c99 -vec-report3
CXXFLAGS = $(CFLAGS)
FCFLAGS  = 
LFLAGS   =  -fopenmp   							# this is necessary for mkl library
DEFINES   =
INCLUDES  = -I $(EMPIRE_MAPPER_LIB_INC_ON_MACHINE)			# MapperLib.h can be provided manually
INCLUDES += -I /opt/intel/composer_xe_2013_sp1.1.106/mkl/include
LIBS      = $(EMPIRE_MAPPER_LIB_ON_MACHINE) 				# libEMPIRE_MapperLib.a can also be provided manually
LIBS     += -Wl,--start-group /opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64/libmkl_intel_thread.a /opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread

# $(EMPIRE_MAPPER_LIB_ON_MACHINE) will only work after startEMPIRE
# second line of LIBS can be different depending on system