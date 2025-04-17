CC  = icx
CXX = icpx
FC  = ifx
LINKER = $(CXX)
CPPFLAGS = -O3

# Libraries to be linked (the order MATTERS!)
LIBS =
LIBS += -L$(DATACREATOR_HOME)/lib -lDataCreator
LIBS += -L$(TINYXML_HOME)/lib -lTinyxml
LIBS += -L$(MORTAR_HOME)/lib -lMortarMapper
LIBS += -L$(GID_HOME)/lib -lGiDParser
LIBS += -L$(ANN_HOME)/lib -lANN

# using LAPACK
LIBS += -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -qopenmp -lpthread

# All Includes
# difference between := and = is that := expands variable instantly,
# which is wrong if SOLVER is defined later than INCLUDES
INCLUDES =
INCLUDES += -I $(ANN_HOME)/include
INCLUDES += -I $(GID_HOME)/include
INCLUDES += -I $(TINYXML_HOME)/include
INCLUDES += -I $(MORTAR_HOME)/include
INCLUDES += -I $(DATACREATOR_HOME)/include

# Libraries' pathes
ANN_HOME = ../ann_1.1.2
GID_HOME = ../GiDParser
TINYXML_HOME = ../tinyxml
MORTAR_HOME = ../mortarMapper
DATACREATOR_HOME = ../dataCreator
