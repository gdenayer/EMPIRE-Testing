CXX = icpc
CPPFLAGS = -O3 -c -DANN -DINTEL_MKL

# Compile objects
SOURCES =
SOURCES += src/AbstractMortarMapper.cpp
SOURCES += src/QuadMortarMapper.cpp 
SOURCES += src/TriangleMortarMapper.cpp 
SOURCES += src/math/MortarMath.cpp

OBJECTS = $(SOURCES:.cpp=.o)

%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(INCLUDES) $< -o $@

# Publish the lib and the header files
install : lib/libMortarMapper.a

lib/libMortarMapper.a: $(OBJECTS)
	rm -rf lib
	rm -rf include
	mkdir lib
	mkdir include
	ar -cvq $@ $^
	cp src/*.h       include

# Clean
clean:
	rm -rf $(OBJECTS) lib include

# All Includes
INCLUDES =
INCLUDES += -I src/math
INCLUDES += -I $(ANN_HOME)/include

# Libraries' pathes
ANN_HOME = ../ann_1.1.2

# File dependency
src/AbstractMortarMapper.o : src/AbstractMortarMapper.h src/math/MortarMath.h
src/QuadMortarMapper.o : src/QuadMortarMapper.h  src/AbstractMortarMapper.h src/math/MortarMath.h
src/TriangleMortarMapper.o : src/TriangleMortarMapper.h  src/AbstractMortarMapper.h src/math/MortarMath.h
src/math/MortarMath.o : src/math/MortarMath.h
