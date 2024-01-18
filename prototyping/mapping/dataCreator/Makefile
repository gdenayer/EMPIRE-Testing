CC = g++
CFLAGS = -c -Wall #-g

# Compile objects
SOURCES =
SOURCES += src/FieldCreator.cpp 
SOURCES += src/FluxCreator.cpp 

OBJECTS = $(SOURCES:.cpp=.o)

# // ".cpp.o:" does not work because then the files in block file dependency
# // must be ordered alphabetically  
%.o : %.cpp
	$(CC) $(CFLAGS) $< -o $@

# Publish the lib and the header files
install : lib/libDataCreator.a

lib/libDataCreator.a : $(OBJECTS)
	ls  $(OBJECTS)
	rm -rf lib
	rm -rf include
	mkdir lib
	mkdir include
	ar -cvq $@ $^
	cp src/FieldCreator.h   include
	cp src/FluxCreator.h   include
	cp src/AbstractDataCreator.h   include

# Clean
clean:
	rm -rf $(OBJECTS) lib include
	
# File dependency
src/FieldCreator.o : src/FieldCreator.h src/AbstractDataCreator.h
src/FluxCreator.o: src/FluxCreator.h src/AbstractDataCreator.h
