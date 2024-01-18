CC = g++
CFLAGS = -c -Wall #-g

# Compile objects
SOURCES =
SOURCES += src/GiDBasicData.cpp
SOURCES += src/GiDResultInventory.cpp 
SOURCES += src/GiDMeshInventory.cpp 
SOURCES += src/GiDFileIO.cpp

OBJECTS = $(SOURCES:.cpp=.o)

# // ".cpp.o:" does not work because then the files in block file dependency
# // must be ordered alphabetically  
%.o : %.cpp
	$(CC) $(CFLAGS) $< -o $@

# Publish the lib and the header files
install : lib/libGiDParser.a

lib/libGiDParser.a : $(OBJECTS)
	ls  $(OBJECTS)
	rm -rf lib
	rm -rf include
	mkdir lib
	mkdir include
	ar -cvq $@ $^
	cp src/GiDFileIO.h   include

# Clean
clean:
	rm -rf $(OBJECTS) lib include
	
# File dependency
src/GiDBasicData.o : src/GiDBasicData.h
src/GiDResultInventory.o: src/GiDResultInventory.h src/GiDBasicData.h
src/GiDMeshInventory.o: src/GiDMeshInventory.h src/GiDBasicData.h
src/GiDFileIO.o: src/GiDFileIO.h src/GiDMeshInventory.h src/GiDResultInventory.h src/GiDBasicData.h