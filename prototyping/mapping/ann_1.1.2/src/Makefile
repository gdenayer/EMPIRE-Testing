#-----------------------------------------------------------------------------
# Variable definition
#-----------------------------------------------------------------------------

BASEDIR	= ..
INCDIR	= $(BASEDIR)/include
LIBDIR	= $(BASEDIR)/lib

SOURCES = ANN.cpp brute.cpp kd_tree.cpp kd_util.cpp kd_split.cpp \
	kd_dump.cpp kd_search.cpp kd_pr_search.cpp kd_fix_rad_search.cpp \
	bd_tree.cpp bd_search.cpp bd_pr_search.cpp bd_fix_rad_search.cpp \
	perf.cpp

HEADERS = kd_tree.h kd_split.h kd_util.h kd_search.h \
	kd_pr_search.h kd_fix_rad_search.h perf.h pr_queue.h pr_queue_k.h

OBJECTS = $(SOURCES:.cpp=.o)


ANNLIB = libANN.a 
C++ = g++ 
CFLAGS = -O3 
MAKELIB = ar ruv 

#-----------------------------------------------------------------------------
# Make the library
#-----------------------------------------------------------------------------

install: $(LIBDIR)/$(ANNLIB)

$(LIBDIR)/$(ANNLIB): $(OBJECTS)
	rm -rf $(LIBDIR) 
	mkdir $(LIBDIR)
	$(MAKELIB) $(ANNLIB) $(OBJECTS)
	mv $(ANNLIB) $(LIBDIR)

#-----------------------------------------------------------------------------
# Make object files
#-----------------------------------------------------------------------------

ANN.o: ANN.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) ANN.cpp

brute.o: brute.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) brute.cpp

kd_tree.o: kd_tree.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_tree.cpp

kd_util.o: kd_util.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_util.cpp

kd_split.o: kd_split.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_split.cpp

kd_search.o: kd_search.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_search.cpp

kd_pr_search.o: kd_pr_search.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_pr_search.cpp

kd_fix_rad_search.o: kd_fix_rad_search.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_fix_rad_search.cpp

kd_dump.o: kd_dump.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) kd_dump.cpp

bd_tree.o: bd_tree.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) bd_tree.cpp

bd_search.o: bd_search.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) bd_search.cpp

bd_pr_search.o: bd_pr_search.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) bd_pr_search.cpp

bd_fix_rad_search.o: bd_fix_rad_search.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) bd_fix_rad_search.cpp

perf.o: perf.cpp
	$(C++) -c -I$(INCDIR) $(CFLAGS) perf.cpp

#-----------------------------------------------------------------------------
# Cleaning
#-----------------------------------------------------------------------------

clean:
	-rm -f *.o core

