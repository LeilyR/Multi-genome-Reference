

CXX = g++
CFLAGS = -Wall -fmessage-length=50 -fopenmp 
CXXFLAGS = -Wall -fmessage-length=50 -fopenmp -pedantic 
#OPTFLAGS = -O3 -DNDEBUG
OPTFLAGS = -O3 -DNDEBUG 
DEBUGFLAGS = -g3 -O0 -gdwarf-2 

DEBUG = 1
#OPT = 1

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

CXXFLAGS += -Wlong-long


STATGEN_INC = libStatGen-1.0.13/include
INCLUDEARGS = -I$(STATGEN_INC)
LIBARGS =   apclusterunix64.so -L../dlib_build -ldlib  libStatGen-1.0.13/libStatGen.a -lz


# include libStatGen-1.0.13/Makefiles/Makefile.include

PROGS = graph
OBJS = pw_alignment.o data.o encoder.o test.o

#model.o


# .DEFAULT_GOAL := graph 
all: $(PROGS)




graph: main.cpp $(OBJS) model.cpp model.hpp pw_alignment.hpp data.hpp graph.hpp needleman.hpp
	$(CXX) -o graph $(CXXFLAGS) $(INCLUDEARGS) main.cpp $(OBJS) $(LIBARGS)
	
graph.o: graph.cpp data.hpp needleman.hpp graph.hpp
	$(CXX) -o graph.o $(CXXFLAGS) $(INCLUDEARGS) -c graph.cpp

pw_alignment.o: pw_alignment.cpp pw_alignment.hpp
	$(CXX) -o pw_alignment.o $(CXXFLAGS) $(INCLUDEARGS) -c pw_alignment.cpp 

data.o: data.cpp pw_alignment.hpp data.hpp
	$(CXX) -o data.o $(CXXFLAGS) $(INCLUDEARGS) -c data.cpp 

encoder.o: encoder.cpp encoder.hpp data.hpp
	$(CXX) -o encoder.o $(CXXFLAGS) $(INCLUDEARGS) -c encoder.cpp

test.o: test.cpp test.hpp 

needleman.o: needleman.cpp needleman.hpp
	$(CXX) -o needleman.o $(CXXFLAGS) $(INCLUDEARGS) -c needleman.cpp

#model.o: model.cpp pw_alignment.hpp data.hpp model.hpp
#model.o: model.cpp pw_alignment.hpp data.hpp model.hpp
#	$(CXX) -o model.o $(CXXFLAGS) $(INCLUDEARGS) -c model.cpp

clean:
	@-rm -f $(PROGS) *.o  *.a *~





