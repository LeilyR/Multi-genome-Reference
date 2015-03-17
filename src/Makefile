

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

INCLUDEDIRS =
INCLUDEARGS =
LIBARGS =  apclusterunix64.so -L../dlib_build -ldlib



PROGS = graph
OBJS = pw_alignment.o data.o encoder.o 

#model.o


all: $(PROGS)

graph: main.cpp $(OBJS) model.cpp model.hpp pw_alignment.hpp data.hpp
	$(CXX) -o graph $(CXXFLAGS) $(INCLUDEARGS) main.cpp $(OBJS) $(LIBARGS)
	
pw_alignment.o: pw_alignment.cpp pw_alignment.hpp
	$(CXX) -o pw_alignment.o $(CXXFLAGS) $(INCLUDEARGS) -c pw_alignment.cpp 

data.o: data.cpp pw_alignment.hpp data.hpp
	$(CXX) -o data.o $(CXXFLAGS) $(INCLUDEARGS) -c data.cpp 
encoder.o: encoder.cpp encoder.hpp data.hpp
	$(CXX) -o encoder.o $(CXXFLAGS) $(INCLUDEARGS) -c encoder.cpp

#model.o: model.cpp pw_alignment.hpp data.hpp model.hpp
#model.o: model.cpp pw_alignment.hpp data.hpp model.hpp
#	$(CXX) -o model.o $(CXXFLAGS) $(INCLUDEARGS) -c model.cpp

clean:
	@-rm -f $(PROGS) *.o  *.a *~





