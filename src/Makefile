

CXX = g++
CFLAGS = -Wall -fmessage-length=50 -fopenmp 
CXXFLAGS = -Wall -fmessage-length=50 -fopenmp -pedantic 
#OPTFLAGS = -O3 -DNDEBUG
OPTFLAGS = -O3 -DNDEBUG 
DEBUGFLAGS = -g -O0  

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
LIBARGS = 



PROGS = graph
OBJS = pw_alignment.o


all: $(PROGS)

graph: main.cpp $(OBJS) 
	$(CXX) -o graph $(CXXFLAGS) $(INCLUDEARGS) main.cpp $(OBJS) $(LIBARGS)
	
pw_alignment.o: pw_alignment.cpp
	$(CXX) -o pw_alignment.o $(CXXFLAGS) $(INCLUDEARGS) -c pw_alignment.cpp 

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~





