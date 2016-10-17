

CXX = g++ -g -std=c++0x
CFLAGS = -Wall -fmessage-length=50 -fopenmp #-march=bdver2
CXXFLAGS = -Wall -fmessage-length=50 -fopenmp -pedantic #-march=bdver2
#OPTFLAGS = -O3 -DNDEBUG
#OPTFLAGS = -O3 -DNDEBUG 
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
#OBJS = pw_alignment.o data.o test.o graph.o needleman.o dynamic_mc.o suffix_tree.o merging_centers.o  read_daligner.o needleman_wunsch.o map_read.o 
OBJS = pw_alignment.o data.o test.o graph.o needleman.o dynamic_mc.o suffix_tree.o merging_centers.o  read_daligner.o needleman_wunsch.o map_read.o 

# .DEFAULT_GOAL := graph 
all: $(PROGS)




#graph: main.cpp $(OBJS) model.cpp model.hpp pw_alignment.hpp data.hpp graph.hpp needleman.hpp encoder.hpp dynamic_encoder.hpp dynamic_decoder.hpp IntervalTree.hpp

graph: main.cpp $(OBJS) model.cpp model.hpp pw_alignment.hpp data.hpp graph.hpp needleman.hpp needleman_wunsch.hpp encoder.hpp dynamic_encoder.hpp dynamic_decoder.hpp suffix_tree.hpp alignment_index.hpp IntervalTree.hpp overlap.hpp connected_component.hpp merging_centers.hpp filter_alignments.hpp read_daligner.hpp map_read.hpp connected_component.hpp connected_component.cpp overlap.cpp model.cpp encoder.cpp intervals.cpp dynamic_encoder.cpp dynamic_decoder.cpp needleman_wunsch.cpp filter_alignments.cpp alignment_index.cpp
	$(CXX) -o graph $(CXXFLAGS) $(INCLUDEARGS) main.cpp $(OBJS) $(LIBARGS)
	
graph.o: graph.hpp graph.cpp
	$(CXX) -o graph.o $(CXXFLAGS) $(INCLUDEARGS) -c graph.cpp

pw_alignment.o: pw_alignment.cpp pw_alignment.hpp
	$(CXX) -o pw_alignment.o $(CXXFLAGS) $(INCLUDEARGS) -c pw_alignment.cpp 

data.o: data.cpp pw_alignment.hpp data.hpp
	$(CXX) -o data.o $(CXXFLAGS) $(INCLUDEARGS) -c data.cpp 

#alignment_index.o: alignment_index.cpp alignment_index.hpp intervals.hpp intervals.cpp
#	$(CXX) -o alignment_index.o $(CXXFLAGS) $(INCLUDEARGS) -c alignment_index.cpp
#overlap.o: overlap.cpp overlap.hpp
#	$(CXX) -o overlap.o $(CXXFLAGS) $(INCLUDEARGS) -c overlap.cpp


dynamic_mc.o: dynamic_mc.cpp dynamic_mc.hpp  data.hpp
	$(CXX) -o dynamic_mc.o $(CXXFLAGS) $(INCLUDEARGS) -c dynamic_mc.cpp

#encoder.o: encoder.cpp encoder.hpp data.hpp 
#	$(CXX) -o encoder.o $(CXXFLAGS) $(INCLUDEARGS) -c encoder.cpp

test.o: test.cpp test.hpp 

needleman.o: needleman.cpp needleman.hpp
	$(CXX) -o needleman.o $(CXXFLAGS) $(INCLUDEARGS) -c needleman.cpp

needleman_wunsch.o: needleman_wunsch.cpp needleman_wunsch.hpp
	$(CXX) -o needleman_wunsch.o $(CXXFLAGS) $(INCLUDEARGS) -c needleman_wunsch.cpp

map_read.o: map_read.cpp map_read.hpp
	$(CXX) -o map_read.o $(CXXFLAGS) $(INCLUDEARGS) -c map_read.cpp

suffix_tree.o: suffix_tree.hpp suffix_tree.cpp
	$(CXX) -o suffix_tree.o $(CXXFLAGS) $(INCLUDEARGS) -c suffix_tree.cpp

merging_centers.o: merging_centers.hpp merging_centers.cpp
	$(CXX) -o merging_centers.o $(CXXFLAGS) $(INCLUDEARGS) -c merging_centers.cpp

#filter_alignments.o: filter_alignments.hpp filter_alignments.cpp 
#	$(CXX) -o filter_alignments.o $(CXXFLAGS) $(INCLUDEARGS) -c filter_alignments.cpp

read_daligner.o: read_daligner.hpp read_daligner.cpp
	$(CXX) -o read_daligner.o $(CXXFLAGS) $(INCLUDEARGS) -c read_daligner.cpp



#model.o: model.cpp model.hpp
#	$(CXX) -o model.o $(CXXFLAGS) $(INCLUDEARGS) -c model.cpp

clean:
	@-rm -f $(PROGS) *.o  *.a *~





