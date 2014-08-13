# required ENV variables: 
# DELPHES

ROOTCINT      = $(ROOTSYS)/bin/rootcint
ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

#CXXFLAGS      = -g -O3 -Wall -fPIC -pipe -Wuninitialized
CXXFLAGS      = -g -O0 -Wall -fPIC -pipe -Wuninitialized -O 
LD            = $(CXX)
LDFLAGS       = -g 
SOFLAGS       = -g -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) 
GLIBS         = $(filter-out -lz, $(ROOTGLIBS)) # -lTMVA -lRooFitCore -lRooFit -lRooStats

EXTHEADERS    = -I../util
LIBPATH       = $(shell pwd)/lib

READER = anaLq.o ExRootTreeReader.o runLq.o
ANA = plotLq.o 

DICTFILES = ${ANA:.o=Dict.o}
DICTHEADERS = ${ANA:.o=Dict.h}


# -- Default rules
$(addprefix obj/,%.o) : %.cc 
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c $< -o $@

%Dict.cc : %.hh %LinkDef.h
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $^ 

%Dict.cc : %.hh
	$(ROOTCINT) -f $@ -c $(EXTHEADERS) $< 

.PHONY: prep all clean vars

# ================================================================================
all: vars prep lib bin
# -----------------------------------------------------------------------

# -- library
lib/libLq.so: $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) 
	$(CXX) $(SOFLAGS) $(GLIBS) $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) $(LIBPATH)/libDelphes.so $(LIBPATH)/libutil.so -o lib/libLq.so 

# -- binaries
bin: lib/libLq.so obj/runLq.o 
	$(LD) $(LDFLAGS) -o bin/runLq $(GLIBS) obj/runLq.o $(LIBPATH)/libLq.so $(LIBPATH)/libDelphes.so $(LIBPATH)/libutil.so

# -- other stuff
obj/ExRootTreeReader.o: 
	$(CXX) $(CXXFLAGS) $(EXTHEADERS) -c delphes/ExRootTreeReader.cc -o $@


# -- preparatory setup
prep:
	mkdir -p obj bin lib
	cd lib && ln -f -s ../../util/lib/libutil.so && cd - 
	cd delphes && ln -f -s  $(DELPHES)/classes && cd - 
	cd lib && ln -f -s $(DELPHES)/libDelphes.so && cd - 

# -- clean up
clean:
	rm -f $(addprefix obj/,$(ANA) $(READER) $(DICTFILES)) 
	rm -f $(DICTHEADERS) 
	rm -f delphes/classes
	rm -f bin/runLq
	rm -f lib/*

# -- ensure that the environment variable DELPHES is set
vars:
ifndef DELPHES
    $(error DELPHES is undefined, please set it to your local Delphes installation)
endif
ifndef CMSSW_BASE
    $(error CMSSW_BASE is undefined, please run cmsenv somewhere)
endif
