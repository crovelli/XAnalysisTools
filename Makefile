ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./
CXX	         += -I$(INCLUDEDIR) -I.
OUTLIB	         = $(INCLUDEDIR)/lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/


$(OUTLIB)XBase.o: $(INCLUDEDIR)/src/XBase.C \
	$(INCLUDEDIR)/src/X.cc \
	$(INCLUDEDIR)/src/BasicPlots.cc \
	$(INCLUDEDIR)/src/Selection.cc \
	$(INCLUDEDIR)/src/Trigger.cc \
	$(INCLUDEDIR)/src/GenLevelStudy.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)XBase.o $<
$(OUTLIB)X.o: $(INCLUDEDIR)/src/X.cc 
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)X.o $<
$(OUTLIB)BasicPlots.o: $(INCLUDEDIR)/src/BasicPlots.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)BasicPlots.o $<
$(OUTLIB)Selection.o: $(INCLUDEDIR)/src/Selection.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Selection.o $<
$(OUTLIB)Trigger.o: $(INCLUDEDIR)/src/Trigger.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)Trigger.o $<
$(OUTLIB)GenLevelStudy.o: $(INCLUDEDIR)/src/GenLevelStudy.cc
	$(CXX) $(CXXFLAGS) -c -I$(INCLUDEDIR) -o $(OUTLIB)GenLevelStudy.o $<

# ==================== XApp =============================================
XApp:  $(INCLUDEDIR)/src/XApp.C \
	$(OUTLIB)XBase.o \
	$(OUTLIB)X.o 
	$(CXX) $(CXXFLAGS) -ldl -o XApp $(OUTLIB)/*.o $(GLIBS) $(LDFLAGS) $ $<
XApp.clean:
	rm -f XApp

# ==================== reduced trees =============================================

clean:
	rm -f $(OUTLIB)*.o
	rm -f XApp

all:  XApp
