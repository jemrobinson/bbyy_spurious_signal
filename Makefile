HEADER_DIR   := include
#EXTERNAL_DIR := external
BIN_DIR      := bin
SOURCE_DIR   := src
OBJECT_DIR   := obj
DEBUG_DIR    := debug
MKDIR        := mkdir -p

SOURCES := $(shell find $(SOURCE_DIR) -name "[^.]*.cxx")
OBJECTS := $(patsubst $(SOURCE_DIR)/%.cxx, $(OBJECT_DIR)/%.o, $(SOURCES))
DEBUGOBJECTS := $(patsubst $(SOURCE_DIR)/%.cxx, $(DEBUG_DIR)/%.o, $(SOURCES))
#EXTERNAL_PACKAGES := $(shell ls -d $(EXTERNAL_DIR)/*/)

CXX := g++ -m64 -Wall -Wextra
OPTIMIZE := -g -O2

DEBUG := -g -DDEBUG
CXXFLAGSDEBUG := -I/afs/cern.ch/sw/lcg/external/tcmalloc/1.7/x86_64-slc5-gcc43-opt/include/
LIBSDEBUG := -L/afs/cern.ch/sw/lcg/external/tcmalloc/1.7/x86_64-slc5-gcc43-opt/lib/ -lprofiler -ltcmalloc

CXXFLAGS := -I$(HEADER_DIR)
#CXXFLAGS += $(patsubst %,-I%,$(EXTERNAL_PACKAGES))
#CXXFLAGS += -I$(BOOST_INC_DIR)
CXXFLAGS += $(shell root-config --cflags)

LIBS := $(shell root-config --libs) -lRooFit 
#-lmapDict -lTreePlayer -lMinuit
#LIBS += -L$(EXTERNAL_DIR)/GoodRunsLists-00-01-03/StandAlone -lGoodRunsLists
#LIBS += -L$(EXTERNAL_DIR)/JetUncertainties-00-03-03/StandAlone -lJetUncertainties
#LIBS += -L$(EXTERNAL_DIR)/JetResolution2010-00-00-03/StandAlone/ -lJetResolution2010
#LIBS += -L$(EXTERNAL_DIR)/TriggerD3PDHelpers/ -lTriggerD3PDHelpers
LIBS += -L$(BOOST_LIB_DIR) -lboost_regex

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cxx
	${MKDIR} ${OBJECT_DIR} >2 /dev/null
	$(CXX) $(OPTIMIZE) -c -o $@ $< $(CXXFLAGS)

# $(DEBUG_DIR)/%.o: $(SOURCE_DIR)/%.cxx
# 	$(CXX) $(DEBUG) -c -o $@ $< $(CXXFLAGS) $(CXXFLAGSDEBUG)

.PHONY: all clean
	
all: bin/fitFourBodyMass

$(BIN_DIR)/fitFourBodyMass: $(OBJECTS)
	${MKDIR} ${BIN_DIR} >2 /dev/null
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

#debug: $(DEBUGOBJECTS)
#	$(CXX) $(DEBUG) -o bin/debug $^ $(LIBS) $(LIBSDEBUG) $(CXXFLAGS) $(CXXFLAGSDEBUG)

clean:
	rm -rf $(OBJECT_DIR)/*.o $(DEBUG_DIR)/*.o ${BIN_DIR}/* out.prof
