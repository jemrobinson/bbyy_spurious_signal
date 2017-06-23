BIN_DIR      := bin
DEBUG_DIR    := debug
HEADER_DIR   := include
INPUT_DIR    := input
OBJECT_DIR   := obj
OUTPUT_DIR   := output
SOURCE_DIR   := src

SOURCES := $(shell find $(SOURCE_DIR) -name "[^.]*.cxx")
OBJECTS := $(patsubst $(SOURCE_DIR)/%.cxx, $(OBJECT_DIR)/%.o, $(SOURCES))
DEBUGOBJECTS := $(patsubst $(SOURCE_DIR)/%.cxx, $(DEBUG_DIR)/%.o, $(SOURCES))

CXX   := g++ -m64 -Wall -Wextra 
MKDIR := mkdir -p

OPTIMIZE := -O2 -DMSG_LEVEL=3
DEBUG    := -g -DMSG_LEVEL=4

CXXFLAGS := -I$(HEADER_DIR) $(shell root-config --cflags)
CXXFLAGSDEBUG := -I/afs/cern.ch/sw/lcg/external/tcmalloc/1.7/x86_64-slc5-gcc43-opt/include/

LIBS := $(shell root-config --libs) -lRooFit
LIBSDEBUG := -L/afs/cern.ch/sw/lcg/external/tcmalloc/1.7/x86_64-slc5-gcc43-opt/lib/ -ltcmalloc

.PHONY: all clean

all: $(BIN_DIR)/fitFourBodyMass

$(BIN_DIR)/fitFourBodyMass: $(OBJECTS)
	${MKDIR} ${BIN_DIR} ${INPUT_DIR} ${OUTPUT_DIR}/plots/lowMass_0tag ${OUTPUT_DIR}/plots/lowMass_1tag ${OUTPUT_DIR}/plots/lowMass_2tag ${OUTPUT_DIR}/plots/highMass_0tag ${OUTPUT_DIR}/plots/highMass_1tag ${OUTPUT_DIR}/plots/highMass_2tag
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

debug: $(DEBUGOBJECTS)
	${MKDIR} ${BIN_DIR} ${DEBUG_DIR}
	$(CXX) $(DEBUG) -o bin/debug $^ $(LIBS) $(LIBSDEBUG) $(CXXFLAGS) $(CXXFLAGSDEBUG)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cxx
	${MKDIR} ${OBJECT_DIR}
	$(CXX) $(OPTIMIZE) -c -o $@ $< $(CXXFLAGS)

$(DEBUG_DIR)/%.o: $(SOURCE_DIR)/%.cxx
	${MKDIR} ${DEBUG_DIR}
	$(CXX) $(DEBUG) -c -o $@ $< $(CXXFLAGS) $(CXXFLAGSDEBUG)

clean:
	rm -rf $(OBJECT_DIR)/*.o $(DEBUG_DIR)/*.o ${BIN_DIR}/* out.prof
