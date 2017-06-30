BIN_DIR      := bin
HEADER_DIR   := include
INPUT_DIR    := input
OBJECT_DIR   := obj
OUTPUT_DIR   := output
PLOT_DIR     := plots
SOURCE_DIR   := src

SOURCES := $(shell find $(SOURCE_DIR) -name "[^.]*.cxx")
FOUR_BODY_MASS_COMPONENTS := fitFourBodyMass PlotStyle FitMassPoint ParameterSet SignalModel
SIGNAL_SHAPE_COMPONENTS := fitSignalShape PlotStyle SignalModel 
CXX   := g++ -m64 -Wall -Wextra 
MKDIR := mkdir -p

OPTIMIZE := -g -O2 -DMSG_LEVEL=3
CXXFLAGS := -I$(HEADER_DIR) $(shell root-config --cflags)
LIBS := $(shell root-config --libs) -lRooFit

.PHONY: all clean

all: $(BIN_DIR)/fitFourBodyMass  $(BIN_DIR)/fitSignalShape

clean:
	rm -rf $(OBJECT_DIR)/*.o ${BIN_DIR}/*
	@echo $(patsubst %, $(OBJECT_DIR)/%.o, $(FOUR_BODY_MASS_COMPONENTS))


$(BIN_DIR)/fitFourBodyMass: $(patsubst %, $(OBJECT_DIR)/%.o, $(FOUR_BODY_MASS_COMPONENTS))
	${MKDIR} ${BIN_DIR} ${INPUT_DIR} ${OUTPUT_DIR}
	${MKDIR} ${PLOT_DIR}/lowMass_0tag ${PLOT_DIR}/lowMass_1tag ${PLOT_DIR}/lowMass_2tag ${PLOT_DIR}/highMass_0tag ${PLOT_DIR}/highMass_1tag ${PLOT_DIR}/highMass_2tag
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

$(BIN_DIR)/fitSignalShape: $(patsubst %, $(OBJECT_DIR)/%.o, $(SIGNAL_SHAPE_COMPONENTS))
	${MKDIR} ${BIN_DIR} ${INPUT_DIR}
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cxx
	${MKDIR} ${OBJECT_DIR}
	$(CXX) $(OPTIMIZE) -c -o $@ $< $(CXXFLAGS)
