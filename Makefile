BIN_DIR      := bin
HEADER_DIR   := include
INPUT_DIR    := input
OBJECT_DIR   := obj
OUTPUT_DIR   := output
PLOT_DIR     := plots
SOURCE_DIR   := src

SOURCES := $(shell find $(SOURCE_DIR) -name "[^.]*.cxx")
FOUR_BODY_MASS_COMPONENTS := fitFourBodyMass ParameterSet PDFModelFitter PlotStyle SignalModel
SIGNAL_SHAPE_COMPONENTS := fitSignalShape PlotStyle ExpGausExpPDF ExpGausExpPDFDict SignalModel
SIGNAL_BIAS_COMPONENTS := extractSignalBias PlotStyle SignalModel ToyGenerator
CXX   := g++ -m64 -Wall -Wextra# -Werror
MKDIR := mkdir -p

OPTIMIZE := -g -O2 -DMSG_LEVEL=3
CXXFLAGS := -std=c++11 -I$(HEADER_DIR) $(shell root-config --cflags)
LIBS := $(shell root-config --libs) -lRooFit

.PHONY: all clean directories

all: directories $(BIN_DIR)/extractSignalBias $(BIN_DIR)/fitFourBodyMass  $(BIN_DIR)/fitSignalShape

clean:
	rm -rf $(OBJECT_DIR)/*.o $(BIN_DIR)/* $(SOURCE_DIR)/*PDFDict*

directories:
	${MKDIR} ${BIN_DIR} ${INPUT_DIR} ${OBJECT_DIR} ${OUTPUT_DIR}

$(BIN_DIR)/extractSignalBias: $(patsubst %, $(OBJECT_DIR)/%.o, $(SIGNAL_BIAS_COMPONENTS))
	${MKDIR} ${PLOT_DIR}/signal_bias/asimov_fits ${PLOT_DIR}/signal_bias/toy_fits
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

$(BIN_DIR)/fitFourBodyMass: $(patsubst %, $(OBJECT_DIR)/%.o, $(FOUR_BODY_MASS_COMPONENTS))
	${MKDIR} ${OUTPUT_DIR}/csv/mass_points ${OUTPUT_DIR}/csv/bkg_only
	${MKDIR} ${PLOT_DIR}/background_model
	${MKDIR} ${PLOT_DIR}/signal_plus_background_fits/lowMass_0tag ${PLOT_DIR}/signal_plus_background_fits/highMass_0tag
	${MKDIR} ${PLOT_DIR}/signal_plus_background_fits/lowMass_1tag ${PLOT_DIR}/signal_plus_background_fits/highMass_1tag
	${MKDIR} ${PLOT_DIR}/signal_plus_background_fits/lowMass_2tag ${PLOT_DIR}/signal_plus_background_fits/highMass_2tag
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

$(BIN_DIR)/fitSignalShape: $(patsubst %, $(OBJECT_DIR)/%.o, $(SIGNAL_SHAPE_COMPONENTS))
	${MKDIR} ${PLOT_DIR}/signal_model/overall
	${MKDIR} ${PLOT_DIR}/signal_model/EGE
	${MKDIR} ${PLOT_DIR}/signal_model/CBGA
	$(CXX) $(OPTIMIZE) -o $@ $^ $(LIBS) $(CXXFLAGS)

$(OBJECT_DIR)/%.o: $(SOURCE_DIR)/%.cxx
	$(CXX) $(OPTIMIZE) -c -o $@ $< $(CXXFLAGS)

$(SOURCE_DIR)/DoubleSidedCrystalBallPDFDict.cxx:
	rootcling -f src/DoubleSidedCrystalBallPDFDict.cxx -s bin/DoubleSidedCrystalBallPDFDict -I./include DoubleSidedCrystalBallPDF.h

$(SOURCE_DIR)/ExpGausExpPDFDict.cxx:
	rootcling -f src/ExpGausExpPDFDict.cxx -s bin/ExpGausExpPDFDict -I./include ExpGausExpPDF.h

