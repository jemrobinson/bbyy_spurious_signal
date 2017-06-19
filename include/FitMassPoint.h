#pragma once
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include <string>
#include "TFile.h"

namespace SpuriousSignal {
  class FitMassPoint {
  public:
    /**
     * FitMassPoint constructor
     */
    FitMassPoint(RooDataSet& data, std::vector<RooAbsPdf*> fit_functions, const std::string& mass_category);
    
    int fit(const bool& verbose=false);
    void plot(RooPlot* frame, const int& resonance_mass, TFile& output_file);
    
  private:
    RooDataSet m_data;
    std::vector<RooAbsPdf*> m_fit_functions;
    std::string m_mass_category;
    std::vector<int> m_colours;
  };
}
