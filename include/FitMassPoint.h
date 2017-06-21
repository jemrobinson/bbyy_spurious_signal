#pragma once
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include <string>

namespace SpuriousSignal {
  class FitMassPoint {
  public:
    /**
     * FitMassPoint constructor
     */
    FitMassPoint(RooDataSet& data, std::vector<RooAbsPdf*> fit_functions, const std::string& mass_category, const std::string& tag_category, const bool& verbose = false);

    void fit();
    void plot(RooPlot* frame, const int& resonance_mass);
    void write(const std::string& f_output_ROOT, const std::string& f_output_text) const;

  private:
    std::string m_mass_category;
    std::string m_tag_category;
    RooDataSet m_data;
    std::vector<RooAbsPdf*> m_fit_functions;
    std::map<std::string, TGraph*> m_fit_graphs;
    std::vector<double> m_nSig;
    std::vector<double> m_nSigError;
    std::vector<double> m_nSigError_withSumW2;
    std::vector<double> m_chi2;
    std::vector<int> m_ndof;
    int m_resonance_mass;
    std::vector<int> m_colours;
    std::map<std::string, std::string> m_fn_names;
    bool m_verbose;

    bool bkg_only() const;
    std::string bkg_name(RooAbsPdf* fit_function) const;
  };
}
