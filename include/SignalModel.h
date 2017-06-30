#pragma once
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include <string>

class RooCategory;
class RooDataSet;
class RooWorkspace;
class TFile;

namespace SpuriousSignal {
  class SignalModel {
  public:
    /**
     * SignalModel constructor
     */
    // SignalModel(RooDataSet& data, RooSimultaneous& sim_PDF, const std::string& mass_category, const std::string& tag_category);
    SignalModel(const std::string& mass_category, const std::string& tag_category);

    void build_simultaneous_PDF(RooRealVar& mass);
    void fit(RooDataSet& data);
    void plot();
    void write(const std::string& output_file_name);

    RooCategory* mass_points();

  private:
    std::string m_mass_category;
    std::string m_tag_category;
    RooWorkspace* m_wk;
    RooDataSet* m_data;

    void add_mass_point(const int& resonance_mass);
  };
}
