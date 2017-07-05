#pragma once
// STL
#include <string>
// #include <utility>
// ROOT and RooFit
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"

class RooCategory;
class RooWorkspace;
class TFile;

namespace SpuriousSignal {
  class SignalModel {
  public:
    /**
     * SignalModel constructor
     */
    SignalModel(const std::string& mass_category, const std::string& tag_category);

    /**
     * Build the simultaneous PDF for this category
     */
    void build_simultaneous_PDF(RooRealVar& mass);

    /**
     * Fit to dataset
     */
    void fit(RooDataSet& data);

    /**
     * Plot fit results
     */
    void plot();

    /**
     * Write workspace to file
     */
    void write(const std::string& output_file_name);

    /**
     * Expose mass points publicly, since these are needed externally to construct the combined dataset
     */
    RooCategory* mass_points();

  private:
    std::string m_mass_category;
    std::string m_tag_category;
    RooWorkspace* m_wk;
    RooDataSet* m_data;

    void add_mass_point(const int& resonance_mass);

    // std::pair<double, double> two_sigma_window(const std::string& mass_category, const std::string& tag_category);
  };
}
