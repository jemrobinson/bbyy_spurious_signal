#pragma once
// STL
#include <string>
// ROOT and RooFit
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include <map>
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
    // void plot(std::map<std::string, RooDataSet*> dataset_map);

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

    // void set_initial_values(const int& resonance_mass);
  };
}
