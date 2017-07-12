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
    RooCategory* mass_points() const;

  private:
    std::string m_mass_category;
    std::string m_tag_category;
    RooWorkspace* m_wk;
    RooDataSet* m_data;
    std::vector<std::string> m_models;
    std::map< std::string, std::vector<std::string> > m_model_parameters;
    std::map< std::string, std::vector<std::string> > m_model_metaparameters;

    /**
     * Add a mass point to the simultaneous PDF
     */
    void add_mass_point(const int& resonance_mass);

    void add_parameterisation(const std::string& name);

    void append(std::vector<double>& target, std::vector<double> source);

    double get_initial_value(const std::string& name);

  };
}
