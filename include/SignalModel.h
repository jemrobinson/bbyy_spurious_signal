#pragma once
#include <string>

namespace SpuriousSignal {
  class SignalModel {
  public:
    /**
     * SignalModel constructor
     */
    SignalModel(const std::string& input_file, const std::string& workspace_name);
    // void build(RooRealVar& mass_resonance, RooRealVar& mass);
    // RooWorkspace workspace();

    double m_a_muGANom_SM_c2;
    double m_b_muGANom_SM_c2;
    double m_c_muGANom_SM_c2;
    double m_a_sigmaGANom_SM_c2;
    double m_b_sigmaGANom_SM_c2;
    // Crystal Ball
    double m_a_muCBNom_SM_c2;
    double m_b_muCBNom_SM_c2;
    double m_c_muCBNom_SM_c2;
    double m_a_sigmaCBNom_SM_c2;
    double m_b_sigmaCBNom_SM_c2;
    double m_nCB_SM_c2;
    double m_a_alphaCB_SM_c2;
    double m_b_alphaCB_SM_c2;
    // Combination
    double m_fracCB_SM_c2;
  };
}
