// #include "RooAddPdf.h"
// #include "RooCBShape.h"
// #include "RooConstVar.h"
// #include "RooFormulaVar.h"
// #include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "SignalModel.h"
#include "TFile.h"


namespace SpuriousSignal {
  /**
   * SignalModel constructor
   */
  SignalModel::SignalModel(const std::string& input_file, const std::string& workspace_name) {
    TFile f_inputWS(input_file.c_str());
    RooWorkspace* wk = dynamic_cast<RooWorkspace*>(f_inputWS.Get(workspace_name.c_str()));

    // Gaussian
    m_a_muGANom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("a_muGANom_SM_c2"))->getVal();
    m_b_muGANom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("b_muGANom_SM_c2"))->getVal();
    m_c_muGANom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("c_muGANom_SM_c2"))->getVal();
    m_a_sigmaGANom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("a_sigmaGANom_SM_c2"))->getVal();
    m_b_sigmaGANom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("b_sigmaGANom_SM_c2"))->getVal();
    // Crystal Ball
    m_a_muCBNom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("a_muCBNom_SM_c2"))->getVal();
    m_b_muCBNom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("b_muCBNom_SM_c2"))->getVal();
    m_c_muCBNom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("c_muCBNom_SM_c2"))->getVal();
    m_a_sigmaCBNom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("a_sigmaCBNom_SM_c2"))->getVal();
    m_b_sigmaCBNom_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("b_sigmaCBNom_SM_c2"))->getVal();
    m_nCB_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("nCB_SM_c2"))->getVal();
    m_a_alphaCB_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("a_alphaCB_SM_c2"))->getVal();
    m_b_alphaCB_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("b_alphaCB_SM_c2"))->getVal();
    // Combination
    m_fracCB_SM_c2 = dynamic_cast<RooRealVar*>(wk->obj("fracCB_SM_c2"))->getVal();
    f_inputWS.Close();
  }

  // void SignalModel::build(RooRealVar& mass_resonance, RooRealVar& mass) {
  //   RooConstVar a_muGANom_SM_c2("a_muGANom_SM_c2", "a_muGANom_SM_c2", m_a_muGANom_SM_c2);
  //   RooConstVar b_muGANom_SM_c2("b_muGANom_SM_c2", "b_muGANom_SM_c2", m_b_muGANom_SM_c2);
  //   RooConstVar c_muGANom_SM_c2("c_muGANom_SM_c2", "c_muGANom_SM_c2", m_c_muGANom_SM_c2);
  //   RooConstVar a_sigmaGANom_SM_c2("a_sigmaGANom_SM_c2", "a_sigmaGANom_SM_c2", m_a_sigmaGANom_SM_c2);
  //   RooConstVar b_sigmaGANom_SM_c2("b_sigmaGANom_SM_c2", "b_sigmaGANom_SM_c2", m_b_sigmaGANom_SM_c2);
  //   // Crystal Ball
  //   RooConstVar a_muCBNom_SM_c2("a_muCBNom_SM_c2", "a_muCBNom_SM_c2", m_a_muCBNom_SM_c2);
  //   RooConstVar b_muCBNom_SM_c2("b_muCBNom_SM_c2", "b_muCBNom_SM_c2", m_b_muCBNom_SM_c2);
  //   RooConstVar c_muCBNom_SM_c2("c_muCBNom_SM_c2", "c_muCBNom_SM_c2", m_c_muCBNom_SM_c2);
  //   RooConstVar a_sigmaCBNom_SM_c2("a_sigmaCBNom_SM_c2", "a_sigmaCBNom_SM_c2", m_a_sigmaCBNom_SM_c2);
  //   RooConstVar b_sigmaCBNom_SM_c2("b_sigmaCBNom_SM_c2", "b_sigmaCBNom_SM_c2", m_b_sigmaCBNom_SM_c2);
  //   RooConstVar nCB_SM_c2("nCB_SM_c2", "nCB_SM_c2", m_nCB_SM_c2);
  //   RooConstVar a_alphaCB_SM_c2("a_alphaCB_SM_c2", "a_alphaCB_SM_c2", m_a_alphaCB_SM_c2);
  //   RooConstVar b_alphaCB_SM_c2("b_alphaCB_SM_c2", "b_alphaCB_SM_c2", m_b_alphaCB_SM_c2);
  //   // Combination
  //   RooConstVar fracCB_SM_c2("fracCB_SM_c2", "fracCB_SM_c2", m_fracCB_SM_c2);

  //   // Construct PDF
  //   RooFormulaVar signal_gaus_mean("signal_gaus_mean", "a_muGANom_SM_c2 + b_muGANom_SM_c2 * (mass_resonance / 100. - 1) + c_muGANom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muGANom_SM_c2, b_muGANom_SM_c2, c_muGANom_SM_c2, mass_resonance));
  //   RooFormulaVar signal_gaus_sigma("signal_gaus_sigma", "a_sigmaGANom_SM_c2 + b_sigmaGANom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaGANom_SM_c2, b_sigmaGANom_SM_c2, mass_resonance));
  //   RooGaussian signal_gaus("signal_gaus", "signal gaus", mass, signal_gaus_mean, signal_gaus_sigma);
  //   RooFormulaVar signal_CB_mean("signal_CB_mean", "a_muCBNom_SM_c2 + b_muCBNom_SM_c2 * (mass_resonance / 100. - 1) + c_muCBNom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muCBNom_SM_c2, b_muCBNom_SM_c2, c_muCBNom_SM_c2, mass_resonance));
  //   RooFormulaVar signal_CB_sigma("signal_CB_sigma", "a_sigmaCBNom_SM_c2 + b_sigmaCBNom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaCBNom_SM_c2, b_sigmaCBNom_SM_c2, mass_resonance));
  //   RooFormulaVar signal_CB_alpha("signal_CB_alpha", "a_alphaCB_SM_c2 + b_alphaCB_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_alphaCB_SM_c2, b_alphaCB_SM_c2, mass_resonance));
  //   RooFormulaVar signal_CB_n("signal_CB_n", "nCB_SM_c2", RooArgList(nCB_SM_c2));
  //   RooCBShape signal_CB("signal_CB", "signal_CB", mass, signal_CB_mean, signal_CB_sigma, signal_CB_alpha, signal_CB_n);
  //   RooFormulaVar signal_frac_CB("signal_frac_CB", "fracCB_SM_c2", RooArgList(fracCB_SM_c2));
  //   RooAddPdf model("signal", "CB + Gaussian", RooArgList(signal_CB, signal_gaus), signal_frac_CB);
  //   m_wk.import(model);
  // }

  // RooWorkspace SignalModel::workspace() {
  //   return m_wk;
  // }
}
