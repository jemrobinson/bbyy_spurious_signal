// STL
#include <fstream>
#include <string>
#include <utility>
#include <vector>
// ROOT
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TLegend.h"
#include "TMath.h"
// RooFit
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooCFunction3Binding.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooNovosibirsk.h"
#include "RooPlot.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
// Local
#include "ParameterSet.h"
#include "Logger.h"

int main(int /*argc*/, char** /*argv*/) {  
  using namespace SpuriousSignal;
  
  int MASS_STEP = 20;
  int NTRIALS = 1;
    
  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL); //WARNING FATAL
  gErrorIgnoreLevel = kBreak;
  
  // Setup colours and mass categories
  std::vector<int> colours = {kViolet, kGreen + 3, kBlue, kRed};
  std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)},
                                                                  {"high", std::make_pair<int, int>(335, 1140)}};
  // std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)}};

  TFile f_inputWS("/afs/cern.ch/user/a/andari/public/Leo/SigParam_low_cbga_cut5/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root");
  RooWorkspace wk(*(RooWorkspace*)f_inputWS.Get("signalWS"));
  // Gaussian
  RooRealVar a_muGANom_SM_c2(*(RooRealVar*)wk.obj("a_muGANom_SM_c2")); a_muGANom_SM_c2.setConstant(true);
  RooRealVar b_muGANom_SM_c2(*(RooRealVar*)wk.obj("b_muGANom_SM_c2")); b_muGANom_SM_c2.setConstant(true);
  RooRealVar c_muGANom_SM_c2(*(RooRealVar*)wk.obj("c_muGANom_SM_c2")); c_muGANom_SM_c2.setConstant(true);
  RooRealVar a_sigmaGANom_SM_c2(*(RooRealVar*)wk.obj("a_sigmaGANom_SM_c2")); a_sigmaGANom_SM_c2.setConstant(true);
  RooRealVar b_sigmaGANom_SM_c2(*(RooRealVar*)wk.obj("b_sigmaGANom_SM_c2")); b_sigmaGANom_SM_c2.setConstant(true);
  // Crystal Ball
  RooRealVar a_muCBNom_SM_c2(*(RooRealVar*)wk.obj("a_muCBNom_SM_c2")); a_muCBNom_SM_c2.setConstant(true);
  RooRealVar b_muCBNom_SM_c2(*(RooRealVar*)wk.obj("b_muCBNom_SM_c2")); b_muCBNom_SM_c2.setConstant(true);
  RooRealVar c_muCBNom_SM_c2(*(RooRealVar*)wk.obj("c_muCBNom_SM_c2")); c_muCBNom_SM_c2.setConstant(true);
  RooRealVar a_sigmaCBNom_SM_c2(*(RooRealVar*)wk.obj("a_sigmaCBNom_SM_c2")); a_sigmaCBNom_SM_c2.setConstant(true);
  RooRealVar b_sigmaCBNom_SM_c2(*(RooRealVar*)wk.obj("b_sigmaCBNom_SM_c2")); b_sigmaCBNom_SM_c2.setConstant(true);
  RooRealVar nCB_SM_c2(*(RooRealVar*)wk.obj("nCB_SM_c2")); nCB_SM_c2.setConstant(true);
  RooRealVar a_alphaCB_SM_c2(*(RooRealVar*)wk.obj("a_alphaCB_SM_c2")); a_alphaCB_SM_c2.setConstant(true);
  RooRealVar b_alphaCB_SM_c2(*(RooRealVar*)wk.obj("b_alphaCB_SM_c2")); b_alphaCB_SM_c2.setConstant(true);
  // Combination
  RooRealVar fracCB_SM_c2(*(RooRealVar*)wk.obj("fracCB_SM_c2")); fracCB_SM_c2.setConstant(true);

  // std::vector<double> 

  for (auto mass_kv : mass_categories) {
    // Define data parameters
    RooRealVar mass("mass", "m_{yyjj}", mass_kv.second.first, mass_kv.second.second);
    RooRealVar weight("weight", "event weight", 1e-10, 1e10);

    // Get the data
    std::string input_directory("/afs/cern.ch/user/j/jrobinso/HGamma/PyAnalysis/text/");
    // RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    RooDataSet* ptr_raw_data = RooDataSet::read((input_directory + "m_yyjj_SM_bkg_" + mass_kv.first + "Mass_0tag_tightIsolated.txt").c_str(), RooArgList(mass, weight));
    // RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
    RooDataSet data("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));

    // Signal: CB+Gaus
    RooRealVar mass_resonance("mass_resonance", "mass of resonance", 305);
    RooFormulaVar signal_gaus_mean("signal_gaus_mean", "a_muGANom_SM_c2 + b_muGANom_SM_c2 * (mass_resonance / 100. - 1) + c_muGANom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muGANom_SM_c2, b_muGANom_SM_c2, c_muGANom_SM_c2, mass_resonance));
    RooFormulaVar signal_gaus_sigma("signal_gaus_sigma", "a_sigmaGANom_SM_c2 + b_sigmaGANom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaGANom_SM_c2, b_sigmaGANom_SM_c2, mass_resonance));
    RooGaussian signal_gaus("signal_gaus", "signal gaus", mass, signal_gaus_mean, signal_gaus_sigma);
    RooFormulaVar signal_CB_mean("signal_CB_mean", "a_muCBNom_SM_c2 + b_muCBNom_SM_c2 * (mass_resonance / 100. - 1) + c_muCBNom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muCBNom_SM_c2, b_muCBNom_SM_c2, c_muCBNom_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_sigma("signal_CB_sigma", "a_sigmaCBNom_SM_c2 + b_sigmaCBNom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaCBNom_SM_c2, b_sigmaCBNom_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_alpha("signal_CB_alpha", "a_alphaCB_SM_c2 + b_alphaCB_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_alphaCB_SM_c2, b_alphaCB_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_n("signal_CB_n", "nCB_SM_c2", RooArgList(nCB_SM_c2));
    RooCBShape signal_CB("signal_CB", "signal_CB", mass, signal_CB_mean, signal_CB_sigma, signal_CB_alpha, signal_CB_n);
    RooFormulaVar signal_frac_CB("signal_frac_CB", "fracCB_SM_c2", RooArgList(fracCB_SM_c2));
    RooAddPdf signal_PDF("signal_PDFa", "CB + Gaussian", RooArgList(signal_CB, signal_gaus), signal_frac_CB);
    
    // Novosibirsk
    RooRealVar novosibirsk_peak("novosibirsk_peak", "peak of Novosibirsk", 260, 200, 500);
    RooRealVar novosibirsk_width("novosibirsk_width", "width of Novosibirsk", 30, 0, 100);
    RooRealVar novosibirsk_tail("novosibirsk_tail", "tail of Novosibirsk", -1, -10, 10);
    RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", mass, novosibirsk_peak, novosibirsk_width, novosibirsk_tail);
    
    // // Modified Cauchy
    // RooRealVar cauchy_mu("cauchy_mu", "mean of Cauchy", -500, -1e4, 1e4);
    // RooRealVar cauchy_t0("cauchy_t0", "t0 of Cauchy", 1000, -1e5, 1e5);
    // RooRealVar cauchy_t1("cauchy_t1", "t1 of Cauchy", 0, -10, 10);
    // RooFormulaVar cauchy_t("cauchy_t", "cauchy_t0 + cauchy_t1 * mass", RooArgList(cauchy_t0, cauchy_t1, mass));
    // RooRealVar cauchy_s0("cauchy_s0", "s0 of Cauchy", -150, -1e5, 1e5);
    // RooRealVar cauchy_s1("cauchy_s1", "s1 of Cauchy", 1, -100, 100);
    // RooFormulaVar cauchy_s("cauchy_s", "cauchy_s0 + cauchy_s1 * mass", RooArgList(cauchy_s0, cauchy_s1, mass));
    // RooGenericPdf modified_cauchy("modified_cauchy", "modified_cauchy", "TMath::CauchyDist(mass - cauchy_mu, cauchy_t, cauchy_s)", RooArgList(mass, cauchy_mu, cauchy_t, cauchy_s));
    // 
    // // Modified Gamma
    // RooRealVar gamma_Gamma("gamma_Gamma", "gamma of Gamma", 5, 0, 20);
    // RooRealVar gamma_beta0("gamma_beta0", "beta0 of Gamma", -20, -20, 20);
    // RooRealVar gamma_beta1("gamma_beta1", "beta1 of Gamma", 0.1, -1, 1);
    // RooFormulaVar gamma_beta("gamma_beta", "TMath::Abs(gamma_beta0 + gamma_beta1 * mass)", RooArgList(gamma_beta0, gamma_beta1, mass));
    // RooRealVar gamma_mu("gamma_mu", "mean of Gamma", 245, 220, 260);
    // RooGenericPdf modified_gamma("modified_gamma", "modified_gamma", "TMath::GammaDist(mass, gamma_Gamma, gamma_mu, TMath::Abs(gamma_beta0 + gamma_beta1 * mass))", RooArgList(mass, gamma_Gamma, gamma_mu, gamma_beta0, gamma_beta1));
    
    // Setup fit functions
    std::vector<RooAbsPdf*> bkg_functions;
    bkg_functions.push_back(&novosibirsk);
    // bkg_functions.push_back(&modified_gamma);
    // bkg_functions.push_back(&modified_cauchy);

    // Number of signal and background events
    RooRealVar nSig("nSig", "number of signal events", 0, 1000);
    RooRealVar nBkg("nBkg", "number of background events", 0, 1e6);

    // Construct parameter sets that need to be remembered
    std::vector<ParameterSet> parameter_sets;
    parameter_sets.push_back(ParameterSet("Novosibirsk", {&novosibirsk_peak, &novosibirsk_width, &novosibirsk_tail, &nSig, &nBkg, &mass_resonance}));

    // Setup fit functions
    std::vector<RooAbsPdf*> splusb_functions;
    for (auto bkg_function : bkg_functions) {
      RooAddPdf* splusb = new RooAddPdf("signal_plus_" + bkg_function->getTitle(), "Signal + " + bkg_function->getTitle(), RooArgList(signal_PDF, *bkg_function), RooArgList(nSig, nBkg));
      splusb_functions.push_back(splusb);
    }

    
    // Recreate output text file
    std::ofstream f_output;
    f_output.open("output/spurious_signal_" + mass_kv.first + ".csv", std::ios::trunc);
    f_output.close();
    
    // Selection criterion
    double Z_spurious(0), Z_spurious_max(-99);
    
    for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
      MSG_INFO("Using fit function " << splusb_functions.at(idx)->getTitle());
      // Iterate over mass points
      for (int iMass = mass_kv.second.first; iMass < mass_kv.second.second; iMass += MASS_STEP) {
        mass_resonance.setVal(iMass); mass_resonance.setConstant(true); int fit_status(0);
        for (int iTrial = 0; iTrial < NTRIALS; ++iTrial) {
          splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::PrintLevel(-1));
          splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migrad"), RooFit::PrintLevel(-1));
          splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "improve"), RooFit::PrintLevel(-1));
          splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "hesse"), RooFit::PrintLevel(-1));
          fit_status = splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "minos"), RooFit::PrintLevel(-1), RooFit::Save(true))->status();
          if (fit_status == 0) { break; }
          MSG_ERROR("Fit did not converge: status = " << fit_status);
        }
        // Record these parameters if Z > max(Z) so far
        Z_spurious = nSig.getVal() / nSig.getError();
        MSG_DEBUG("... resonance mass = " << mass_resonance.getVal() << " => nSig: " << nSig.getVal() << ", nBkg: " << nBkg.getVal() << " [" << fit_status << "] => Z = " << Z_spurious);
        if (Z_spurious > Z_spurious_max) {
          Z_spurious_max = Z_spurious;
          MSG_INFO("Recording new max(z) = " << Z_spurious_max << " at " << mass_resonance.getVal());
          parameter_sets.at(idx).record_values();
        }
        // Write parameters to output file
        f_output.open("output/spurious_signal_" + mass_kv.first + ".csv", std::ios::app);
        f_output << bkg_function->getTitle() << " " << mass_resonance.getVal() << " " << nSig.getVal() << " " << nSig.getError() << " " << nBkg.getVal() << " " << nBkg.getError() << " " << fit_status << std::endl;
        f_output.close();
      }
    }

    // MSG_INFO("Initial");
    // MSG_INFO("novosibirsk_peak: " << novosibirsk_peak.getVal());
    // MSG_INFO("novosibirsk_width: " << novosibirsk_width.getVal());
    // MSG_INFO("novosibirsk_tail: " << novosibirsk_tail.getVal());
    // novosibirsk_params.restore_values();
    // MSG_INFO("Restored");
    // MSG_INFO("novosibirsk_peak: " << novosibirsk_peak.getVal());
    // MSG_INFO("novosibirsk_width: " << novosibirsk_width.getVal());
    // MSG_INFO("novosibirsk_tail: " << novosibirsk_tail.getVal());
    
    // // Fit PDFs to data
    // // std::vector<RooAbsPdf*> splusb_functions;
    // std::vector< std::vector<double> > parameters;
    // for (auto bkg_function : bkg_functions) {
    //   std::cout << "   INFO: " << "Using fit function " << bkg_function->getTitle() << std::endl;
    //   double nSigMax(-1e6);
    //   // RooAbsPdf* ptrSplusB = 0;
    //   // Iterate over mass points
    //   for (int iMass = mass_kv.second.first; iMass < mass_kv.second.second; iMass += MASS_STEP) {
    //     mass_resonance.setVal(iMass); mass_resonance.setConstant(true);
    //     // nSig.setVal(10); nBkg.setVal(10);
    //     RooAddPdf SplusB(bkg_function->getTitle(), "Signal + background", RooArgList(signal_PDF, *bkg_function), RooArgList(nSig, nBkg));
    //     RooFitResult* fit_result(0);
    //     for (int iTrial = 0; iTrial < NTRIALS; ++iTrial) {
    //       SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::PrintLevel(-1));
    //       SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migrad"), RooFit::PrintLevel(-1));
    //       SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "improve"), RooFit::PrintLevel(-1));
    //       SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "hesse"), RooFit::PrintLevel(-1));
    //       fit_result = SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "minos"), RooFit::PrintLevel(-1), RooFit::Save(true));
    //       if (fit_result->status() == 0) { break; }
    //       std::cout << "  ERROR: " << "Fit did not converge: status = " << fit_result->status() << std::endl;
    //     }
    //     std::cout << "   INFO: " << "... resonance mass = " << mass_resonance.getVal() << " => nSig: " << nSig.getVal() << ", nBkg: " << nBkg.getVal() << ". Fit status = " << fit_result->status() << std::endl;
    //     if (nSig.getVal() > nSigMax) {
    //       nSigMax = nSig.getVal();
    //       // delete ptrSplusB; ptrSplusB = new RooAddPdf(SplusB);
    //     }
    //     // Write parameters to output file
    //     f_output.open("output/spurious_signal_" + mass_kv.first + ".csv", std::ios::app);
    //     f_output << bkg_function->getTitle() << " " << mass_resonance.getVal() << " " << nSig.getVal() << " " << nSig.getError() << " " << nBkg.getVal() << " " << nBkg.getError() << " " << fit_result->status() << std::endl;
    //     f_output.close();
    //   }
    //   std::cout << "   INFO: " << "=> for fit function " << bkg_function->getTitle() << " maximum signal is " << nSigMax << std::endl;
    //   // splusb_functions.push_back(ptrSplusB);
    // }
        
    // Recreate output graph file and open a new canvas for plotting
    TFile f_root_output((std::string("output/fit_functions_") + mass_kv.first + ".root").c_str(), "RECREATE");
    TCanvas canvas("canvas", "canvas", 800, 600);
    
    // Plot PDF fits to the data (MC in this case)
    RooPlot* frame = mass.frame();
    data.plotOn(frame);
    f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "data");
    for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
      parameter_sets.at(idx).restore_values();
      std::cout << "Restored nSig = " << nSig.getVal() << ", nBkg = " << nBkg.getVal() << ", mass_resonance = " << mass_resonance.getVal() << std::endl;
      // Background
      splusb_functions.at(idx)->plotOn(frame, RooFit::Components(bkg_functions.at(idx)->getTitle()), RooFit::LineColor(colours.at(idx)), RooFit::LineStyle(kDashed));
      f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "bkg_" + bkg_functions.at(idx)->getTitle());
      // Signal + background
      splusb_functions.at(idx)->plotOn(frame, RooFit::LineColor(colours.at(idx)));
      f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "splusb_" + bkg_functions.at(idx)->getTitle());
    }
    frame->Draw();
    for (int i = 0; i < frame->numItems(); ++i) {
      std::cout << "   INFO: " <<frame->getObject(i)->GetName() << std::endl;
    }
    TLegend legend(0.6, 0.7, 0.9, 0.9);
    legend.AddEntry(frame->findObject("signal_plus_novosibirsk_Norm[mass]"), "Novosibirsk", "L");
    // legend.AddEntry(frame->findObject("signal_plus_modified_cauchy_Norm[mass]"), "Cauchy", "L");
    // legend.AddEntry(frame->findObject("signal_plus_modified_gamma_Norm[mass]"), "Gamma", "L");
    legend.Draw();    
    canvas.Print((std::string("output/m_yyjj_") + mass_kv.first + "Mass.pdf").c_str());
    f_root_output.Close();
  }
  return 0;
}

//  INFO: ... resonance mass = 245 => nSig: 2.74822e-06, nBkg: 91286.7. Fit status = 0
//  INFO: ... resonance mass = 295 => nSig: 0.000239165, nBkg: 91286.7. Fit status = 0
//  INFO: ... resonance mass = 345 => nSig: 239.726, nBkg: 91049.4. Fit status = 0
//  INFO: ... resonance mass = 395 => nSig: 2.2286, nBkg: 91281.6. Fit status = 0
//  INFO: ... resonance mass = 445 => nSig: 2.37313e-05, nBkg: 91286.8. Fit status = 0
//  INFO: => for fit function novosibirsk maximum signal is 239.726 (with 91049.4 bkg events) at 345
//  INFO: Using fit function modified_gamma
//  INFO: ... resonance mass = 245 => nSig: 36.1847, nBkg: 91254. Fit status = 3
//  INFO: ... resonance mass = 295 => nSig: 100.397, nBkg: 91185.1. Fit status = 0
//  INFO: ... resonance mass = 345 => nSig: 399.611, nBkg: 90886.2. Fit status = 0
//  INFO: ... resonance mass = 395 => nSig: 248.82, nBkg: 91034.8. Fit status = 0
//  INFO: ... resonance mass = 445 => nSig: 107.264, nBkg: 91179.7. Fit status = 0
//  INFO: => for fit function modified_gamma maximum signal is 399.611 (with 90886.2 bkg events) at 345
//  INFO: novosibirsk: 239.726 / 91049.4
//  INFO: modified_gamma: 399.611 / 90886.2
