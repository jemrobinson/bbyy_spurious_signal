// STL
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
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
  int PRINT_LEVEL = 1;
    
  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  gErrorIgnoreLevel = kBreak;
  
  // Setup colours and mass categories
  std::vector<int> colours = {kViolet, kGreen + 3, kBlue, kRed};
  // std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)},
  //                                                                 {"high", std::make_pair<int, int>(335, 1140)}};
  std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)}};

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
    // Define parameter ranges
    double mass_dn((mass_kv.first == "low") ? 220 : 0);
    double mass_up((mass_kv.first == "low") ? 320 : 1000);
    
    // Define data parameters
    RooRealVar mass("mass", "m_{yyjj}", mass_kv.second.first, mass_kv.second.second);
    RooRealVar weight("weight", "event weight", 1e-10, 1e10);

    // Get the data
    RooDataSet* ptr_raw_data = RooDataSet::read(("input/m_yyjj_SM_bkg_" + mass_kv.first + "Mass_0tag_tightIsolated.csv").c_str(), RooArgList(mass, weight));
    RooDataSet data("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));
    MSG_INFO("Loaded " << data.numEntries() << " events");

    // Construct vector of mass points
    int mass_point = mass_kv.second.first - MASS_STEP;
    std::vector<double> mass_points((mass_kv.second.second - mass_kv.second.first) / MASS_STEP + 1); 
    std::generate(mass_points.begin(), mass_points.end(), [&mass_point, &MASS_STEP]{ return mass_point += MASS_STEP; }); 
    MSG_INFO("Will perform fits for " << mass_points.size() << " mass points between " << mass_kv.second.first << " and " << mass_kv.second.second << " GeV");
    
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
    RooRealVar novosibirsk_peak("novosibirsk_peak", "peak of Novosibirsk", 270, mass_dn, mass_up);
    RooRealVar novosibirsk_tail("novosibirsk_tail", "tail of Novosibirsk", -1, -20, 20);
    RooRealVar novosibirsk_width("novosibirsk_width", "width of Novosibirsk", 30, 0, 500);
    RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", mass, novosibirsk_peak, novosibirsk_width, novosibirsk_tail);
    // // Low 245
    // peak of Novosibirsk [271.254 +/- 11.6352]
    // width of Novosibirsk [33.7789 +/- 8.80857]
    // tail of Novosibirsk [-1.26493 +/- 0.605308]
    // // High 335
    // peak of Novosibirsk [348.826 +/- 713.726]
    // width of Novosibirsk [72.8156 +/- 101.824]
    // tail of Novosibirsk [-0.965015 +/- 2.66738]
    // // High 635
    // peak of Novosibirsk [317.717 +/- 832.046]
    // width of Novosibirsk [78.9108 +/- 375.908]
    // tail of Novosibirsk [-0.879161 +/- 11.6147]
    // // High 935
    // peak of Novosibirsk [317.573 +/- 831.01]
    // width of Novosibirsk [78.8959 +/- 376.094]
    // tail of Novosibirsk [-0.87915 +/- 11.6104]
    
    // Modified Gamma
    RooRealVar gamma_Gamma("gamma_Gamma", "gamma of Gamma", 8, 0, 1000);
    RooRealVar gamma_beta0("gamma_beta0", "beta0 of Gamma", 0.07, 0, 100);
    RooRealVar gamma_beta1("gamma_beta1", "beta1 of Gamma", 0.01, -0.1, 0.1);
    RooRealVar gamma_mu("gamma_mu", "mean of Gamma", 270, mass_dn, mass_up);
    RooFormulaVar gamma_beta("gamma_beta", "@1 * exp(@2 * @0)", RooArgList(mass, gamma_beta0, gamma_beta1)); // exponential as we require beta > 0 everywhere
    RooGamma modified_gamma("modified_gamma", "modified_gamma", mass, gamma_Gamma, gamma_beta, gamma_mu);
    // // Low 245
    // gamma_Gamma	  = 8	 +/-  nan	(limited)
    // gamma_beta0	  = 50	 +/-  nan	(limited)
    // gamma_beta1	  = 0	 +/-  nan	(limited)
    // gamma_mu	  = 270	 +/-  nan	(limited)
    // // High 335
    // gamma of Gamma [13.6365 +/- 0.0103797]
    // mean of Gamma [17.2259 +/- 165.527]
    // beta0 of Gamma [0.111925 +/- 0.00788861]
    // beta1 of Gamma [0.054172 +/- 4.03044e-05]
    // // High 635
    // gamma of Gamma [13.6365 +/- 976.716]
    // mean of Gamma [17.2259 +/- 0.123598]
    // beta0 of Gamma [0.111925 +/- 0.00311876]
    // beta1 of Gamma [0.054172 +/- 0.00171638]
    // // High 935
    // gamma of Gamma [13.6365 +/- 976.672]
    // mean of Gamma [17.2259 +/- 0.123555]
    // beta0 of Gamma [0.111925 +/- 0.0031181]
    // beta1 of Gamma [0.054172 +/- 0.00171601]

        
    // Modified Cauchy
    RooRealVar cauchy_mu("cauchy_mu", "location of Cauchy", 270, mass_dn, mass_up);
    RooRealVar cauchy_s0("cauchy_s0", "s0 of Cauchy", 10, 0, 1000);
    RooRealVar cauchy_s1("cauchy_s1", "s1 of Cauchy", 0.01, -0.1, 0.1);
    RooFormulaVar cauchy_s("cauchy_s", "@1 * exp(@2 * @0)", RooArgList(mass, cauchy_s0, cauchy_s1)); // exponential as we require s > 0 everywhere
    RooGenericPdf modified_cauchy("modified_cauchy", "modified_cauchy", "TMath::CauchyDist(@0, @1, @2)", RooArgList(mass, cauchy_mu, cauchy_s));
    // // Low 245
    // cauchy_mu	  = 293.993	 +/-  1.69854	(limited)
    // cauchy_s0	  = 2.74147	 +/-  0.176684	(limited)
    // cauchy_s1	  = 0.00999998	 +/-  0.000369937	(limited)
    // // High 335
    // location of Cauchy [396.704 +/- 1649.82]
    // s0 of Cauchy [13.7118 +/- 90.0009]
    // s1 of Cauchy [0.0061196 +/- 0.0178356]
    // // High 635
    // location of Cauchy [341.904 +/- 1557.38]
    // s0 of Cauchy [81.4045 +/- 82.4699]
    // s1 of Cauchy [0.00620036 +/- 0.00417543]
    // cauchy_mu	  = 1.69941e-06	 +/-  396.122	(limited)
    // cauchy_s0	  = 100	 +/-  51.8684	(limited)
    // cauchy_s1	  = -0.00269309	 +/-  0.000508089	(limited)
    // // High 935
    // location of Cauchy [342.081 +/- 1467.87]
    // s0 of Cauchy [78.7587 +/- 77.1952]
    // s1 of Cauchy [0.00620047 +/- 0.0043522]
            
    // Number of signal and background events
    RooRealVar nSig("nSig", "number of signal events", 0, 10);
    RooRealVar nBkg("nBkg", "number of background events", 0, 10000);

    // Setup fit functions
    std::vector<RooAbsPdf*> bkg_functions;
    bkg_functions.push_back(&novosibirsk);
    bkg_functions.push_back(&modified_gamma);
    bkg_functions.push_back(&modified_cauchy);

    // Construct parameter sets that need to be remembered
    std::vector<ParameterSet> parameter_sets;
    parameter_sets.push_back(ParameterSet("Novosibirsk", {&novosibirsk_peak, &novosibirsk_width, &novosibirsk_tail, &nSig, &nBkg, &mass_resonance}));
    parameter_sets.push_back(ParameterSet("Modified Gamma", {&gamma_Gamma, &gamma_mu, &gamma_beta0, &gamma_beta1, &nSig, &nBkg, &mass_resonance}));
    parameter_sets.push_back(ParameterSet("Modified Cauchy", {&cauchy_mu, &cauchy_s0, &cauchy_s1, &nSig, &nBkg, &mass_resonance}));

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
    for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
      // Selection criterion
      double Z_spurious(0), Z_spurious_max(-99);
      MSG_INFO("Using fit function " << splusb_functions.at(idx)->getTitle());

      // Iterate over mass points
      for (int iMass = mass_kv.second.first; iMass < mass_kv.second.second; iMass += MASS_STEP) {
        mass_resonance.setVal(iMass); mass_resonance.setConstant(true); int fit_status(0);
        fit_status = splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(PRINT_LEVEL), RooFit::Save(true))->status();
        if (fit_status != 0) { MSG_ERROR("... simplex fit did not converge: status = " << fit_status); }
        MSG_DEBUG("... simplex: " << parameter_sets.at(idx));
        fit_status = splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(true), RooFit::PrintLevel(PRINT_LEVEL), RooFit::Save(true))->status();
        if (fit_status != 0) { MSG_ERROR("... migrad + improve + hesse + minos fit did not converge: status = " << fit_status); }
        MSG_DEBUG("... migrad+improve+hesse+minos: " << parameter_sets.at(idx));
        // Record these parameters if Z > max(Z) so far
        Z_spurious = nSig.getVal() / nSig.getError(); // (N_S +/- errS_SumW2_true) / (errS_SumW2_false)
        MSG_INFO("... resonance mass = " << mass_resonance.getVal() << " => nSig: " << nSig.getVal() << " +/- " << nSig.getError() << ", nBkg: " << nBkg.getVal() << " +/- " << nBkg.getError() << " [" << fit_status << "] => Z = " << Z_spurious);
        if (Z_spurious > Z_spurious_max) {
          Z_spurious_max = Z_spurious;
          MSG_INFO("Recording new max(z) = " << Z_spurious_max << " at " << mass_resonance.getVal());
          parameter_sets.at(idx).record_values();
        }
        // Write parameters to output file
        f_output.open("output/spurious_signal_" + mass_kv.first + ".csv", std::ios::app);
        f_output << bkg_functions.at(idx)->getTitle() << " " << mass_resonance.getVal() << " " << nSig.getVal() << " " << nSig.getError() << " " << nBkg.getVal() << " " << nBkg.getError() << " " << fit_status << " " << Z_spurious << std::endl;
        f_output.close();
      }
    }

    // Recreate output graph file and open a new canvas for plotting
    TFile f_root_output((std::string("output/fit_functions_") + mass_kv.first + ".root").c_str(), "RECREATE");
    TCanvas canvas("canvas", "canvas", 800, 600);
    
    // Plot PDF fits to the data (MC in this case)
    RooPlot* frame = mass.frame();
    data.plotOn(frame);
    f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "data");
    for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
      parameter_sets.at(idx).restore_values();
      MSG_INFO("Restored " << bkg_functions.at(idx)->getTitle() << " parameters at resonant mass = " << mass_resonance.getVal() << ": nSig = " << nSig.getVal() << ", nBkg = " << nBkg.getVal());
      // Background
      splusb_functions.at(idx)->plotOn(frame, RooFit::Components(bkg_functions.at(idx)->getTitle()), RooFit::LineColor(colours.at(idx)), RooFit::LineStyle(kDashed));
      f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "bkg_" + bkg_functions.at(idx)->getTitle());
      // Signal + background
      splusb_functions.at(idx)->plotOn(frame, RooFit::LineColor(colours.at(idx)));
      f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "splusb_" + bkg_functions.at(idx)->getTitle());
    }
    // for (int i = 0; i < frame->numItems(); ++i) { MSG_INFO(frame->getObject(i)->GetName()); }
    frame->Draw();
    TLegend legend(0.6, 0.7, 0.9, 0.9);
    legend.AddEntry(frame->findObject("signal_plus_novosibirsk_Norm[mass]"), "Novosibirsk", "L");
    legend.AddEntry(frame->findObject("signal_plus_modified_gamma_Norm[mass]"), "Gamma", "L");
    legend.AddEntry(frame->findObject("signal_plus_modified_cauchy_Norm[mass]"), "Cauchy", "L");
    legend.Draw();    
    canvas.Print((std::string("output/m_yyjj_") + mass_kv.first + "Mass.pdf").c_str());
    MSG_INFO("Saved " + mass_kv.first + " mass plot");
    f_root_output.Close();
  }
  return 0;
}

