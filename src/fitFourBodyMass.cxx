// STL
#include <fstream>
#include <string>
#include <utility>
#include <vector>
// #include <algorithm>
// // ROOT
// #include "TCanvas.h"
// #include "TError.h"
#include "TFile.h"
// #include "TLegend.h"
// #include "TMath.h"
#include "TStyle.h"
// // RooFit
#include "RooAddPdf.h"
#include "RooConstVar.h"
// #include "RooArgList.h"
// // #include "RooArgusBG.h"
#include "RooCBShape.h"
// // #include "RooCFunction3Binding.h"
// #include "RooDataSet.h"
// // #include "RooExponential.h"
// #include "RooFFTConvPdf.h"
#include "RooFormulaVar.h"
// #include "RooGamma.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
// #include "RooLandau.h"
// #include "RooLognormal.h"
// #include "RooMsgService.h"
// #include "RooNLLVar.h"
#include "RooNovosibirsk.h"
// #include "RooPlot.h"
// #include "RooPlot.h"
// #include "RooPoisson.h"
// #include "RooRealVar.h"
// #include "RooConstVar.h"
// Local
#include "ParameterSet.h"
#include "Logger.h"
#include "FitMassPoint.h"
#include "SignalModel.h"

std::map<std::string, std::vector<std::string> > process_args(int argc, char** argv) {
  // Read in command line arguments
  std::string current_argument("");
  std::map<std::string, std::vector<std::string> > output_map;
  for (auto parameter : std::vector<std::string>(argv, argv+argc)) {
    if (parameter.substr(0, 2) == "--") {
      current_argument = parameter.substr(2, parameter.size());
      output_map[current_argument] = std::vector<std::string>(0);
    } else {
      if (!current_argument.empty()) { output_map[current_argument].push_back(parameter); }
    }
  }
  // Add defaults
  if (output_map.find("mass") == output_map.end()) { output_map["mass"] = std::vector<std::string>({"low", "high"}); }
  if (output_map.find("tag") == output_map.end()) { output_map["tag"] = std::vector<std::string>({"0", "1", "2"}); }
  return output_map;
}


int main(int argc, char** argv) {
  using namespace SpuriousSignal;
  // Read in parse arguments
  std::map< std::string, std::vector<std::string> > args = process_args(argc, argv);
  std::vector<std::string> mass_categories = args["mass"];
  std::vector<std::string> tag_categories = args["tag"];

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  gErrorIgnoreLevel = kBreak;
  gStyle->SetOptTitle(0);

  int MASS_STEP = 10;

  // Setup colours and mass categories
  std::vector<int> colours = {kViolet, kGreen + 3, kBlue, kRed, kMagenta, kCyan};

  // Set up signal model
  SignalModel sm("/afs/cern.ch/user/a/andari/public/Leo/SigParam_low_cbga_cut5/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root", "signalWS");
  RooConstVar a_muGANom_SM_c2("a_muGANom_SM_c2", "a_muGANom_SM_c2", sm.m_a_muGANom_SM_c2);
  RooConstVar b_muGANom_SM_c2("b_muGANom_SM_c2", "b_muGANom_SM_c2", sm.m_b_muGANom_SM_c2);
  RooConstVar c_muGANom_SM_c2("c_muGANom_SM_c2", "c_muGANom_SM_c2", sm.m_c_muGANom_SM_c2);
  RooConstVar a_sigmaGANom_SM_c2("a_sigmaGANom_SM_c2", "a_sigmaGANom_SM_c2", sm.m_a_sigmaGANom_SM_c2);
  RooConstVar b_sigmaGANom_SM_c2("b_sigmaGANom_SM_c2", "b_sigmaGANom_SM_c2", sm.m_b_sigmaGANom_SM_c2);
  // Crystal Ball
  RooConstVar a_muCBNom_SM_c2("a_muCBNom_SM_c2", "a_muCBNom_SM_c2", sm.m_a_muCBNom_SM_c2);
  RooConstVar b_muCBNom_SM_c2("b_muCBNom_SM_c2", "b_muCBNom_SM_c2", sm.m_b_muCBNom_SM_c2);
  RooConstVar c_muCBNom_SM_c2("c_muCBNom_SM_c2", "c_muCBNom_SM_c2", sm.m_c_muCBNom_SM_c2);
  RooConstVar a_sigmaCBNom_SM_c2("a_sigmaCBNom_SM_c2", "a_sigmaCBNom_SM_c2", sm.m_a_sigmaCBNom_SM_c2);
  RooConstVar b_sigmaCBNom_SM_c2("b_sigmaCBNom_SM_c2", "b_sigmaCBNom_SM_c2", sm.m_b_sigmaCBNom_SM_c2);
  RooConstVar nCB_SM_c2("nCB_SM_c2", "nCB_SM_c2", sm.m_nCB_SM_c2);
  RooConstVar a_alphaCB_SM_c2("a_alphaCB_SM_c2", "a_alphaCB_SM_c2", sm.m_a_alphaCB_SM_c2);
  RooConstVar b_alphaCB_SM_c2("b_alphaCB_SM_c2", "b_alphaCB_SM_c2", sm.m_b_alphaCB_SM_c2);
  // Combination
  RooConstVar fracCB_SM_c2("fracCB_SM_c2", "fracCB_SM_c2", sm.m_fracCB_SM_c2);

  for (auto mass_category : mass_categories) {
    // Define mass ranges
    std::pair<int, int> mass_range = (mass_category == "low" ? std::make_pair<int, int>(245, 485) : std::make_pair<int, int>(335, 1140));
    std::pair<int, int> peak_range = (mass_category == "low" ? std::make_pair<int, int>(260, 280) : std::make_pair<int, int>(0, 1000));

    // Define data parameters
    RooRealVar mass("mass", "m_{yyjj}", mass_range.first, mass_range.second, "GeV");
    RooRealVar weight("weight", "event weight", 1e-10, 1e10);

    // Iterate over tag categories
    for (auto tag_category : tag_categories) {
      // Get the data
      RooDataSet* ptr_raw_data = RooDataSet::read(("input/m_yyjj_SM_bkg_" + mass_category + "Mass_" + tag_category + "tag_tightIsolated.csv").c_str(), RooArgList(mass, weight));
      RooDataSet data("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));
      MSG_INFO("Loaded " << data.numEntries() << " events for " << tag_category << "-tag category, corresponding to " << data.sumEntries() << " data events");

      // Construct vectors of mass points
      int mass_point(mass_range.first - MASS_STEP);
      std::vector<double> mass_points((mass_range.second - mass_range.first) / MASS_STEP + 1);
      std::generate(mass_points.begin(), mass_points.end(), [&mass_point, &MASS_STEP]{ return mass_point += MASS_STEP; });

      // Signal: CB+Gaus
      RooRealVar mass_resonance("mass_resonance", "mass of resonance", 300, 0, 2000);
      RooFormulaVar signal_gaus_mean("signal_gaus_mean", "a_muGANom_SM_c2 + b_muGANom_SM_c2 * (mass_resonance / 100. - 1) + c_muGANom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muGANom_SM_c2, b_muGANom_SM_c2, c_muGANom_SM_c2, mass_resonance));
      RooFormulaVar signal_gaus_sigma("signal_gaus_sigma", "a_sigmaGANom_SM_c2 + b_sigmaGANom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaGANom_SM_c2, b_sigmaGANom_SM_c2, mass_resonance));
      RooGaussian signal_gaus("signal_gaus", "signal gaus", mass, signal_gaus_mean, signal_gaus_sigma);
      RooFormulaVar signal_CB_mean("signal_CB_mean", "a_muCBNom_SM_c2 + b_muCBNom_SM_c2 * (mass_resonance / 100. - 1) + c_muCBNom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muCBNom_SM_c2, b_muCBNom_SM_c2, c_muCBNom_SM_c2, mass_resonance));
      RooFormulaVar signal_CB_sigma("signal_CB_sigma", "a_sigmaCBNom_SM_c2 + b_sigmaCBNom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaCBNom_SM_c2, b_sigmaCBNom_SM_c2, mass_resonance));
      RooFormulaVar signal_CB_alpha("signal_CB_alpha", "a_alphaCB_SM_c2 + b_alphaCB_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_alphaCB_SM_c2, b_alphaCB_SM_c2, mass_resonance));
      RooFormulaVar signal_CB_n("signal_CB_n", "nCB_SM_c2", RooArgList(nCB_SM_c2));
      RooCBShape signal_CB("signal_CB", "signal_CB", mass, signal_CB_mean, signal_CB_sigma, signal_CB_alpha, signal_CB_n);
      RooFormulaVar signal_frac_CB("signal_frac_CB", "fracCB_SM_c2", RooArgList(fracCB_SM_c2));
      RooAddPdf signal_PDF("signal_PDF", "CB + Gaussian", RooArgList(signal_CB, signal_gaus), signal_frac_CB);

      // Novosibirsk
      RooRealVar novosibirsk_peak("novosibirsk_peak", "peak of Novosibirsk", 270, peak_range.first, peak_range.second);
      RooRealVar novosibirsk_tail("novosibirsk_tail", "tail of Novosibirsk", -1, -20, 20);
      RooRealVar novosibirsk_width("novosibirsk_width", "width of Novosibirsk", 30, 0, 500);
      RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", mass, novosibirsk_peak, novosibirsk_width, novosibirsk_tail);

      // Modified Gamma
      RooRealVar gamma_alpha0("gamma_alpha0", "alpha0 of Gamma", (mass_category == "low" ? 1.5 : 0.08), 0, 1000);
      RooRealVar gamma_alpha1("gamma_alpha1", "alpha1 of Gamma", 0.003, -0.01, 0.01);
      RooFormulaVar gamma_alpha("gamma_alpha", "gamma_alpha0 + gamma_alpha1 * mass", RooArgList(mass, gamma_alpha0, gamma_alpha1));
      RooRealVar gamma_theta0("gamma_theta0", "theta0 of Gamma", (mass_category == "low" ? 0.07 : 190), 0, 1000);
      RooRealVar gamma_theta1("gamma_theta1", "theta1 of Gamma", 0.1, -1, 1);
      RooFormulaVar gamma_theta("gamma_theta", "gamma_theta0 + gamma_theta1 * mass", RooArgList(mass, gamma_theta0, gamma_theta1));
      RooRealVar gamma_mu("gamma_mu", "minimum of Gamma", mass_range.first, 0, mass_range.second);
      RooGenericPdf modified_gamma("modified_gamma", "modified_gamma","TMath::GammaDist(mass, gamma_alpha, gamma_mu, gamma_theta)", RooArgList(mass, gamma_alpha, gamma_theta, gamma_mu));

      // Modified Landau
      RooRealVar landau_mean("landau_mean", "mean of Landau", 270, peak_range.first, peak_range.second);
      RooRealVar landau_sigma0("landau_sigma0", "sigma0 of Landau", (mass_category == "low" ? 6.5 : 100), 0, 1000);
      RooRealVar landau_sigma1("landau_sigma1", "sigma1 of Landau", (mass_category == "low" ? 0.05 : 0), -1, 1);
      RooFormulaVar landau_sigma("landau_sigma", "landau_sigma0 + landau_sigma1 * mass", RooArgList(mass, landau_sigma0, landau_sigma1));
      RooGenericPdf modified_landau("modified_landau", "modified_landau", "TMath::Landau(mass, landau_mean, landau_sigma)", RooArgList(mass, landau_mean, landau_sigma));

      // Setup fit functions
      std::vector<RooAbsPdf*> bkg_functions;
      bkg_functions.push_back(&novosibirsk);
      bkg_functions.push_back(&modified_gamma);
      bkg_functions.push_back(&modified_landau);

      // Number of signal and background events
      RooRealVar nSig("nSig", "number of signal events", 0, -(0.5 * data.sumEntries()), (0.5 * data.sumEntries()));
      RooRealVar nBkg("nBkg", "number of background events", 1, 0, 1.5 * data.sumEntries());

      // Construct parameter sets that need to be remembered
      std::vector<ParameterSet> parameter_sets;
      parameter_sets.push_back(ParameterSet("Novosibirsk", {&novosibirsk_peak, &novosibirsk_width, &novosibirsk_tail, &nSig, &nBkg, &mass_resonance}));
      parameter_sets.push_back(ParameterSet("Modified Gamma", {&gamma_alpha0, &gamma_alpha1, &gamma_theta0, &gamma_theta1, &gamma_mu, &nSig, &nBkg, &mass_resonance}));
      parameter_sets.push_back(ParameterSet("Modified Landau", {&landau_mean, &landau_sigma0, &landau_sigma1, &nSig, &nBkg, &mass_resonance}));

      // Recreate output ROOT file
      std::string f_output_ROOT("output/fit_functions_" + mass_category + "Mass_" + tag_category + "tag.root");
      TFile f_ROOT(f_output_ROOT.c_str(), "RECREATE");
      f_ROOT.Close();

      // Recreate output text file
      std::string f_output_text("output/spurious_signal_" + mass_category + "Mass_" + tag_category + "tag.csv");
      std::ofstream f_text;
      f_text.open(f_output_text, std::ios::trunc);
      f_text.close();

      // Do background-only fits
      MSG_INFO("Performing background-only fits for " << bkg_functions.size() << " fit functions.");
      FitMassPoint fits_bkg_only(data, bkg_functions, mass_category, tag_category);
      fits_bkg_only.fit();
      fits_bkg_only.plot(mass.frame(), -1);
      fits_bkg_only.write(f_output_ROOT, f_output_text);
      for (auto& parameter_set : parameter_sets) { parameter_set.record_values(); }

      // Do S+B fits for different backgrounds
      MSG_INFO("Performing signal + background fits for " << bkg_functions.size() << " fit functions.");
      std::vector<RooAbsPdf*> splusb_functions;
      for (auto bkg_function : bkg_functions) {
        splusb_functions.push_back(new RooAddPdf("signal_plus_" + bkg_function->getTitle(), "signal + " + bkg_function->getTitle(), RooArgList(signal_PDF, *bkg_function), RooArgList(nSig, nBkg)));
      }
      FitMassPoint fits_splusb(data, splusb_functions, mass_category, tag_category);
      for (auto mass_point : mass_points) {
        // Set mass and restore bkg parameters to best bkg-only fit
        for (auto& parameter_set : parameter_sets) { parameter_set.restore_values(); }
        mass_resonance.setVal(mass_point); mass_resonance.setConstant(true);
        // Fit, plot and output results
        fits_splusb.fit();
        fits_splusb.plot(mass.frame(), mass_point);
        fits_splusb.write(f_output_ROOT, f_output_text);
      }
    }
  }
  return 0;
}

