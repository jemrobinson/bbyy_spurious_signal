// Local
#include "ParameterSet.h"
#include "Logger.h"
#include "PDFModelFitter.h"
#include "PlotStyle.h"
#include "SignalModel.h"
// STL
#include <fstream>
#include <string>
#include <utility>
#include <vector>
// ROOT and RooFit
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooNovosibirsk.h"
#include "RooWorkspace.h"
#include "TFile.h"



std::map<std::string, std::vector<std::string> > process_args(int argc, char** argv)
{
  // Read in command line arguments
  std::string current_argument("");
  std::map<std::string, std::vector<std::string> > output_map;
  // Iterate over arguments (as vector)
  for (auto parameter : std::vector<std::string>(argv, argv + argc)) {
    if (parameter.substr(0, 2) == "--") {
      current_argument = parameter.substr(2, parameter.size());
      output_map[current_argument] = std::vector<std::string>(0);
    } else if (!current_argument.empty()) {
      output_map[current_argument].push_back(parameter);
    }
  }
  // Add defaults
  if (output_map.find("mass") == output_map.end()) { output_map["mass"] = std::vector<std::string>({"low", "high"}); }
  if (output_map.find("tag") == output_map.end()) { output_map["tag"] = std::vector<std::string>({"0", "1", "2"}); }
  // Return map
  return output_map;
}


int main(int argc, char** argv)
{
  using namespace SpuriousSignal;
  int MASS_STEP = 10;

  // Read in parse arguments
  std::map< std::string, std::vector<std::string> > args = process_args(argc, argv);
  std::vector<std::string> mass_categories = args["mass"];
  std::vector<std::string> tag_categories = args["tag"];

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  gErrorIgnoreLevel = kBreak;

  // Open workspace file
  TFile f_input("output/signal_model_workspace.root", "READ");

  // Define weight parameter
  RooRealVar weight("weight", "event weight", -1e10, 1e10);

  for (auto mass_category : mass_categories) {
    // Define mass ranges for this category
    std::pair<int, int> mass_range = (mass_category == "low" ? std::make_pair<int, int>(245, 485) : std::make_pair<int, int>(335, 1140));
    std::pair<int, int> peak_range = (mass_category == "low" ? std::make_pair<int, int>(260, 280) : std::make_pair<int, int>(0, 1000));

    // Iterate over tag categories
    for (auto tag_category : tag_categories) {
      // Load the workspace
      RooWorkspace* wk(0);
      f_input.GetObject(("signal_model_" + mass_category + "Mass_" + tag_category + "tag").c_str(), wk);
      if (!wk) { MSG_ERROR("Could not open workspace!"); }
      wk->var("mass")->setRange(mass_range.first, mass_range.second);

      // Get the data
      RooDataSet* ptr_raw_data = RooDataSet::read(("input/m_yyjj_SM_bkg_" + mass_category + "Mass_" + tag_category + "tag_tightIsolated.csv").c_str(), RooArgList(*wk->var("mass"), weight));
      RooDataSet data("data", "data", RooArgSet(*wk->var("mass"), weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));
      MSG_INFO("Loaded " << data.numEntries() << " events for " << tag_category << "-tag category, corresponding to " << data.sumEntries() << " data events");
      // Construct vectors of mass points: either full range or specified set
      bool appendToFile(false);
      std::vector<double> mass_points;
      if (args.find("mX") != args.end()) {
        for (auto mass_point : args["mX"]) { mass_points.push_back(std::stoi(mass_point)); }
        appendToFile = true;
      } else {
        int mass_point(mass_range.first - MASS_STEP);
        mass_points.resize((mass_range.second - mass_range.first) / MASS_STEP + 1);
        std::generate(mass_points.begin(), mass_points.end(), [&mass_point, &MASS_STEP] { return mass_point += MASS_STEP; });
      }

      // Novosibirsk: 0 if  -tail < width / ( peak - x ) => (peak - x) > width / (-tail)
      RooRealVar novosibirsk_peak("novosibirsk_peak", "peak of Novosibirsk", 270, peak_range.first, peak_range.second);
      RooRealVar novosibirsk_tail("novosibirsk_tail", "tail of Novosibirsk", -1, -2, 10);
      RooRealVar novosibirsk_width("novosibirsk_width", "width of Novosibirsk", 30, 0, 500);
      RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", *wk->var("mass"), novosibirsk_peak, novosibirsk_width, novosibirsk_tail);
      // Modified Gamma
      RooRealVar gamma_alpha0("gamma_alpha0", "alpha0 of Gamma", (mass_category == "low" ? 1.5 : 0.08), 0, 1000);
      RooRealVar gamma_alpha1("gamma_alpha1", "alpha1 of Gamma", 0.003, -0.001, 0.01);
      RooFormulaVar gamma_alpha("gamma_alpha", "gamma_alpha0 + gamma_alpha1 * mass", RooArgList(*wk->var("mass"), gamma_alpha0, gamma_alpha1));
      RooRealVar gamma_theta0("gamma_theta0", "theta0 of Gamma", (mass_category == "low" ? 0.07 : 190), 0, 1000);
      RooRealVar gamma_theta1("gamma_theta1", "theta1 of Gamma", 0.1, -0.1, 1);
      RooFormulaVar gamma_theta("gamma_theta", "gamma_theta0 + gamma_theta1 * mass", RooArgList(*wk->var("mass"), gamma_theta0, gamma_theta1));
      RooRealVar gamma_mu("gamma_mu", "minimum of Gamma", mass_range.first, 0, mass_range.second + 100);
      RooGenericPdf modified_gamma("modified_gamma", "modified_gamma", "TMath::GammaDist(mass, gamma_alpha, gamma_mu, gamma_theta)", RooArgList(*wk->var("mass"), gamma_alpha, gamma_theta, gamma_mu));
      // Modified Landau
      RooRealVar landau_mean("landau_mean", "mean of Landau", 270, peak_range.first, peak_range.second);
      RooRealVar landau_sigma0("landau_sigma0", "sigma0 of Landau", (mass_category == "low" ? 6.5 : 100), 0, 1000);
      RooRealVar landau_sigma1("landau_sigma1", "sigma1 of Landau", (mass_category == "low" ? 0.05 : 0), -1, 1);
      RooFormulaVar landau_sigma("landau_sigma", "landau_sigma0 + landau_sigma1 * mass", RooArgList(*wk->var("mass"), landau_sigma0, landau_sigma1));
      RooGenericPdf modified_landau("modified_landau", "modified_landau", "TMath::Landau(mass, landau_mean, landau_sigma)", RooArgList(*wk->var("mass"), landau_mean, landau_sigma));
      // Exponential polynominal: degree-1
      RooRealVar exppoly1_p1("exppoly1_p1", "exppoly1_p1", -6e-03, -1.0, 0.0);
      RooGenericPdf exppoly1("exppoly1", "exppoly1", "TMath::Exp(exppoly1_p1 * mass)", RooArgList(*wk->var("mass"), exppoly1_p1));
      // Exponential polynominal: degree-2
      RooRealVar exppoly2_p1("exppoly2_p1", "exppoly2_p1", -4e-06, -1.0, 0.0);
      RooGenericPdf exppoly2("exppoly2", "exppoly2", "TMath::Exp(exppoly2_p1 * mass * mass)", RooArgList(*wk->var("mass"), exppoly2_p1));
      // Inverse polynomial: degree-2
      RooRealVar invpoly2_p1("invpoly2_p1", "invpoly2_p1", 1e19 , 1e6, 1e35);
      RooGenericPdf invpoly2("invpoly2", "invpoly2", "1 + invpoly2_p1 / (mass * mass)", RooArgList(*wk->var("mass"), invpoly2_p1));
      // Inverse polynomial: degree-2
      RooRealVar invpoly3_p1("invpoly3_p1", "invpoly3_p1", 1e15, 1e9, 1e40);
      RooGenericPdf invpoly3("invpoly3", "invpoly3", "1 + invpoly3_p1 / (mass * mass * mass)", RooArgList(*wk->var("mass"), invpoly3_p1));

      // Setup fit functions
      std::vector<RooAbsPdf*> bkg_functions;
      if (mass_category == "low") {
        bkg_functions.push_back(&novosibirsk);
        bkg_functions.push_back(&modified_gamma);
        bkg_functions.push_back(&modified_landau);
      } else if (mass_category == "high") {
        bkg_functions.push_back(&exppoly1);
        bkg_functions.push_back(&exppoly2);
        bkg_functions.push_back(&invpoly2);
        bkg_functions.push_back(&invpoly3);
      }

      // Number of signal and background events
      RooRealVar nSig("nSig", "number of signal events", 0, -(0.5 * data.sumEntries()), (0.5 * data.sumEntries()));
      RooRealVar nBkg("nBkg", "number of background events", data.sumEntries(), 0, 1.5 * data.sumEntries());

      // Construct parameter sets that need to be remembered
      std::vector<ParameterSet> parameter_sets;
      if (mass_category == "low") {
        parameter_sets.push_back(ParameterSet(PlotStyle::label("novosibirsk"), {&novosibirsk_peak, &novosibirsk_width, &novosibirsk_tail, &nSig, &nBkg}));
        parameter_sets.push_back(ParameterSet(PlotStyle::label("modified_gamma"), {&gamma_alpha0, &gamma_alpha1, &gamma_theta0, &gamma_theta1, &gamma_mu, &nSig, &nBkg}));
        parameter_sets.push_back(ParameterSet(PlotStyle::label("modified_landau"), {&landau_mean, &landau_sigma0, &landau_sigma1, &nSig, &nBkg}));
      } else if (mass_category == "high") {
        parameter_sets.push_back(ParameterSet(PlotStyle::label("exppoly1"), {&exppoly1_p1, &nSig, &nBkg}));
        parameter_sets.push_back(ParameterSet(PlotStyle::label("exppoly2"), {&exppoly2_p1, &nSig, &nBkg}));
        parameter_sets.push_back(ParameterSet(PlotStyle::label("invpoly2"), {&invpoly2_p1, &nSig, &nBkg}));
        parameter_sets.push_back(ParameterSet(PlotStyle::label("invpoly3"), {&invpoly3_p1, &nSig, &nBkg}));
      }

      // Recreate output ROOT file
      std::string fileSuffix(appendToFile ? "_mX" + args["mX"][0] : "");
      std::string f_output_ROOT("output/fit_functions_" + mass_category + "Mass_" + tag_category + "tag" + fileSuffix + ".root");
      TFile f_ROOT(f_output_ROOT.c_str(), (appendToFile ? "WRITE" : "RECREATE"));
      f_ROOT.Close();

      // Recreate output text file
      std::string f_output_text("output/csv/spurious_signal_" + mass_category + "Mass_" + tag_category + "tag" + fileSuffix + ".csv");
      std::ofstream f_text;
      f_text.open(f_output_text, (appendToFile ? std::ios::app : std::ios::trunc));
      f_text.close();

      // Do background-only fits
      MSG_INFO("Performing background-only fits for " << bkg_functions.size() << " fit functions.");
      PDFModelFitter fits_bkg_only(data, bkg_functions, mass_category, tag_category, true);
      fits_bkg_only.fit();
      fits_bkg_only.plot(wk->var("mass")->frame(), -1);
      fits_bkg_only.write(f_output_ROOT, f_output_text);

      // Record values for each set of parameters
      for (auto& parameter_set : parameter_sets) {
        parameter_set.record_values();
      }

      // Load signal model
      RooAbsPdf* signal_model = wk->pdf("signal_PDF");

      // Do S+B fits for different backgrounds
      MSG_INFO("Performing signal + background fits for " << bkg_functions.size() << " fit functions.");

      // Construct S+B PDFs
      std::vector<RooAbsPdf*> splusb_functions;
      for (auto bkg_function : bkg_functions) {
        splusb_functions.push_back(new RooAddPdf("signal_plus_" + bkg_function->getTitle(), "signal + " + bkg_function->getTitle(), RooArgList(*signal_model, *bkg_function), RooArgList(nSig, nBkg)));
      }

      PDFModelFitter fits_splusb(data, splusb_functions, mass_category, tag_category, true);
      for (auto mass_point : mass_points) {
        MSG_INFO("Fitting mass point \033[1m" << mass_point << "\033[0m GeV.");

        // Set mass and restore bkg parameters to best bkg-only fit
        for (auto& parameter_set : parameter_sets) { parameter_set.restore_values(); }
        wk->var("mass_resonance")->setVal(mass_point);
        wk->var("mass_resonance")->setConstant(true);

        // Fit, plot and output results
        fits_splusb.fit();
        fits_splusb.plot(wk->var("mass")->frame(), mass_point);
        fits_splusb.write(f_output_ROOT, f_output_text);
      }
    }
  }
  f_input.Close();
  return 0;
}

