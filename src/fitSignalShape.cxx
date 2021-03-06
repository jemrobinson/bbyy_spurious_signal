// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "SignalModel.h"
// STL
#include <string>
#include <utility>
#include <vector>
// ROOT and RooFit
#include "RooDataSet.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "TFile.h"

#include "RooPlot.h"
#include "TCanvas.h"

int main(int argc, char** argv)
{
  using namespace SpuriousSignal;

  // Determine which signal models to use
  std::vector<std::string> model_names;
  std::vector<std::string> input_args(argv, argv + argc);
  for (auto arg: input_args){
     if (arg == "EGE") { model_names.push_back("EGE"); }
     if (arg == "CBGA") { model_names.push_back("CBGA"); }
  }
  if (model_names.size() == 0) { model_names.push_back("EGE"); }

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gErrorIgnoreLevel = kBreak;

  // Construct mass and tag categories
  std::vector<std::string> mass_categories({"low", "high"});
  std::vector<std::string> tag_categories({"0", "1", "2"});

  // Define data parameters
  RooRealVar weight("weight", "event weight", -1e10, 1e10);

  // Create output file
  std::string output_file_name("output/signal_model_workspace.root");
  TFile f_output(output_file_name.c_str(), "RECREATE");
  f_output.Close();

  // Iterate over mass and tag categories
  for (auto mass_category : mass_categories) {
    for (auto tag_category : tag_categories) {
      // Define one mass variable per workspace
      RooRealVar mass("mass", "m_{yyjj}", 10, 10000, "GeV");

      // Load data for each mass point
      std::map<std::string, RooDataSet*> dataset_map;
      for (auto resonance_mass : PlotStyle::resonance_masses(mass_category)) {
        std::string mX(std::to_string(resonance_mass));
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooDataSet* ptr_raw_data = RooDataSet::read(("input/m_yyjj_Xhh_m" + mX + "_" + mass_category + "Mass_" + tag_category + "tag_tightIsolated.csv").c_str(), RooArgList(mass, weight));
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooDataSet* _data = new RooDataSet("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight), RooFit::Cut(("0.9 * " + mX + " < mass && mass < 1.1 * " + mX).c_str()));
        MSG_INFO("Loaded " << _data->numEntries() << " mX = " << resonance_mass << " events for " << mass_category << " mass, " << tag_category << "-tag category, corresponding to " << _data->sumEntries() << " data events");
        dataset_map[mX] = _data;
      }

      // Construct signal model
      MSG_INFO("Constructing simultaneous PDF");
      SignalModel model(mass_category, tag_category, model_names);
      model.build_simultaneous_PDF(mass);

      // Construct combined dataset in terms of  (mass, mass_point)
      RooDataSet data("data_combined", "data_combined", RooArgSet(mass, weight), RooFit::Index(*model.mass_points()), RooFit::Import(dataset_map), RooFit::WeightVar(weight));

      // Fit PDFs to data
      MSG_INFO("Preparing to fit PDFs to events");
      model.fit(data);

      // Plot output
      MSG_INFO("Preparing to plot results");
      model.plot();

      // Write model to output file
      model.write(output_file_name);
    }
  }
  return 0;
}