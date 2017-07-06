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

int main(int /*argc*/, char** /*argv*/)
{
  using namespace SpuriousSignal;

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  gErrorIgnoreLevel = kBreak;

  // Construct mass and tag categories
  std::vector<std::string> mass_categories({"low"});
  std::vector<std::string> tag_categories({"1"});
  // std::vector<std::string> mass_categories({"low", "high"});
  // std::vector<std::string> tag_categories({"0", "1", "2"});

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
      RooRealVar mass("mass", "m_{yyjj}", (mass_category == "low" ? 245 : 335), (mass_category == "low" ? 485 : 1140), "GeV");

      // Construct signal model
      MSG_INFO("Constructing simultaneous PDF");
      SignalModel model(mass_category, tag_category);
      model.build_simultaneous_PDF(mass);

      // Load data for each mass point
      std::map<std::string, RooDataSet*> dataset_map;
      for (auto resonance_mass : PlotStyle::resonance_masses(mass_category)) {
        std::string mX(std::to_string(resonance_mass));
        RooDataSet* ptr_raw_data = RooDataSet::read(("input/m_yyjj_Xhh_m" + mX + "_" + mass_category + "Mass_" + tag_category + "tag_tightIsolated.csv").c_str(), RooArgList(mass, weight));
        RooDataSet* _data = new RooDataSet("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));
        MSG_INFO("Loaded " << _data->numEntries() << " mX = " << resonance_mass << " events for " << mass_category << " mass, " << tag_category << "-tag category, corresponding to " << _data->sumEntries() << " data events");
        dataset_map[mX] = _data;
      }

      // Construct combined dataset in terms of  (mass, mass_point)
      RooDataSet data("data_combined", "data_combined", RooArgSet(mass, weight), RooFit::Index(*model.mass_points()), RooFit::Import(dataset_map), RooFit::WeightVar(weight));

      // Fit PDFs to data
      MSG_INFO("Preparing to fit PDF");
      model.fit(data);

      // Plot output
      MSG_INFO("Preparing to plot results");
      // model.plot();
      model.plot(dataset_map);

      // Write model to output file
      model.write(output_file_name);
    }
  }
}