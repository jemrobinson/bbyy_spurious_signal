// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "SignalModel.h"
// STL
#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <utility>
#include <vector>
// ROOT and RooFit/RooStats
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRandom.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"

int main(int /*argc*/, char** /*argv*/)
{
  using namespace SpuriousSignal;

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gErrorIgnoreLevel = kBreak;

  RooRandom::randomGenerator()->SetSeed(20170720);

  // Define signal injection points
  std::vector<int> masses = {260, 275, 300, 325, 350, 400, 450, 500, 750, 1000};
  std::map< std::string, std::vector<double> > injected_signals = {{"ATLAS_Run_1_combination", {1.7, 1.7, 1.6, 1.6, 1.5, 0.8, 0.7, 0.6, 0.35, 0.01}},
                                                                   {"ATLAS_Run_2_bbyy", {7, 7, 6.1, 5.6, 5.1, 4.1, 4, 4, 4, 4}},
                                                                   {"ATLAS_Run_2_bbbb", {99, 99, 99, 99, 99, 1.48, 0.89, 0.73, 0.16, 0.11}},
                                                                   {"CMS_Run_2_bbyy", {1.23, 1.40, 1.04, 0.60, 0.80, 0.42, 0.48, 0.23, 0.10, 0.31}},
                                                                   {"Expected_bbyy", {1.15, 1.0, 0.9, 0.8, 0.7, 0.55, 0.5, 0.38, 0.18, 0.13}}};
                                                                  //  {"Small_xs", {0.00115, 0.0010, 0.0009, 0.0008, 0.0007, 0.00055, 0.0005, 0.00038, 0.00018, 0.00013}}};
  std::map< std::string, std::vector<std::string> > bkg_model_parameters = {{"novosibirsk", {"novosibirsk_width", "novosibirsk_peak", "novosibirsk_tail"}},
                                                                            {"exppoly1", {"exppoly1_p1"}}};

                                                                          // Construct mass and tag categories
  std::vector<std::string> mass_categories({"low", "high"});
  std::vector<std::string> tag_categories({"1", "2"});

  // Open workspace files
  TFile f_bkg("output/background_model_workspace.root", "READ");
  TFile f_signal("output/signal_model_workspace.root", "READ");

  // Recreate output text file
  std::string f_output_text("output/signal_bias.csv");
  std::ofstream f_output;
  f_output.open(f_output_text, std::ios::trunc);
  f_output.close();

  // Define weight parameter
  RooRealVar weight("weight", "event weight", -1e10, 1e10);

  // Iterate over mass and tag categories
  for (auto mass_category : mass_categories) {
    std::string fitfn_name(mass_category == "low" ? "novosibirsk" : "exppoly1");

    for (auto tag_category : tag_categories) {
      MSG_INFO("Loading workspaces for \033[1m" << mass_category << " mass " << tag_category << "-tag\033[0m category");

      // Load workspaces
      RooWorkspace* wk_bkg(0);
      f_bkg.GetObject(("background_model_" + mass_category + "Mass_" + tag_category + "tag").c_str(), wk_bkg);
      if (!wk_bkg) { MSG_ERROR("Could not open background workspace!"); }
      RooWorkspace* wk_sig(0);
      f_signal.GetObject(("signal_model_" + mass_category + "Mass_" + tag_category + "tag").c_str(), wk_sig);
      if (!wk_sig) { MSG_ERROR("Could not open signal workspace!"); }

      // Load background PDF
      RooDataSet* bkg_dataset = dynamic_cast<RooDataSet*>(wk_bkg->data("data"));
      MSG_INFO("Loaded background dataset with " << bkg_dataset->numEntries() << " entries, corresponding to " << bkg_dataset->sumEntries() << " events");
      wk_bkg->factory(("ExtendPdf::bkg_model(" + fitfn_name + ", norm[" + std::to_string(bkg_dataset->sumEntries()) + "])").c_str());
      RooDataHist* asimov_bkg = wk_bkg->pdf("bkg_model")->generateBinned(*wk_bkg->var("mass"), RooFit::ExpectedData());

      // Iterate over mass points
      for (auto resonance_mass : PlotStyle::resonance_masses(mass_category)) {
        // Load signal data
        std::string mX(std::to_string(resonance_mass));
        RooDataSet* ptr_signal_data = RooDataSet::read(("input/m_yyjj_Xhh_m" + mX + "_" + mass_category + "Mass_" + tag_category + "tag_tightIsolated.csv").c_str(), RooArgList(*wk_bkg->var("mass"), weight));
        MSG_INFO("... loading \033[1m" << resonance_mass << " GeV\033[0m dataset containing " << ptr_signal_data->numEntries() << " entries");

        // Load signal PDF
        RooAbsPdf* signal_model = wk_sig->pdf("signal_PDF");
        wk_sig->var("mass_resonance")->setVal(resonance_mass);
        wk_sig->var("mass_resonance")->setConstant(true);
        wk_sig->var("mass")->setRange(wk_bkg->var("mass")->getMin(), wk_bkg->var("mass")->getMax());

        int iSignal = std::distance(masses.begin(), std::find(masses.begin(), masses.end(), resonance_mass));
        for (auto& kv : injected_signals) {
          double injected_signal_pb(kv.second.at(iSignal));
          MSG_INFO("  ... using limits from \033[1m" << kv.first << "\033[0m");

          // Construct formula to calculate scaled weight for events
          RooFormulaVar weight_function("scaled_weight", "scaled_weight", ("weight * " + std::to_string(injected_signal_pb) + " / 5.0").c_str(), weight);
          RooDataSet signal_dataset_clone(*ptr_signal_data);
          RooRealVar* scaled_weight = dynamic_cast<RooRealVar*>(signal_dataset_clone.addColumn(weight_function));
          RooDataSet signal_dataset("data", "data", RooArgSet(*wk_bkg->var("mass"), *scaled_weight), RooFit::Import(signal_dataset_clone), RooFit::WeightVar(*scaled_weight), RooFit::Cut(("0.9 * " + mX + " < mass && mass < 1.1 * " + mX).c_str()));
          double nInjectedEvents(signal_dataset.sumEntries());
          MSG_INFO("  ... injecting " << injected_signal_pb << " pb signal => corresponding to " << nInjectedEvents << " events");

          // Combine Asimov background with signal
          RooDataHist combined_dataset(*asimov_bkg);
          combined_dataset.add(signal_dataset);

          // Construct combined PDF
          RooRealVar nSig("nSig", "number of signal events", nInjectedEvents, 0, 5 * combined_dataset.sumEntries());
          RooRealVar nBkg("nBkg", "number of background events", bkg_dataset->sumEntries(), 0.5 * bkg_dataset->sumEntries(), 1.5 * bkg_dataset->sumEntries());
          RooAddPdf combined_PDF(("signal_plus_" + fitfn_name).c_str(), ("signal_plus_" + fitfn_name).c_str(), RooArgList(*signal_model, *wk_bkg->pdf(fitfn_name.c_str())), RooArgList(nSig, nBkg));

          // Force the background parameters to be constant
          MSG_INFO("  ... setting " << bkg_model_parameters[fitfn_name].size() << " background model parameters as constant.");
          // for (auto parameter : bkg_model_parameters[fitfn_name]) {
          //   wk_bkg->var(parameter.c_str())->setConstant();
          // }

          // std::cout << "nBkg initial: " << nBkg.getVal() << " vs. " << bkg_dataset->sumEntries() << std::endl;
          // for (auto parameter : bkg_model_parameters[fitfn_name]) {
          //   std::cout << parameter << " initial: " << wk_bkg->var(parameter.c_str())->getVal() << std::endl;
          // }

          // Fit S+B PDF to combined data
          combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
          combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(false), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(false), RooFit::Minos(true), RooFit::PrintLevel(-1));

          // std::cout << "nBkg final: " << nBkg.getVal() << " vs. " << bkg_dataset->sumEntries() << std::endl;
          // std::cout << " -> bkg bias = " << 100 * (nBkg.getVal() - bkg_dataset->sumEntries()) / bkg_dataset->sumEntries() << "%" << std::endl;

          // std::cout << "nSig final: " << nSig.getVal() << " vs. " << nInjectedEvents << std::endl;
          // for (auto parameter : bkg_model_parameters[fitfn_name]) {
          //   std::cout << parameter << " final: " << wk_bkg->var(parameter.c_str())->getVal() << std::endl;
          // }

          // Plot output
          if (kv.first == "Expected_bbyy") {
            PlotStyle::EnsureAtlasStyle();
            RooPlot* frame = wk_bkg->var("mass")->frame();
            TLegend legend(0.5, 0.8, 0.93, 0.93);
            combined_dataset.plotOn(frame);
            legend.AddEntry((TGraph*)frame->getObject(frame->numItems() - 1), "Asimov bkg + injected signal", "P");
            combined_PDF.plotOn(frame);
            legend.AddEntry((TGraph*)frame->getObject(frame->numItems() - 1), ("S+B (" + PlotStyle::label(fitfn_name) + ") fit").c_str(), "L");
            legend.AddEntry((TObject*)(0), kv.first.c_str());
            TCanvas canvas("canvas", "canvas", 600, 600);
            frame->Draw();
            TLatex textBox; textBox.SetNDC(); textBox.SetTextFont(42);
            textBox.DrawLatex(0.5, 0.74, (mX + " GeV signal: " + PlotStyle::to_string(injected_signal_pb, 1) + " pb (hh)").c_str());
            legend.SetBorderSize(0);
            legend.SetFillStyle(0);
            legend.Draw();
            canvas.Print(("plots/signal_bias/signal_bias_" + mass_category + "Mass_" + tag_category + "tag_mX_" + mX + "_" + kv.first + ".pdf").c_str());
          }

          // Write text output
          std::ofstream f_output;
          f_output.open(f_output_text, std::ios::app);
          f_output << mass_category << " " << tag_category << " " << mX << " " << kv.first << " " << nInjectedEvents << " " << nSig.getVal() << " " << nSig.getError() << std::endl;
          f_output.close();
        }
      }
    }
  }
  f_bkg.Close();
  return 0;
}