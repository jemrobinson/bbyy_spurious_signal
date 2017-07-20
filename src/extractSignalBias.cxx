// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "SignalModel.h"
// STL
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
// // ROOT and RooFit/RooStats
#include "RooMsgService.h"
#include "RooWorkspace.h"
#include "RooRandom.h"
#include "TFile.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TLatex.h"

int main(int /*argc*/, char** /*argv*/)
{
  using namespace SpuriousSignal;

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gErrorIgnoreLevel = kBreak;

  RooRandom::randomGenerator()->SetSeed(20170720);

  // Define signal injection points
  std::vector<double> injected_signals_pb = {5, 1, 0.5, 0.1, 0.01};

  // Construct mass and tag categories
  std::vector<std::string> mass_categories({"low", "high"});
  std::vector<std::string> tag_categories({"0", "1", "2"});

  // Number of iterations
  unsigned int nIterations = 10; //100;

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
    std::string fitfn_name(mass_category == "low" ? "novosibirsk" : "invpoly3");

    for (auto tag_category : tag_categories) {
      MSG_INFO("Loading workspaces for \033[1m" << mass_category << " mass " << tag_category << "-tag\033[0m category");

      // Load workspaces
      RooWorkspace* wk_bkg(0);
      f_bkg.GetObject(("background_model_" + mass_category + "Mass_" + tag_category + "tag").c_str(), wk_bkg);
      if (!wk_bkg) { MSG_ERROR("Could not open workspace!"); }
      RooWorkspace* wk_sig(0);
      f_signal.GetObject(("signal_model_" + mass_category + "Mass_" + tag_category + "tag").c_str(), wk_sig);
      if (!wk_sig) { MSG_ERROR("Could not open workspace!"); }

      // Load background PDF
      RooDataSet* bkg_dataset = dynamic_cast<RooDataSet*>(wk_bkg->data("data"));
      RooAbsPdf* bkg_PDF = wk_bkg->pdf(fitfn_name.c_str());
      MSG_INFO("Loaded background dataset with " << bkg_dataset->sumEntries() << " entries, corresponding to " << bkg_dataset->sumEntries() << " events");

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

        for (auto injected_signal_pb : injected_signals_pb) {
          // Construct formula to calculate (fake) weight for events
          RooFormulaVar weight_function("scaled_weight", "scaled_weight", ("weight * " + std::to_string(injected_signal_pb) + " / 5.0").c_str(), weight);
          RooDataSet signal_dataset_clone(*ptr_signal_data);
          RooRealVar* scaled_weight = dynamic_cast<RooRealVar*>(signal_dataset_clone.addColumn(weight_function));
          RooDataSet signal_dataset("data", "data", RooArgSet(*wk_bkg->var("mass"), *scaled_weight), RooFit::Import(signal_dataset_clone), RooFit::WeightVar(*scaled_weight), RooFit::Cut(("0.9 * " + mX + " < mass && mass < 1.1 * " + mX).c_str()));
          double nInjectedEvents(signal_dataset.sumEntries());
          MSG_INFO("  ... injecting " << injected_signal_pb << " pb signal => corresponding to " << nInjectedEvents << " events");

          std::vector<double> values, errors;
          for (unsigned int iteration = 0; iteration < nIterations; ++iteration) {
            // Generate Asimov background and combine with signal
            RooDataSet* asimov_data = bkg_PDF->generate(*wk_bkg->var("mass"), RooFit::NumEvents(int(bkg_dataset->sumEntries())));
            RooDataSet combined_dataset(signal_dataset);
            combined_dataset.append(*asimov_data);

            // Construct combined PDF
            RooRealVar nSig("nSig", "number of signal events", combined_dataset.sumEntries(), 0, 5 * combined_dataset.sumEntries());
            RooRealVar nBkg("nBkg", "number of background events", asimov_data->numEntries(), 0, 1.5 * asimov_data->numEntries());
            RooAddPdf combined_PDF("signal_plus_" + bkg_PDF->getTitle(), "signal + " + bkg_PDF->getTitle(), RooArgList(*signal_model, *bkg_PDF), RooArgList(nSig, nBkg));

            // Fit S+B PDF to combined data
            combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
            combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1));
            MSG_ERROR(mX << " " << iteration << " " << injected_signal_pb << " " << nInjectedEvents << " " << nSig.getVal() << " +/- " << nSig.getError());
            values.push_back(nSig.getVal()); errors.push_back(nSig.getError());

            // Plot output
            if (injected_signal_pb == 5 && iteration == 0) {
              PlotStyle::EnsureAtlasStyle();
              RooPlot* frame = wk_bkg->var("mass")->frame();
              TLegend legend(0.5, 0.8, 0.93, 0.93);
              combined_dataset.plotOn(frame);
              legend.AddEntry((TGraph*)frame->getObject(frame->numItems() - 1), "Asimov bkg + injected signal", "P");
              combined_PDF.plotOn(frame);
              legend.AddEntry((TGraph*)frame->getObject(frame->numItems() - 1), ("S+B (" + PlotStyle::label(fitfn_name) + ") fit").c_str(), "L");
              TCanvas canvas("canvas", "canvas", 600, 600);
              frame->Draw();
              TLatex textBox; textBox.SetNDC(); textBox.SetTextFont(42); //textBox.SetTextSize(0.02);
              textBox.DrawLatex(0.5, 0.74, (mX + " GeV signal: " + PlotStyle::to_string(injected_signal_pb, 0) + " pb (hh)").c_str());
              legend.SetBorderSize(0);
              legend.SetFillStyle(0);
              legend.Draw();
              canvas.Print(("plots/" + mass_category + "Mass_" + tag_category + "tag/signal_bias_" + mass_category + "Mass_" + tag_category + "tag_mX_" + mX + "_" + PlotStyle::to_string(injected_signal_pb, 0) + "pb.pdf").c_str());
            }
          }

          // Combine weighted values
          double nExtracted(0.0), nExtractedError(0.0);
          for (unsigned int iteration = 0; iteration < nIterations; ++iteration) {
            if (errors.at(iteration) <= 0.0) { continue; }
            double reciprocal_sqrd_error(1.0 / (errors.at(iteration) * errors.at(iteration)));
            nExtracted += values.at(iteration) * reciprocal_sqrd_error;
            nExtractedError += reciprocal_sqrd_error;
          }
          nExtracted /= nExtractedError;
          nExtractedError = 1.0 / sqrt(nExtractedError);
          MSG_ERROR(" -> " << mX << " " << injected_signal_pb << " " << nInjectedEvents << " " << nExtracted << " +/- " << nExtractedError);

          // Write text output
          std::ofstream f_output;
          f_output.open(f_output_text, std::ios::app);
          f_output << mass_category << " " << tag_category << " " << mX << " " << injected_signal_pb << " " << nInjectedEvents << " " << nExtracted << " " << nExtractedError << std::endl;
          f_output.close();

          // // Write text output
          // std::ofstream f_text;
          // f_text.open(f_output_text, std::ios::app);
          // for (unsigned int idx = 0; idx < m_fit_functions.size(); ++idx) {
          //   double nSig(m_nSig.size() == 0 ? 0 : m_nSig.at(idx));
          //   double Z(m_nSig.size() == 0 ? 0 : nSig / m_nSigError.at(idx));
          //   double Z_uncertainty(m_nSig.size() == 0 ? 0 : m_nSigError_withSumW2.at(idx) / m_nSigError.at(idx));
          //   f_text << bkg_name(m_fit_functions.at(idx)) << " " << m_resonance_mass << " " << nSig << " " << Z << " " << Z_uncertainty << " " << m_chi2.at(idx) << " " << m_ndof.at(idx) << std::endl;
          // }
          // f_text.close();

        }
      }
    }
  }
  f_bkg.Close();
  return 0;
}