// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "ToyGenerator.h"
// STL
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
// ROOT and RooFit/RooStats
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooConstVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooMsgService.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#define PI (3.14159265358979323846264338327950288)


int main(int argc, char** argv)
{
  using namespace SpuriousSignal;

  // Determine whether to use toys and/or Asimov
  bool useToys(true), useAsimov(true); int nToys(1000);
  std::vector<std::string> input_args(argv, argv + argc);
  for (auto itArg = input_args.begin(); itArg != input_args.end(); ++itArg) {
    if (*itArg == "--disable-toys") { useToys = false; }
    if (*itArg == "--disable-asimov") { useAsimov = false; }
    if (*itArg == "--nToys" && (itArg + 1) != input_args.end()) { nToys = std::stoi(*(itArg + 1)); }
  }
  if (nToys < 1) { useToys = false; }

  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gErrorIgnoreLevel = kBreak;

  // Define signal injection points
  // std::vector<int> masses = {260};//, 275, 300, 325, 350, 400, 450, 500, 750, 1000};
  // std::map< std::string, std::vector<double> > injected_signals = {{"ATLAS_Run_1_combination", {1.7, 1.7, 1.6, 1.6, 1.5, 0.8, 0.7, 0.6, 0.35, 0.01}},
  //                                                                  {"ATLAS_Run_2_bbyy", {7, 7, 6.1, 5.6, 5.1, 4.1, 4, 4, 4, 4}},
  //                                                                  {"ATLAS_Run_2_bbbb", {99, 99, 99, 99, 99, 1.48, 0.89, 0.73, 0.16, 0.11}},
  //                                                                  {"CMS_Run_2_bbyy", {1.23, 1.40, 1.04, 0.60, 0.80, 0.42, 0.48, 0.23, 0.10, 0.31}},
  //                                                                  {"Expected_bbyy", {1.15, 1.0, 0.9, 0.8, 0.7, 0.55, 0.5, 0.38, 0.18, 0.13}},
  //                                                                  {"Small_xs", {0.00115, 0.0010, 0.0009, 0.0008, 0.0007, 0.00055, 0.0005, 0.00038, 0.00018, 0.00013}}};
  std::map< std::string, std::vector<double> > injected_signals = {{"Expected_bbyy", {1.15, 1.0, 0.9, 0.8, 0.7, 0.55, 0.5, 0.38, 0.18, 0.13}}};
  std::map< std::string, std::vector<std::string> > bkg_model_parameters = {{"novosibirsk", {"novosibirsk_width", "novosibirsk_peak", "novosibirsk_tail"}},
                                                                            {"exppoly1", {"exppoly1_p1"}}};

  // Construct mass and tag categories
  std::vector<std::string> mass_categories({"low", "high"});
  std::vector<std::string> tag_categories({"1", "2"});

  // Open workspace files
  TFile f_bkg("output/background_model_workspace.root", "READ");
  TFile f_signal("output/signal_model_workspace.root", "READ");

  // Recreate output text file
  std::string f_output_text_base("output/csv/signal_bias.csv");
  std::string f_output_text(f_output_text_base + ".inprogress");
  std::ofstream f_output;
  f_output.open(f_output_text, std::ios::trunc);
  f_output.close();

  // Define weight parameter
  RooRealVar weight("weight", "event weight", -1e10, 1e10);

  // Iterate over mass and tag categories
  for (auto mass_category : mass_categories) {
    std::string fitfn_name(mass_category == "low" ? "novosibirsk" : "exppoly1");
    std::vector<int> masses = PlotStyle::resonance_masses(mass_category);

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
      double nBkgEvents(bkg_dataset->sumEntries());
      MSG_INFO("Loaded background dataset with " << bkg_dataset->numEntries() << " entries, corresponding to " << nBkgEvents << " events");
      wk_bkg->factory(("ExtendPdf::bkg_model(" + fitfn_name + ", norm[" + std::to_string(nBkgEvents) + "])").c_str());
      RooDataHist* asimov_bkg = wk_bkg->pdf("bkg_model")->generateBinned(*wk_bkg->var("mass"), RooFit::ExpectedData());

      // Take a snapshot of the bkg parameters
      RooArgSet* bkg_params = wk_bkg->pdf(fitfn_name.c_str())->getParameters(*wk_bkg->var("mass"));
      wk_bkg->saveSnapshot("bkg_parameters", *bkg_params);

      // Novosibirsk peak is constrained to be within 10% of the best fit value using the data sidebands
      if (fitfn_name == "novosibirsk") {
        wk_bkg->var("novosibirsk_peak")->setRange(wk_bkg->var("novosibirsk_peak")->getVal() - (tag_category == 1 ? 6.55 : 4.99), wk_bkg->var("novosibirsk_peak")->getVal() + (tag_category == 1 ? 6.55 : 4.99));
        wk_bkg->var("novosibirsk_width")->setRange(wk_bkg->var("novosibirsk_width")->getVal() - (tag_category == 1 ? 6.03 : 4.79), wk_bkg->var("novosibirsk_width")->getVal() + (tag_category == 1 ? 6.03 : 4.79));
        wk_bkg->var("novosibirsk_tail")->setRange(wk_bkg->var("novosibirsk_tail")->getVal() - (tag_category == 1 ? 0.39 : 2.34), wk_bkg->var("novosibirsk_tail")->getVal() + (tag_category == 1 ? 0.39 : 2.34));
      }

      // Iterate over mass points
      for (unsigned int idxMass = 0; idxMass < masses.size(); ++idxMass) {
        // Set mass point
        int mX(masses.at(idxMass));
        std::string s_mX(std::to_string(mX));

        // Load signal data
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
        RooDataSet* ptr_signal_data = RooDataSet::read(("input/m_yyjj_Xhh_m" + s_mX + "_" + mass_category + "Mass_" + tag_category + "tag_tightIsolated.csv").c_str(), RooArgList(*wk_bkg->var("mass"), weight));
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        MSG_INFO("... loading \033[1m" << mX << " GeV\033[0m dataset containing " << ptr_signal_data->numEntries() << " entries");

        // Load signal PDF
        RooAbsPdf* signal_model = wk_sig->pdf("signal_PDF");
        wk_sig->var("mass_resonance")->setVal(mX);
        wk_sig->var("mass_resonance")->setConstant(true);
        wk_sig->var("mass")->setRange(wk_bkg->var("mass")->getMin(), wk_bkg->var("mass")->getMax());

        for (auto& kv : injected_signals) {
          if (idxMass >= kv.second.size()) {
            MSG_WARNING("No cross-section is defined for \033[1m" << kv.first << "\033[0m at this mass point");
            continue;
          }
          double injected_signal_pb(kv.second.at(idxMass));
          MSG_INFO("  ... using limits from \033[1m" << kv.first << "\033[0m");

          // Construct formula to calculate scaled weight for events
          RooFormulaVar weight_function("scaled_weight", "scaled_weight", ("weight * " + std::to_string(injected_signal_pb) + " / 5.0").c_str(), weight);
          RooDataSet signal_dataset_clone(*ptr_signal_data);
          RooRealVar* scaled_weight = dynamic_cast<RooRealVar*>(signal_dataset_clone.addColumn(weight_function));
          RooDataSet signal_dataset("data", "data", RooArgSet(*wk_bkg->var("mass"), *scaled_weight), RooFit::Import(signal_dataset_clone), RooFit::WeightVar(*scaled_weight), RooFit::Cut(("0.9 * " + s_mX + " < mass && mass < 1.1 * " + s_mX).c_str()));
          double nInjectedEvents(signal_dataset.sumEntries());
          MSG_INFO("  ... injecting " << injected_signal_pb << " pb signal => corresponding to " << nInjectedEvents << " events");

          // Construct combined PDF
          RooRealVar nSig("nSig", "number of signal events", nInjectedEvents, 0, 10 * nInjectedEvents);
          RooRealVar nBkg("nBkg", "number of background events", nBkgEvents, 0, 10 * nBkgEvents);
          RooAddPdf combined_PDF(("signal_plus_" + fitfn_name).c_str(), ("signal_plus_" + fitfn_name).c_str(), RooArgList(*signal_model, *wk_bkg->pdf(fitfn_name.c_str())), RooArgList(nSig, nBkg));

          // Run with toys if specified
          double nSig_toys(-1), nSig_toys_err(-1); int n_median(nToys / 2);
          if (useToys) {
            MSG_INFO("  ... calculating limits with \033[1m" << nToys << "\033[0m toys");
            std::vector<double> nSigFit(0), nSigInj(0), _deviations(nToys);

            for (int idxToy = 0; idxToy < nToys; ++idxToy) {
              // Reset background parameters
              wk_bkg->loadSnapshot("bkg_parameters");

              // Construct signal and background toys
              RooDataSet* toy_signal = ToyGenerator::instance()->generate(&signal_dataset);
              RooDataSet* toy_bkg = wk_bkg->pdf("bkg_model")->generate(RooArgSet(*wk_bkg->var("mass")), RooFit::NumEvents(nBkgEvents), RooFit::Extended());
              RooDataSet toy_combined(*toy_signal);
              toy_combined.append(*toy_bkg);

              // Set sensible initial parameters
              nBkg.setVal(toy_bkg->numEntries());
              nSig.setVal(toy_signal->numEntries());

              // Constrain number of background events (would be possible from sidebands in reality)
              RooGaussian nBkg_constraint("nBkg_constraint", "nBkg_constraint", nBkg, RooFit::RooConst(toy_bkg->numEntries()), RooFit::RooConst(0.1 * toy_bkg->numEntries()));

              // Fit with combined PDF
              combined_PDF.fitTo(toy_combined, RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::ExternalConstraints(nBkg_constraint));
              combined_PDF.fitTo(toy_combined, RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::ExternalConstraints(nBkg_constraint));

              // Store final values
              nSigFit.push_back(nSig.getVal());
              nSigInj.push_back(toy_signal->numEntries());
            }
            // Median
            std::nth_element(nSigFit.begin(), nSigFit.begin() + n_median, nSigFit.end());
            nSig_toys = nSigFit.at(n_median);
            // Standard deviation
            std::transform(nSigFit.begin(), nSigFit.end(), _deviations.begin(), [nSig_toys](double x) { return x - nSig_toys; });
            double _std_dev = std::sqrt(std::inner_product(_deviations.begin(), _deviations.end(), _deviations.begin(), 0.0) / (_deviations.size() - 1));
            // Error on the median https://physics.stackexchange.com/questions/292871/uncertainty-propagation-upon-taking-the-median
            nSig_toys_err = _std_dev * std::sqrt(PI * (2 * nToys + 1) / (4 * nToys * nToys)) ;

            // Plot output
            if (kv.first == "Expected_bbyy") {
              PlotStyle::EnsureAtlasStyle();
              int _range_max(4 * nInjectedEvents);
              TH1D toys("toys", "toys", _range_max, 0, _range_max);
              TH1D injected("injected", "injected", _range_max, 0, _range_max);
              for (auto nSig : nSigFit ) { toys.Fill(nSig); }
              for (auto nInj : nSigInj ) { injected.Fill(nInj); }
              TCanvas canvas("canvas", "canvas", 600, 600);
              TLegend legend(0.5, 0.8, 0.93, 0.93);
              TLatex textBox; textBox.SetNDC(); textBox.SetTextFont(42);
              textBox.DrawLatex(0.5, 0.74, (s_mX + " GeV signal").c_str());
              toys.Draw("e hist");
              injected.Draw("hist same");
              injected.SetLineColor(kRed);
              toys.GetXaxis()->SetTitle("Number of signal events");
              toys.GetYaxis()->SetTitle("Number of toys");
              toys.SetMinimum(0);
              canvas.Update();
              TLine l_fit_signal(nSig_toys, canvas.GetUymin(), nSig_toys, canvas.GetUymax());
              l_fit_signal.SetLineColor(kBlack); l_fit_signal.SetLineStyle(kDashed); l_fit_signal.Draw("same");
              TLine l_inj_signal(nInjectedEvents, canvas.GetUymin(), nInjectedEvents, canvas.GetUymax());
              l_inj_signal.SetLineColor(kRed); l_inj_signal.SetLineStyle(kDashed); l_inj_signal.Draw("same");
              legend.AddEntry(&toys, ("Fit result: (" + PlotStyle::to_string(nToys, 0) + " toys)").c_str(), "P");
              legend.AddEntry(&l_fit_signal, "Median of toys", "L");
              legend.AddEntry(&injected, "Injected signal", "L");
              legend.AddEntry(&l_inj_signal, "Nominal signal", "L");
              legend.SetBorderSize(0);
              legend.SetFillStyle(0);
              legend.Draw();
              canvas.Print(("plots/signal_bias/toy_fits/signal_bias_" + mass_category + "Mass_" + tag_category + "tag_mX_" + s_mX + "_toys.pdf").c_str());
            }
          }

          // Run with Asimov if specified
          if (useAsimov) {
            MSG_INFO("  ... calculating limits with \033[1mAsimov background\033[0m");

            // Combine Asimov background with signal
            RooDataHist combined_dataset(*asimov_bkg);
            combined_dataset.add(signal_dataset);

            // Set sensible initial parameters
            nBkg.setVal(nBkgEvents);
            nSig.setVal(nInjectedEvents);
            wk_bkg->loadSnapshot("bkg_parameters");

            // Constrain number of background events (would be possible from sidebands in reality)
            RooGaussian nBkg_constraint("nBkg_constraint", "nBkg_constraint", nBkg, RooFit::RooConst(nBkgEvents), RooFit::RooConst(0.1 * nBkgEvents));

            // Fit S+B PDF to combined data
            combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::ExternalConstraints(nBkg_constraint));
            // combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(false), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(false), RooFit::Minos(true), RooFit::PrintLevel(-1), RooFit::ExternalConstraints(nBkg_constraint));
            // std::cout << "- sumW2Error: " << nSig.getVal() << " from " << nInjectedEvents << ", " << nBkg.getVal() << " from " << nBkgEvents << std::endl;
            combined_PDF.fitTo(combined_dataset, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(false), RooFit::Minos(true), RooFit::PrintLevel(-1), RooFit::ExternalConstraints(nBkg_constraint));
            // std::cout << "+ sumW2Error: " << nSig.getVal() << " from " << nInjectedEvents << ", " << nBkg.getVal() << " from " << nBkgEvents << std::endl;

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
              textBox.DrawLatex(0.5, 0.74, (s_mX + " GeV signal: " + PlotStyle::to_string(injected_signal_pb, 2) + " pb (hh)").c_str());
              legend.SetBorderSize(0);
              legend.SetFillStyle(0);
              legend.Draw();
              canvas.Print(("plots/signal_bias/asimov_fits/signal_bias_" + mass_category + "Mass_" + tag_category + "tag_mX_" + s_mX + "_asimov.pdf").c_str());
            }
          }

          // Write text output
          std::ofstream f_output;
          f_output.open(f_output_text, std::ios::app);
          f_output << mass_category << " " << tag_category << " " << s_mX << " " << kv.first << " " << nInjectedEvents << " " << nSig.getVal() << " " << nSig.getError() << " " << nSig_toys << " " << nSig_toys_err << " " << nToys << std::endl;
          f_output.close();
          MSG_INFO(mass_category << " " << tag_category << " " << s_mX << " " << kv.first << " " << nInjectedEvents << ": asimov = " << nSig.getVal() << " +/- " << nSig.getError() << ", toys = " << nSig_toys << " +/- " << nSig_toys_err);
        }
      }
    }
  }
  // Close ROOT files
  f_signal.Close();
  f_bkg.Close();

  // Replace any previous csv file
  remove(f_output_text_base.c_str());
  rename(f_output_text.c_str(), f_output_text_base.c_str());

  return 0;
}