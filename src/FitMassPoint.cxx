#include "FitMassPoint.h"
#include "Logger.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include "RooChi2Var.h"

namespace SpuriousSignal {
  /**
   * FitMassPoint constructor
   */
  FitMassPoint::FitMassPoint(RooDataSet& data, std::vector<RooAbsPdf*> fit_functions, const std::string& mass_category)
    : m_data(data)
    , m_fit_functions(fit_functions)
    , m_mass_category(mass_category)
    , m_colours({kViolet, kGreen + 3, kBlue, kRed, kMagenta, kCyan})
  {}

  int FitMassPoint::fit(const bool& verbose) {
    // mass_resonance.setVal(iMass); mass_resonance.setConstant(true); int fit_status(0);
    RooFitResult* fit_result(0);
    for (auto fit_fn : m_fit_functions) {
      MSG_INFO("... fitting \033[1m" << m_mass_category << " mass\033[0m fit function \033[1m" << fit_fn->getTitle() << "\033[0m");
      // Simplex fit
      fit_result = fit_fn->fitTo(m_data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Save(true));
      if (fit_result->status() != 0) { MSG_ERROR("... simplex fit did not converge: status = " << fit_result->status()); }
      // Migrad + Improve + Hesse
      fit_result = fit_fn->fitTo(m_data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Save(true));
      if (fit_result->status() != 0) { MSG_ERROR("... migrad + improve + hesse fit did not converge: status = " << fit_result->status()); }
      if (verbose) {
        fit_result->Print("v");
      }
    }
    return fit_result->status();
  }

  void FitMassPoint::plot(RooPlot* frame, const int& resonance_mass, TFile& output_file) {
    bool bkg_only(resonance_mass < 0);
    TCanvas canvas("canvas", "canvas", 800, 600);
    // Plot data and write to file
    m_data.plotOn(frame);
    TGraph* g_data = (TGraph*)frame->getObject(frame->numItems() - 1);
    int ndof(g_data->GetN());
    output_file.WriteObject(g_data, "data");
    // Plot backgrounds and write to file
    for (unsigned int idx = 0; idx < m_fit_functions.size(); ++idx) {
      if (bkg_only) {
        // Background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::LineColor(m_colours.at(idx)));
        // Chi^2 of fit to data
        double chi2(frame->chiSquare() * ndof);
        MSG_INFO(std::setw(15) << m_fit_functions.at(idx)->getTitle() << " (bkg-only): chi2 / ndof =  " << chi2 << " / " << ndof - m_fit_functions.at(idx)->getParameters(m_data)->getSize());
        std::cout << "cf. " << frame->chiSquare(m_fit_functions.at(idx)->getParameters(m_data)->getSize()) << std::endl;
      } else {
        // Background
        std::string bkg_name(m_fit_functions.at(idx)->getTitle());
        bkg_name.erase(bkg_name.find("signal + "), sizeof("signal + ") - 1);
        m_fit_functions.at(idx)->plotOn(frame, RooFit::Components(bkg_name.c_str()), RooFit::LineColor(m_colours.at(idx)), RooFit::LineStyle(kDashed));
        output_file.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "bkg_" + m_fit_functions.at(idx)->getTitle());
        // Signal + background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::LineColor(m_colours.at(idx)));
        output_file.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "splusb_" + m_fit_functions.at(idx)->getTitle());
        // Chi^2 of fit to data
        double chi2(frame->chiSquare() * ndof);
        MSG_INFO(std::setw(15) << m_fit_functions.at(idx)->getTitle() << " (S+B): chi2 / ndof =  " << chi2 << " / " << ndof - m_fit_functions.at(idx)->getParameters(m_data)->getSize());
      }
    }
    frame->Draw();
    // Add legend
    TLegend legend(0.6, 0.7, 0.9, 0.9);
    if (bkg_only) {
      legend.AddEntry(frame->findObject("novosibirsk_Norm[mass]"), "Novosibirsk", "L");
      // legend.AddEntry(frame->findObject("gammagaus_Norm[mass]"), "Gamma + Gaussian", "L");
      legend.AddEntry(frame->findObject("invgaus_Norm[mass]"), "Inverse Gaussian", "L");
      legend.AddEntry(frame->findObject("modified_landau_Norm[mass]"), "Modified Landau", "L");
      legend.AddEntry(frame->findObject("modified_gamma_Norm[mass]"), "Modified Gamma", "L");
      // legend.AddEntry(frame->findObject("modified_cauchy_Norm[mass]"), "Modified Cauchy", "L");
      // legend.AddEntry(frame->findObject("log_normal_Norm[mass]"), "Log Normal", "L");
      // legend.AddEntry(frame->findObject("rayleigh_Norm[mass]"), "Rayleigh", "L");
      // legend.AddEntry(frame->findObject("inverted_argus_Norm[mass]"), "Inverted Argus", "L");
      // legend.AddEntry(frame->findObject("poisson_Norm[mass]"), "Poisson", "L");
      // legend.AddEntry(frame->findObject("gaussexp_Norm[mass]"), "Gaussian * exponential", "L");
    } else {
      legend.AddEntry(frame->findObject("signal_plus_novosibirsk_Norm[mass]"), "Novosibirsk", "L");
      // legend.AddEntry(frame->findObject("signal_plus_gammagaus_Norm[mass]"), "Gamma + Gaussian", "L");
      legend.AddEntry(frame->findObject("signal_plus_invgaus_Norm[mass]"), "Inverse Gaussian", "L");
      legend.AddEntry(frame->findObject("signal_plus_modified_landau_Norm[mass]"), "Modified Landau", "L");
      legend.AddEntry(frame->findObject("signal_plus_modified_gamma_Norm[mass]"), "Modified Gamma", "L");
      // legend.AddEntry(frame->findObject("signal_plus_modified_cauchy_Norm[mass]"), "Modified Cauchy", "L");
      // legend.AddEntry(frame->findObject("signal_plus_log_normal_Norm[mass]"), "Log Normal", "L");
      // legend.AddEntry(frame->findObject("signal_plus_rayleigh_Norm[mass]"), "Rayleigh", "L");
      // legend.AddEntry(frame->findObject("signal_plus_inverted_argus_Norm[mass]"), "Inverted Argus", "L");
      // legend.AddEntry(frame->findObject("signal_plus_poisson_Norm[mass]"), "Poisson", "L");
      // legend.AddEntry(frame->findObject("signal_plus_gaussexp_Norm[mass]"), "Gaussian * exponential", "L");
    }
    legend.Draw();
    if (bkg_only) {
      canvas.Print((std::string("output/m_yyjj_") + m_mass_category + "Mass_bkgOnly.pdf").c_str());
    } else {
      canvas.Print((std::string("output/m_yyjj_") + m_mass_category + "Mass_mX_" + std::to_string(resonance_mass) + ".pdf").c_str());
    }
    MSG_INFO("Saved " + m_mass_category + " mass plot");
  }
}