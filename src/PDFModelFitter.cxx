// Local
#include "PDFModelFitter.h"
#include "Logger.h"
#include "PlotStyle.h"
// STL
#include <fstream>
#include <iomanip>
// ROOT and RooFit
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
// #include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"

namespace SpuriousSignal {
  /**
   * PDFModelFitter constructor
   */
  PDFModelFitter::PDFModelFitter(RooDataSet& data, std::vector<RooAbsPdf*> fit_functions, const std::string& mass_category, const std::string& tag_category, const bool& verbose)
    : m_mass_category(mass_category)
    , m_tag_category(tag_category)
    , m_data(data)
    , m_fit_functions(fit_functions)
    , m_resonance_mass(-1)
    , m_verbose(verbose)
  {}

  void PDFModelFitter::fit()
  {
    RooFitResult* fit_result(0);
    m_nSig.clear();
    m_nSigError.clear();
    m_nSigError_withSumW2.clear();
    m_nBkg.clear();

    for (auto fit_fn : m_fit_functions) {
      // Get a pointer to nSig if this is part of the fit function
      RooAbsArg* nSig = fit_fn->getParameters(m_data)->find("nSig");
      if (nSig != 0) { dynamic_cast<RooRealVar*>(nSig)->setVal(0); }

      // Begin fitting with simplex + migrad + improve + hesse (+ minos)
      MSG_INFO("... fitting \033[1m" << m_mass_category << " mass " << m_tag_category << "-tag \033[0m with function \033[1m" << fit_fn->getTitle() << "\033[0m");
      // Simplex fit
      fit_result = fit_fn->fitTo(m_data, RooFit::SumW2Error(false), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Save(true));
      if (m_verbose && fit_result->status() != 0) {
        MSG_ERROR("... simplex fit did not converge: status = " << fit_result->status());
      }

      // Get signal uncertainty from Migrad + Improve + Hesse (with SumW2)
      if (nSig != 0) {
        fit_fn->fitTo(m_data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1));
        m_nSigError_withSumW2.push_back(dynamic_cast<RooRealVar*>(nSig)->getError());
      }

      // Migrad + Improve + Hesse
      fit_result = fit_fn->fitTo(m_data, RooFit::SumW2Error(false), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Save(true));
      if (fit_result->status() != 0) {
        MSG_ERROR("... migrad + improve + hesse fit did not converge: status = " << fit_result->status());
      }

      // Get signal value and error from Migrad + Improve + Hesse (without SumW2)
      if (nSig != 0) {
        m_nSig.push_back(dynamic_cast<RooRealVar*>(nSig)->getVal());
        m_nSigError.push_back(dynamic_cast<RooRealVar*>(nSig)->getError());
        m_nBkg.push_back(dynamic_cast<RooRealVar*>(fit_fn->getParameters(m_data)->find("nBkg"))->getVal());
      }

      // Print final fit result if requested
      if (m_verbose) { fit_result->Print("v"); }
    }
  }

  void PDFModelFitter::plot(RooPlot* frame, const int& resonance_mass)
  {
    PlotStyle::EnsureAtlasStyle();
    TCanvas canvas("canvas", "canvas", 600, 600);
    m_resonance_mass = resonance_mass;
    m_fit_graphs.clear();
    m_chi2.clear();
    m_ndof.clear();
    // Plot data and write to file
    m_data.plotOn(frame);
    m_fit_graphs["data"] = (TGraph*)frame->getObject(frame->numItems() - 1);
    // Add legend
    TLegend legend(0.4, 0.7, 0.93, 0.93);

    // Plot backgrounds and write to file
    for (unsigned int idx = 0; idx < m_fit_functions.size(); ++idx) {
      if (bkg_only()) {
        // Background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::LineColor(PlotStyle::colour(bkg_name(m_fit_functions.at(idx)))));
        m_fit_graphs[("bkg_" + m_tag_category + "tag" + bkg_name(m_fit_functions.at(idx)))] = (TGraph*)frame->getObject(frame->numItems() - 1);
      }
      else {
        // Restore nSig and nBkg
        dynamic_cast<RooRealVar*>(m_fit_functions.at(idx)->getParameters(m_data)->find("nSig"))->setVal(m_nSig.at(idx));
        dynamic_cast<RooRealVar*>(m_fit_functions.at(idx)->getParameters(m_data)->find("nBkg"))->setVal(m_nBkg.at(idx));
        // Background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::Components(bkg_name(m_fit_functions.at(idx)).c_str()), RooFit::LineColor(PlotStyle::colour(bkg_name(m_fit_functions.at(idx)))), RooFit::LineStyle(kDashed));
        m_fit_graphs["mX_" + std::to_string(resonance_mass) + "_bkg_" + m_tag_category + "tag" + bkg_name(m_fit_functions.at(idx))] = (TGraph*)frame->getObject(frame->numItems() - 1);
        // Signal + background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::LineColor(PlotStyle::colour(bkg_name(m_fit_functions.at(idx)))));
        m_fit_graphs["mX_" + std::to_string(resonance_mass) + "_splusb_" + m_tag_category + "tag" + bkg_name(m_fit_functions.at(idx))] = (TGraph*)frame->getObject(frame->numItems() - 1);
      }

      // Chi^2 of fit to data
      m_chi2.push_back(frame->chiSquare() * m_fit_graphs["data"]->GetN());
      m_ndof.push_back(m_fit_graphs["data"]->GetN() - m_fit_functions.at(idx)->getParameters(m_data)->getSize());
      MSG_INFO(std::setw(15) << m_fit_functions.at(idx)->getTitle() << " " << (bkg_only() ? "(bkg-only)" : "(S+B)") << ": chi2 / ndof =  " << m_chi2.back() << " / " << m_ndof.back());
      std::string s_chi2 = std::to_string(m_chi2.back());
      s_chi2 = s_chi2.substr(0, s_chi2.find(".") + 3);
      legend.AddEntry((TGraph*)frame->getObject(frame->numItems() - 1), (PlotStyle::label(bkg_name(m_fit_functions.at(idx))) + ": #chi^{2} / ndof = " + s_chi2 + " / " + std::to_string(m_ndof.back())).c_str(), "L");
    }
    legend.AddEntry(m_fit_graphs["data"], ("MC bkg: " + m_tag_category + "-tag, " + m_mass_category + " mass").c_str(), "P");
    frame->Draw();
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.Draw();

    if (bkg_only()) {
      canvas.Print(("plots/m_yyjj_" + m_mass_category + "Mass_" + m_tag_category + "tag_bkgOnly.pdf").c_str());
    }
    else {
      canvas.Print(("plots/" + m_mass_category + "Mass_" + m_tag_category + "tag/m_yyjj_" + m_mass_category + "Mass_" + m_tag_category + "tag_mX_" + std::to_string(resonance_mass) + ".pdf").c_str());
    }

    // MSG_INFO("Created \033[1m" + m_mass_category + " mass " + m_tag_category + "-tag\033[0m plot for " << (bkg_only() ? "background-only" : "signal-plus-background with mX = " + std::to_string(resonance_mass)));
    MSG_INFO("Created \033[1m" << m_mass_category << " mass " << m_tag_category << "-tag\033[0m plot for " << (bkg_only() ? "background-only" : "signal-plus-background with mX = " + std::to_string(resonance_mass)));
  }

  bool PDFModelFitter::bkg_only() const
  {
    return m_resonance_mass < 0;
  }

  std::string PDFModelFitter::bkg_name(RooAbsPdf* fit_function) const
  {
    std::string bkg_name(fit_function->getTitle());

    if (!this->bkg_only()) {
      bkg_name.erase(bkg_name.find("signal + "), sizeof("signal + ") - 1);
    }

    return bkg_name;
  }

  void PDFModelFitter::write(const std::string& f_output_ROOT, const std::string& f_output_text) const
  {
    // Write ROOT output
    TFile f_ROOT(f_output_ROOT.c_str(), "WRITE");

    for (auto graph_kv : m_fit_graphs) {
      f_ROOT.WriteObject(graph_kv.second, graph_kv.first.c_str());
    }

    f_ROOT.Close();
    // Write text output
    std::ofstream f_text;
    f_text.open(f_output_text, std::ios::app);

    for (unsigned int idx = 0; idx < m_fit_functions.size(); ++idx) {
      double nSig(m_nSig.size() == 0 ? 0 : m_nSig.at(idx));
      double Z(m_nSig.size() == 0 ? 0 : nSig / m_nSigError.at(idx));
      double Z_uncertainty(m_nSig.size() == 0 ? 0 : m_nSigError_withSumW2.at(idx) / m_nSigError.at(idx));
      f_text << bkg_name(m_fit_functions.at(idx)) << " " << m_resonance_mass << " " << nSig << " " << Z << " " << Z_uncertainty << " " << m_chi2.at(idx) << " " << m_ndof.at(idx) << std::endl;
    }

    f_text.close();
  }
}
