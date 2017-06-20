#include "FitMassPoint.h"
#include "Logger.h"
#include "RooChi2Var.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace SpuriousSignal {
  /**
   * FitMassPoint constructor
   */
  FitMassPoint::FitMassPoint(RooDataSet& data, std::vector<RooAbsPdf*> fit_functions, const std::string& mass_category, const bool& verbose)
    : m_mass_category(mass_category)
    , m_data(data)
    , m_fit_functions(fit_functions)
    , m_resonance_mass(-1)
    , m_colours({kViolet, kGreen + 3, kBlue, kRed, kMagenta, kCyan})
    , m_fn_names({{"novosibirsk", "Novosibirsk"}, {"modified_gamma", "Modified Gamma"}, {"modified_landau", "Modified Landau"}})
    , m_verbose(verbose)
  {}

  void FitMassPoint::fit() {
    RooFitResult* fit_result(0);
    m_nSig.clear(); m_nSigError.clear(); m_nSigError_withSumW2.clear(); 
    for (auto fit_fn : m_fit_functions) {
      // Get a pointer to nSig if this is part of the fit function
      RooAbsArg* nSig = fit_fn->getParameters(m_data)->find("nSig");
      MSG_INFO("... fitting \033[1m" << m_mass_category << " mass\033[0m with function \033[1m" << fit_fn->getTitle() << "\033[0m");
      std::cout << "there are " << m_nSig.size() << ", " << m_nSigError.size()<< " and " << m_nSigError_withSumW2.size() << " signal numbers" << std::endl;
      // Simplex fit
      fit_result = fit_fn->fitTo(m_data, RooFit::SumW2Error(false), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Save(true));
      if (m_verbose && fit_result->status() != 0) { MSG_ERROR("... simplex fit did not converge: status = " << fit_result->status()); }
      // Get signal uncertainty from Migrad + Improve + Hesse (with SumW2)
      if (nSig != 0) {
        fit_fn->fitTo(m_data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1));
        m_nSigError_withSumW2.push_back(dynamic_cast<RooRealVar*>(nSig)->getError());
      }
      // Migrad + Improve + Hesse
      fit_result = fit_fn->fitTo(m_data, RooFit::SumW2Error(false), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Save(true));
      if (fit_result->status() != 0) { MSG_ERROR("... migrad + improve + hesse fit did not converge: status = " << fit_result->status()); }
      // Get signal value and error from Migrad + Improve + Hesse (without SumW2)
      if (nSig != 0) {
        m_nSig.push_back(dynamic_cast<RooRealVar*>(nSig)->getVal());
        m_nSigError.push_back(dynamic_cast<RooRealVar*>(nSig)->getError());
      }
      std::cout << "-> now there are " << m_nSig.size() << ", " << m_nSigError.size() << " and " << m_nSigError_withSumW2.size() << " signal numbers" << std::endl;
      // Print final fit result if requested
      if (m_verbose) { fit_result->Print("v");  }
    }
  }

  // TGraphErrors* nSignal = fitValues("nSignal", pdf, false, pdf.GetName(), "Spurious Signal");
  // TGraphErrors* nSignal_sumW2Errors = fitValues("nSignal", pdf, true, pdf.GetName(), "Spurious Signal");
  // TGraphErrors* signalZ = Z(*nSignal, *nSignal_sumW2Errors, pdf.GetName() + "_Z", "S/#deltaS");
  // TGraphErrors* Selector::Z(const TGraphErrors& graph, const TGraphErrors& uncert, const TString& name, const TString& yTitle)
  // zGraph->SetPoint(i, graph.GetX()[i], graph.GetY()[i]/graph.GetErrorY(i));
  // zGraph->SetPointError(i, 0, uncert.GetErrorY(i)/graph.GetErrorY(i));
  // TGraphErrors* Selector::fitValues(const TString& varName, const Parameterization& pdf, bool sumW2Errors, const TString& name, const TString& yTitle)

  void FitMassPoint::plot(RooPlot* frame, const int& resonance_mass) {
    TCanvas canvas("canvas", "canvas", 800, 600);
    m_resonance_mass = resonance_mass;
    m_fit_graphs.clear();
    m_chi2.clear(); m_ndof.clear(); 
    // Plot data and write to file
    m_data.plotOn(frame);
    m_fit_graphs["data"] = (TGraph*)frame->getObject(frame->numItems() - 1);
    // Add legend
    TLegend legend(0.4, 0.6, 0.9, 0.9);
    // Plot backgrounds and write to file
    for (unsigned int idx = 0; idx < m_fit_functions.size(); ++idx) {
      if (bkg_only()) {
        // Background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::LineColor(m_colours.at(idx)));
        m_fit_graphs[("bkg_" + bkg_name(m_fit_functions.at(idx)))] = (TGraph*)frame->getObject(frame->numItems() - 1);
      } else {
        // Background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::Components(bkg_name(m_fit_functions.at(idx)).c_str()), RooFit::LineColor(m_colours.at(idx)), RooFit::LineStyle(kDashed));
        m_fit_graphs["mX_" + std::to_string(resonance_mass) + "_bkg_" + bkg_name(m_fit_functions.at(idx))] = (TGraph*)frame->getObject(frame->numItems() - 1);
        // Signal + background
        m_fit_functions.at(idx)->plotOn(frame, RooFit::LineColor(m_colours.at(idx)));
        m_fit_graphs["mX_" + std::to_string(resonance_mass) + "_splusb_" + bkg_name(m_fit_functions.at(idx))] = (TGraph*)frame->getObject(frame->numItems() - 1);
      }
      // Chi^2 of fit to data
      m_chi2.push_back(frame->chiSquare() * m_fit_graphs["data"]->GetN());
      m_ndof.push_back(m_fit_graphs["data"]->GetN() - m_fit_functions.at(idx)->getParameters(m_data)->getSize());
      MSG_INFO(std::setw(15) << m_fit_functions.at(idx)->getTitle() << " " << (bkg_only() ? "(bkg-only)" : "(S+B)" ) << ": chi2 / ndof =  " << m_chi2.back() << " / " << m_ndof.back());
      std::string s_chi2 = std::to_string(m_chi2.back()); s_chi2 = s_chi2.substr(0, s_chi2.find(".") + 3);
      legend.AddEntry((TGraph*)frame->getObject(frame->numItems() - 1), (m_fn_names[bkg_name(m_fit_functions.at(idx))] + ": #chi^{2} / ndof = " + s_chi2 + " / " + std::to_string(m_ndof.back())).c_str() , "L");
    }
    frame->Draw();
    legend.Draw();
    if (bkg_only()) {
      canvas.Print((std::string("output/m_yyjj_") + m_mass_category + "Mass_bkgOnly.pdf").c_str());
    } else {
      canvas.Print((std::string("output/m_yyjj_") + m_mass_category + "Mass_mX_" + std::to_string(resonance_mass) + ".pdf").c_str());
    }
    MSG_INFO("Created " + m_mass_category + " mass plot for " << (bkg_only() ? "background-only" : "signal-plus-background with mX = " + std::to_string(resonance_mass)));
  }

  bool FitMassPoint::bkg_only() const {
    return m_resonance_mass < 0;
  }

  std::string FitMassPoint::bkg_name(RooAbsPdf* fit_function) const {
    std::string bkg_name(fit_function->getTitle());
    if (!this->bkg_only()) {
      bkg_name.erase(bkg_name.find("signal + "), sizeof("signal + ") - 1);
    }
    return bkg_name;
  }

  void FitMassPoint::write(const std::string& f_output_ROOT, const std::string& f_output_text) const {
    // Write ROOT output
    TFile f_ROOT(f_output_ROOT.c_str(), "WRITE");
    for (auto graph_kv : m_fit_graphs) {
      f_ROOT.WriteObject(graph_kv.second, graph_kv.first.c_str());
    }
    f_ROOT.Close();

    // Calculate Z
      // Z.setError(nSig.getError());
      // Z.setVal(nSig.getVal() / nSig.getError());
      // Z.setError(Z.getError() / nSig.getError());
      // // std::cout << "Z: " << Z.getVal() << " +/- " << Z.getError() << std::endl;


    // Write text output
    std::ofstream f_text;
    f_text.open(f_output_text, std::ios::app);
    for (unsigned int idx = 0; idx < m_fit_functions.size(); ++idx) {
      double nSig(m_nSig.size() == 0 ? 0 : m_nSig.at(idx));
      double Z(m_nSig.size() == 0 ? 0 : nSig / m_nSigError.at(idx));
      double Z_uncertainty(m_nSig.size() == 0 ? 0 : m_nSigError_withSumW2.at(idx) / m_nSigError.at(idx));
      std::cout << bkg_name(m_fit_functions.at(idx)) << " " << m_resonance_mass << " " << nSig << " " << Z << " " << Z_uncertainty << " " << m_chi2.at(idx) << " " << m_ndof.at(idx) << std::endl;
      f_text << bkg_name(m_fit_functions.at(idx)) << " " << m_resonance_mass << " " << nSig << " " << Z << " " << Z_uncertainty << " " << m_chi2.at(idx) << " " << m_ndof.at(idx) << std::endl;
    }
    f_text.close();
  }
}
