// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "SignalModel.h"
#include "RooDSCB.h"
#include "RooExpGausExp.h"
// ROOT and RooFit
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "TAxis.h"
#include "TCanvas.h"
// #include "TColor.h"
#include <algorithm>
#include "RooHist.h"
#include "TPad.h"
#include "TLine.h"
// #include "TPaveText.h"
#include "TLatex.h"

namespace SpuriousSignal {
  /**
   * SignalModel constructor
   */
  SignalModel::SignalModel(const std::string& mass_category, const std::string& tag_category)
    : m_mass_category(mass_category)
    , m_tag_category(tag_category)
    , m_wk(0)
    , m_data(0)
  {
    m_wk = new RooWorkspace(("signal_model_" + mass_category + "Mass_" + tag_category + "tag").c_str());
    // Crystal Ball
    m_wk->factory("CB_alpha_p0[0.40, 0, 5]");
    m_wk->factory("CB_alpha_p1[-0.003, -0.01, 0.01]");
    m_wk->factory("CB_frac_p0[-2.9, -50, 50]");
    m_wk->factory("CB_frac_p1[-0.003, -1, 1]");
    m_wk->factory("CB_mu_p0[4.7, -20, 20]");
    m_wk->factory("CB_mu_p1[0.001, -0.1, 0.1]");
    m_wk->factory("CB_n_p0[10, 0, 20]");
    m_wk->factory("CB_n_p1[-0.0015, -1, 1]");
    m_wk->factory("CB_sigma_p0[-8.7, -20, 20]");
    m_wk->factory("CB_sigma_p1[0.06, -1, 1]");
    m_wk->factory("CB_sigma_p2[-8e-5, -0.001, 0.001]");
    // Gaussian
    m_wk->factory("gaus_mu_p0[4.7, 0, 10]");
    m_wk->factory("gaus_mu_p1[0.003, -0.01, 0.01]");
    m_wk->factory("gaus_sigma_k_p0[8.7, -20, 20]");
    m_wk->factory("gaus_sigma_k_p1[-0.05, -1, 1]");
    m_wk->factory("gaus_sigma_k_p2[7e-5, -0.001, 0.001]");
    // Set up category to distinguish samples
    RooCategory mass_points("mass_points", "mass_points");
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) { mass_points.defineType(std::to_string(resonance_mass).c_str()); }
    m_wk->import(mass_points);
    m_wk->addClassDeclImportDir("include/");
    // m_wk->importClassCode("RooDSCB*");
  }

  void SignalModel::build_simultaneous_PDF(RooRealVar& mass)
  {
    m_wk->import(mass);
    m_wk->factory("Simultaneous::signal_model(mass_points)");
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category) ) { add_mass_point(resonance_mass); }
    // Also build a simple CB+G shape
    m_wk->factory("single_CB_mu[260.325, 0, +INF]");
    m_wk->factory("single_CB_sigma[6.5, 0, 50]");
    m_wk->factory("single_CB_alpha[-0.1, -5, 5]");
    m_wk->factory("single_CB_n[2, 0, 200]");
    m_wk->factory("single_gaus_mu[0, -10, 10]");
    m_wk->factory("single_gaus_sigma_k[6, 2, 10]");
    m_wk->factory("single_CB_frac[0.9, 0.0, 1.0]");
    m_wk->factory("prod::single_gaus_sigma(single_gaus_sigma_k, single_CB_sigma)");
    // m_wk->factory("CBShape::single_signal_PDF(mass, single_CB_mu, single_CB_sigma, single_CB_alpha, single_CB_n)");
    m_wk->factory("CBShape::single_CB(mass, single_CB_mu, single_CB_sigma, single_CB_alpha, single_CB_n)");
    // m_wk->factory("Gaussian::single_gaus(mass, single_gaus_mu, single_gaus_sigma)");
    // m_wk->factory("Gaussian::single_gaus(mass, single_CB_mu, single_gaus_sigma)");
    // m_wk->factory("SUM::single_signal_PDF(single_CB_frac * single_CB, single_gaus)");
    m_wk->factory("single_DSCB_mu[325.20, 0, +INF]");
    m_wk->factory("single_DSCB_sigma[1.7, 0.01, 8.0]");
    m_wk->factory("single_DSCB_alphaLo[0.2, 0.15, 1.0]");
    m_wk->factory("single_DSCB_nLo[9.0, 0.1, 20]");
    m_wk->factory("single_DSCB_alphaHi[0.2, 0.15, 1.0]");
    m_wk->factory("single_DSCB_nHi[5.0, 0.1, 10]");
    // m_wk->factory("RooDSCB::single_signal_PDF(mass, single_DSCB_mu, single_DSCB_sigma, single_DSCB_alphaLo, single_DSCB_nLo, single_DSCB_alphaHi, single_DSCB_nHi)");
    m_wk->factory("single_EGE_mu[260, 0, +INF]");
    m_wk->factory("single_EGE_sigma[4, 0.0, 10]");
    m_wk->factory("single_EGE_kLo[2, 0.1, +INF]");
    m_wk->factory("single_EGE_kHi[2, 0.1, +INF]");
    m_wk->factory("RooExpGausExp::single_signal_PDF(mass, single_EGE_mu, single_EGE_sigma, single_EGE_kLo, single_EGE_kHi)");
    m_wk->importClassCode("RooExpGausExp*");
    // m_wk->factory("Landau::conv(mass, ml[260, 0, +INF], sl[20, 0, 100])");
    // m_wk->factory("Gaussian::gauss(mass, mg[0], sg[1e-10])");
    // m_wk->factory("FCONV::conv(mass, single_CB, gauss)");
    m_wk->Print("v");
  }

  void SignalModel::add_mass_point(const int& resonance_mass) {
    std::string mX(std::to_string(resonance_mass));
    m_wk->factory(("expr::CB_mu_Xhh_m" + mX + "('TMath::Exp(CB_mu_p0 + CB_mu_p1 * " + mX + ")', CB_mu_p0, CB_mu_p1)").c_str());
    m_wk->factory(("expr::CB_sigma_Xhh_m" + mX + "('TMath::Exp(CB_sigma_p0 + CB_sigma_p1 * " + mX + " + CB_sigma_p2 * " + mX + " * " + mX + ")', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2)").c_str());
    m_wk->factory(("expr::CB_alpha_Xhh_m" + mX + "('TMath::Exp(CB_alpha_p0 + CB_alpha_p1 * " + mX + ")', CB_alpha_p0, CB_alpha_p1)").c_str());
    m_wk->factory(("expr::CB_n_Xhh_m" + mX + "('TMath::Exp(CB_n_p0 + CB_n_p1 * " + mX + ")', CB_n_p0, CB_n_p1)").c_str());
    m_wk->factory(("expr::CB_frac_Xhh_m" + mX + "('0.75 + 0.25 * TMath::TanH(CB_frac_p0 + CB_frac_p1 * " + mX + ")', CB_frac_p0, CB_frac_p1)").c_str());
    m_wk->factory(("expr::gaus_mu_Xhh_m" + mX + "('TMath::Exp(gaus_mu_p0 + gaus_mu_p1 * " + mX + ")', gaus_mu_p0, gaus_mu_p1)").c_str());
    m_wk->factory(("expr::gaus_sigma_k_Xhh_m" + mX + "('TMath::Exp(gaus_sigma_k_p0 + gaus_sigma_k_p1 * " + mX + " + gaus_sigma_k_p2 * " + mX + " * " + mX + ")', gaus_sigma_k_p0, gaus_sigma_k_p1, gaus_sigma_k_p2)").c_str());
    m_wk->factory(("prod::gaus_sigma_Xhh_m" + mX + "(gaus_sigma_k_Xhh_m" + mX + ", CB_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("CBShape::CB_Xhh_m" + mX + "(mass, CB_mu_Xhh_m" + mX + ", CB_sigma_Xhh_m" + mX + ", CB_alpha_Xhh_m" + mX + ", CB_n_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("Gaussian::gaus_Xhh_m" + mX + "(mass, gaus_mu_Xhh_m" + mX + ", gaus_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("SUM::Xhh_m" + mX + "(CB_frac_Xhh_m" + mX +" * CB_Xhh_m" + mX + ", gaus_Xhh_m" + mX + ")").c_str());
    dynamic_cast<RooSimultaneous*>(m_wk->pdf("signal_model"))->addPdf(*m_wk->pdf(("Xhh_m" + mX).c_str()), mX.c_str());
  }

  RooCategory* SignalModel::mass_points() {
    return m_wk->cat("mass_points");
  }

  void SignalModel::fit(RooDataSet& data)
  {
    m_data = &data;
    MSG_INFO("Fitting PDFs to data");
    // m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(1));
    // m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(1));
  }

  void SignalModel::write(const std::string& output_file_name)
  {
    // m_wk->var("CB_alpha_p0")->setConstant();
    // // MSG_INFO("CB_alpha_p0: " << m_wk->var("CB_alpha_p0")->getVal());
    // m_wk->var("CB_alpha_p1")->setConstant();
    // // MSG_INFO("CB_alpha_p1: " << m_wk->var("CB_alpha_p1")->getVal());
    // m_wk->var("CB_alpha_p2")->setConstant();
    // // MSG_INFO("CB_alpha_p2: " << m_wk->var("CB_alpha_p2")->getVal());
    // m_wk->var("CB_mu_p0")->setConstant();
    // // MSG_INFO("CB_mu_p0: " << m_wk->var("CB_mu_p0")->getVal());
    // m_wk->var("CB_mu_p1")->setConstant();
    // // MSG_INFO("CB_mu_p1: " << m_wk->var("CB_mu_p1")->getVal());
    // m_wk->var("CB_n")->setConstant();
    // // MSG_INFO("CB_n: " << m_wk->var("CB_n")->getVal());
    // m_wk->var("CB_sigma_p0")->setConstant();
    // // MSG_INFO("CB_sigma_p0: " << m_wk->var("CB_sigma_p0")->getVal());
    // m_wk->var("CB_sigma_p1")->setConstant();
    // // MSG_INFO("CB_sigma_p1: " << m_wk->var("CB_sigma_p1")->getVal());
    // m_wk->var("gaus_mu_p0")->setConstant();
    // // MSG_INFO("gaus_mu_p0: " << m_wk->var("gaus_mu_p0")->getVal());
    // m_wk->var("gaus_mu_p1")->setConstant();
    // // MSG_INFO("gaus_mu_p1: " << m_wk->var("gaus_mu_p1")->getVal());
    // m_wk->var("gaus_sigma_p0")->setConstant();
    // // MSG_INFO("gaus_sigma_p0: " << m_wk->var("gaus_sigma_p0")->getVal());
    // m_wk->var("gaus_sigma_p1")->setConstant();
    // // MSG_INFO("gaus_sigma_p1: " << m_wk->var("gaus_sigma_p1")->getVal());
    // m_wk->var("CB_frac_p0")->setConstant();
    // // MSG_INFO("CB_frac_p0: " << m_wk->var("CB_frac_p0")->getVal());
    // m_wk->var("CB_frac_p1")->setConstant();
    // // MSG_INFO("CB_frac_p1: " << m_wk->var("CB_frac_p1")->getVal());
    m_wk->factory("mass_resonance[300, 0, 2000]");
    m_wk->factory("expr::CB_mu('CB_mu_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
    m_wk->factory("expr::CB_alpha('CB_alpha_p0 + CB_alpha_p1 * mass_resonance + CB_alpha_p2 * mass_resonance * mass_resonance', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2, mass_resonance)");
    m_wk->factory("expr::CB_sigma('CB_sigma_p0 + CB_sigma_p1 * mass_resonance + CB_sigma_p2 * mass_resonance * mass_resonance', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2, mass_resonance)");
    m_wk->factory("expr::CB_n('CB_n_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
    m_wk->factory("CBShape::CB_PDF(mass, CB_mu, CB_sigma, CB_alpha, CB_n)");
    m_wk->factory("expr::gaus_mu('gaus_mu_p0 + gaus_mu_p1 * mass_resonance', gaus_mu_p0, gaus_mu_p1, mass_resonance)");
    m_wk->factory("expr::gaus_sigma('gaus_sigma_k * CB_sigma', gaus_sigma_k, CB_sigma)");
    m_wk->factory("Gaussian::gaus_PDF(mass, gaus_mu, gaus_sigma)");
    m_wk->factory("expr::CB_frac('TMath::Min(TMath::Max(CB_frac_p0 + CB_frac_p1 * mass_resonance, 0.0), 1.0)', CB_frac_p0, CB_frac_p1, mass_resonance)");
    m_wk->factory("SUM::signal_PDF(CB_frac * CB_PDF, gaus_PDF)");
    MSG_INFO("Preparing to write workspace to " << output_file_name);
    m_wk->writeToFile(output_file_name.c_str(), false);
  }

  void SignalModel::plot()
  // void SignalModel::plot(std::map<std::string, RooDataSet*> dataset_map)
  {
    MSG_INFO("Plotting \033[1m" << m_mass_category << " mass " << m_tag_category << "-tag\033[0m");
    PlotStyle::EnsureAtlasStyle();
    // Plot all on one canvas
    TCanvas canvas("canvas", "canvas", 600, 600);
    RooPlot* frame = m_wk->var("mass")->frame();

    for (unsigned int idx = 0; idx < PlotStyle::resonance_masses(m_mass_category).size(); ++idx) {
      std::string mX(std::to_string(PlotStyle::resonance_masses(m_mass_category).at(idx)));
      int colour(PlotStyle::colours().at(idx));
      m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::MarkerColor(colour));
      m_wk->pdf("signal_model")->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(colour));
    }
    frame->Draw();
    frame->SetMinimum(0);
    canvas.Print(("plots/signal_model_" + m_mass_category + "Mass_" + m_tag_category + "tag_overall.pdf").c_str());
    MSG_INFO("Created plots/signal_model_" << m_mass_category << "Mass_" << m_tag_category << "tag_overall.pdf");

    // Plot each fit individually
    std::vector<std::string> best_fit;
    for (unsigned int idx = 0; idx < 7; ++idx) { best_fit.push_back(""); }
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) {
      // Construct useful variables
      std::string mX(std::to_string(resonance_mass));
      double mass_low(0.9 * resonance_mass), mass_high(1.1 * resonance_mass);
      int nBins(int(mass_high - mass_low));
      // Fit with simple CB+G
      RooDataSet* data_slice = dynamic_cast<RooDataSet*>(m_data->reduce(RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::SelectVars(m_wk->argSet("mass, weight"))));
      // set_initial_values(resonance_mass);
      // m_wk->pdf("single_signal_PDF")->getParameters(*data_slice)->Print("v");
      m_wk->pdf("single_signal_PDF")->fitTo(*data_slice, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::PrintLevel(1));
      m_wk->pdf("single_signal_PDF")->Print("v");
      MSG_INFO("single_EGE_mu: " << m_wk->var("single_EGE_mu")->getVal());
      MSG_INFO("single_EGE_sigma: " << m_wk->var("single_EGE_sigma")->getVal());
      MSG_INFO("single_EGE_kLo: " << m_wk->var("single_EGE_kLo")->getVal());
      MSG_INFO("single_EGE_kHi: " << m_wk->var("single_EGE_kHi")->getVal());
      // MSG_INFO("For resonance with mass \033[1m" << mX << " GeV\033[0m");
      // MSG_INFO("single_CB_mu: " << m_wk->var("single_CB_mu")->getVal());
      // MSG_INFO("single_CB_sigma: " << m_wk->var("single_CB_sigma")->getVal());
      // MSG_INFO("single_CB_alpha: " << m_wk->var("single_CB_alpha")->getVal());
      // MSG_INFO("single_CB_n: " << m_wk->var("single_CB_n")->getVal());
      // MSG_INFO("single_gaus_mu: " << m_wk->var("single_gaus_mu")->getVal());
      // MSG_INFO("single_gaus_sigma_k: " << m_wk->var("single_gaus_sigma_k")->getVal());
      // MSG_INFO("single_CB_frac: " << m_wk->var("single_CB_frac")->getVal());
      // best_fit.at(0) = best_fit.at(0) + ", " + std::to_string(m_wk->var("single_CB_mu")->getVal());
      // best_fit.at(1) = best_fit.at(1) + ", " + std::to_string(m_wk->var("single_CB_sigma")->getVal());
      // best_fit.at(2) = best_fit.at(2) + ", " + std::to_string(m_wk->var("single_CB_alpha")->getVal());
      // best_fit.at(3) = best_fit.at(3) + ", " + std::to_string(m_wk->var("single_CB_n")->getVal());
      // best_fit.at(4) = best_fit.at(4) + ", " + std::to_string(m_wk->var("single_gaus_mu")->getVal());
      // best_fit.at(5) = best_fit.at(5) + ", " + std::to_string(m_wk->var("single_gaus_sigma_k")->getVal());
      // best_fit.at(6) = best_fit.at(6) + ", " + std::to_string(m_wk->var("single_CB_frac")->getVal());
      // Construct frames
      RooPlot* frame = m_wk->var("mass")->frame(mass_low, mass_high, nBins);
      RooPlot* frame_ratio = m_wk->var("mass")->frame(mass_low, mass_high, nBins);
      m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(kBlack));
      m_wk->pdf("signal_model")->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(kRed));
      RooHist* pull_hist_full = frame->pullHist(); pull_hist_full->SetLineColor(kRed); pull_hist_full->SetMarkerColor(kRed);
      data_slice->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::Invisible()); // RooFit::MarkerColor(kViolet)
      m_wk->pdf("single_signal_PDF")->plotOn(frame, RooFit::LineColor(kBlue));
      RooHist* pull_hist_simple = frame->pullHist(); pull_hist_simple->SetLineColor(kBlue); pull_hist_simple->SetMarkerColor(kBlue);
      // Construct a histogram with the pulls of the data w.r.t the curve
      frame_ratio->addPlotable(pull_hist_full, "P");
      frame_ratio->addPlotable(pull_hist_simple, "P");
      // Setup canvas
      TCanvas canvas("canvas", "canvas", 600, 600);
      gPad->SetBottomMargin(0.35);
      frame->SetMinimum(1e-10);
      frame->GetXaxis()->SetLabelOffset(-100);
      frame->Draw();
      TLatex textBox; textBox.SetNDC(); textBox.SetTextFont(42); textBox.SetTextSize(0.02);
      textBox.SetTextColor(kBlue);
      textBox.DrawLatex(0.2, 0.90, ("#mu_{EGE} = " + PlotStyle::to_string(m_wk->var("single_EGE_mu")->getVal(), 2)).c_str());
      textBox.DrawLatex(0.2, 0.86, ("#sigma_{EGE} = " + PlotStyle::to_string(m_wk->var("single_EGE_sigma")->getVal(), 2)).c_str());
      textBox.DrawLatex(0.2, 0.82, ("kLo_{EGE} = " + PlotStyle::to_string(m_wk->var("single_EGE_kLo")->getVal(), 2)).c_str());
      textBox.DrawLatex(0.2, 0.78, ("kHi_{EGE} = " + PlotStyle::to_string(m_wk->var("single_EGE_kHi")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.90, ("#mu_{DSCB} = " + PlotStyle::to_string(m_wk->var("single_DSCB_mu")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.86, ("#sigma_{DSCB} = " + PlotStyle::to_string(m_wk->var("single_DSCB_sigma")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.82, ("#alphaLo_{DSCB} = " + PlotStyle::to_string(m_wk->var("single_DSCB_alphaLo")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.78, ("nLo_{DSCB} = " + PlotStyle::to_string(m_wk->var("single_DSCB_nLo")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.74, ("#alphaHi_{DSCB} = " + PlotStyle::to_string(m_wk->var("single_DSCB_alphaHi")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.70, ("nHi_{DSCB} = " + PlotStyle::to_string(m_wk->var("single_DSCB_nHi")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.90, ("#mu_{CB} = " + PlotStyle::to_string(m_wk->var("single_CB_mu")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.86, ("#sigma_{CB} = " + PlotStyle::to_string(m_wk->var("single_CB_sigma")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.82, ("#alpha_{CB} = " + PlotStyle::to_string(m_wk->var("single_CB_alpha")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.78, ("n_{CB} = " + PlotStyle::to_string(m_wk->var("single_CB_n")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.74, ("f_{CB} = " + PlotStyle::to_string(m_wk->var("single_CB_frac")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.70, ("#mu_{GA} = " + PlotStyle::to_string(m_wk->var("single_gaus_mu")->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.2, 0.66, ("#sigma_{GA} = " + PlotStyle::to_string(m_wk->function("single_gaus_sigma")->getVal(), 2)).c_str());
      // textBox.SetTextColor(kRed);
      // textBox.DrawLatex(0.8, 0.90, ("#mu_{CB} = " + PlotStyle::to_string(m_wk->function(("CB_mu_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.8, 0.86, ("#sigma_{CB} = " + PlotStyle::to_string(m_wk->function(("CB_sigma_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.8, 0.82, ("#alpha_{CB} = " + PlotStyle::to_string(m_wk->function(("CB_alpha_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.8, 0.78, ("n_{CB} = " + PlotStyle::to_string(m_wk->function(("CB_n_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.8, 0.74, ("f_{CB} = " + PlotStyle::to_string(m_wk->function(("CB_frac_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.8, 0.70, ("#mu_{GA} = " + PlotStyle::to_string(m_wk->function(("gaus_mu_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      // textBox.DrawLatex(0.8, 0.66, ("#sigma_{GA} = " + PlotStyle::to_string(m_wk->function(("gaus_sigma_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
      TPad pad_bottom("pad_bottom", "pad_bottom", 0.0, 0.0, 1.0, 1.0);
      pad_bottom.SetTopMargin(0.65);
      pad_bottom.SetFillStyle(0);
      pad_bottom.Draw();
      pad_bottom.cd();
      frame_ratio->GetYaxis()->SetTitle("Pull");
      frame_ratio->GetYaxis()->SetTitleOffset(1.8);
      frame_ratio->GetYaxis()->SetRangeUser(-4.5, 4.5);
      frame_ratio->GetYaxis()->SetNdivisions(5, 4, 0);
      frame_ratio->Draw();
      TLine l(frame_ratio->GetXaxis()->GetXmin(), 0.0, frame_ratio->GetXaxis()->GetXmax(), 0.0);
      l.SetLineColor(kRed); l.SetLineWidth(2); l.SetLineStyle(kDashed); l.Draw();
      canvas.Print(("plots/" + m_mass_category + "Mass_" + m_tag_category + "tag/signal_model_" + m_mass_category + "Mass_" + m_tag_category + "tag_mX_" + mX + ".pdf").c_str());
      MSG_INFO("Created plots/" << m_mass_category << "Mass_" << m_tag_category << "tag/signal_model_" << m_mass_category << "Mass_" << m_tag_category << "tag_mX_" << mX << ".pdf");
    }
    MSG_INFO("CB_mu: " << best_fit.at(0));
    MSG_INFO("CB_sigma: " << best_fit.at(1));
    MSG_INFO("CB_alpha: " << best_fit.at(2));
    MSG_INFO("CB_n: " << best_fit.at(3));
    MSG_INFO("gaus_mu: " << best_fit.at(4));
    MSG_INFO("gaus_sigma_k: " << best_fit.at(5));
    MSG_INFO("CB_frac: " << best_fit.at(6));
  }

  void SignalModel::set_initial_values(const int& resonance_mass) {
    double CB_mu(0.5), CB_sigma(0.5), CB_alpha(0.5), CB_n(0.5), CB_frac(0.5), gaus_mu(0.5), gaus_sigma(0.5);
    if (m_mass_category == "low") {
      if (m_tag_category == "0") {
        if (resonance_mass == 260) { CB_mu = 260.33; CB_sigma = 3.90; CB_alpha = 2.90; CB_n = 10; CB_frac = 0.72; gaus_mu = 273.54; gaus_sigma = 9.15; }
        else if (resonance_mass == 275) { CB_mu = 275.13; CB_sigma = 4.74; CB_alpha = 0.94; CB_n = 10; CB_frac = 0.58; gaus_mu = 282.25; gaus_sigma = 15.22; }
        else if (resonance_mass == 300) { CB_mu = 300.09; CB_sigma = 5.45; CB_alpha = 0.48; CB_n = 10; CB_frac = 0.72; gaus_mu = 314.30; gaus_sigma = 11.59; }
        else if (resonance_mass == 325) { CB_mu = 325.00; CB_sigma = 6.14; CB_alpha = 0.46; CB_n = 10; CB_frac = 0.70; gaus_mu = 339.98; gaus_sigma = 12.39; }
        else if (resonance_mass == 350) { CB_mu = 350.73; CB_sigma = 6.51; CB_alpha = 1.37; CB_n = 10; CB_frac = 0.39; gaus_mu = 348.59; gaus_sigma = 24.40; }
        else if (resonance_mass == 400) { CB_mu = 400.94; CB_sigma = 6.11; CB_alpha = 0.31; CB_n = 10; CB_frac = 0.63; gaus_mu = 416.13; gaus_sigma = 15.24; }
        else { return; }
      }
    }
    m_wk->var("single_CB_mu")->setVal(CB_mu); //m_wk->var("single_CB_mu")->setConstant();
    m_wk->var("single_CB_sigma")->setVal(CB_sigma); //m_wk->var("single_CB_sigma")->setConstant();
    m_wk->var("single_CB_alpha")->setVal(CB_alpha); //m_wk->var("single_CB_alpha")->setConstant();
    m_wk->var("single_CB_n")->setVal(CB_n); //m_wk->var("single_CB_n")->setConstant();
    m_wk->var("single_CB_frac")->setVal(CB_frac); //m_wk->var("single_CB_frac")->setConstant();
    m_wk->var("single_gaus_sigma_k")->setVal(gaus_sigma / CB_sigma); //m_wk->var("single_gaus_sigma_k")->setConstant();
    m_wk->var("single_gaus_mu")->setVal(gaus_mu); //m_wk->var("single_gaus_mu")->setConstant();
  }

}
