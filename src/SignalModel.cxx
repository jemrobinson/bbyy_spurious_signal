// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "SignalModel.h"
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
    // Crystal Ball centre
    m_wk->factory("CB_alpha_p0[95, -INF, +INF]");
    m_wk->factory("CB_alpha_p1[-0.6, -10, 10]");
    m_wk->factory("CB_alpha_p2[0.001, -1, 1]");
    m_wk->factory("CB_mu_p0[4.7, -INF, +INF]");
    m_wk->factory("CB_mu_p1[0.003, 0.0, 10.0]");
    // m_wk->factory("CB_mu_p2[1e-5, -1, 1]");
    m_wk->factory("CB_n_p0[7, -INF, +INF]");
    m_wk->factory("CB_n_p1[-0.01, -10, 10]");
    m_wk->factory("CB_n_p2[0.0, -1, 1]");
    m_wk->factory("CB_sigma_p0[2, -INF, +INF]");
    m_wk->factory("CB_sigma_p1[0.01, -10, 10]");
    m_wk->factory("CB_sigma_p2[0.0, -1, 1]");
    // Gaussian tails
    m_wk->factory("gaus_mu_p0[100, -INF, +INF]");
    m_wk->factory("gaus_mu_p1[0.7, 0.0, 10.0]");
    m_wk->factory("gaus_mu_p2[1e-5, -1, 1]");
    m_wk->factory("gaus_sigma_k[9, 1, +INF]");
    // Combination
    m_wk->factory("CB_frac_p0[0.14, -100, 100]");
    m_wk->factory("CB_frac_p1[-0.002, -10, 10]");
    // Set up category to distinguish samples
    RooCategory mass_points("mass_points", "mass_points");
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) { mass_points.defineType(std::to_string(resonance_mass).c_str()); }
    m_wk->import(mass_points);
  }

  void SignalModel::build_simultaneous_PDF(RooRealVar& mass)
  {
    m_wk->import(mass);
    m_wk->factory("Simultaneous::signal_model(mass_points)");
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category) ) { add_mass_point(resonance_mass); }
  }

  void SignalModel::add_mass_point(const int& resonance_mass) {
    std::string mX(std::to_string(resonance_mass));
    m_wk->factory(("expr::CB_mu_Xhh_m" + mX + "('CB_mu_p0 + CB_mu_p1 * " + mX + "', CB_mu_p0, CB_mu_p1)").c_str());
    // m_wk->factory(("expr::CB_mu_Xhh_m" + mX + "('CB_mu_p0 + CB_mu_p1 * " + mX + " + CB_mu_p2 * " + mX + " * " + mX + "', CB_mu_p0, CB_mu_p1, CB_mu_p2)").c_str());
    m_wk->factory(("expr::CB_alpha_Xhh_m" + mX + "('CB_alpha_p0 + CB_alpha_p1 * " + mX + " + CB_alpha_p2 * " + mX + " * " + mX + "', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2)").c_str());
    m_wk->factory(("expr::CB_sigma_Xhh_m" + mX + "('CB_sigma_p0 + CB_sigma_p1 * " + mX + " + CB_sigma_p2 * " + mX + " * " + mX + "', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2)").c_str());
    // m_wk->factory(("expr::CB_n_Xhh_m" + mX + "('CB_n_p0 + CB_n_p1 * " + mX + "', CB_n_p0, CB_n_p1)").c_str());
    m_wk->factory(("expr::CB_n_Xhh_m" + mX + "('CB_n_p0 + CB_n_p1 * " + mX + " + CB_n_p2 * " + mX + " * " + mX + "', CB_n_p0, CB_n_p1, CB_n_p2)").c_str());
    m_wk->factory(("CBShape::CB_Xhh_m" + mX + "(mass, CB_mu_Xhh_m" + mX + ", CB_sigma_Xhh_m" + mX + ", CB_alpha_Xhh_m" + mX + ", CB_n_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("expr::gaus_mu_Xhh_m" + mX + "('gaus_mu_p0 + gaus_mu_p1 * " + mX + "', gaus_mu_p0, gaus_mu_p1)").c_str());
    // m_wk->factory(("expr::gaus_mu_Xhh_m" + mX + "('gaus_mu_p0 + gaus_mu_p1 * " + mX + " + gaus_mu_p2 * " + mX + " * " + mX + "', gaus_mu_p0, gaus_mu_p1, gaus_mu_p2)").c_str());
    m_wk->factory(("expr::gaus_sigma_Xhh_m" + mX + "('gaus_sigma_k * CB_sigma_Xhh_m" + mX + "', gaus_sigma_k, CB_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("Gaussian::gaus_Xhh_m" + mX + "(mass, gaus_mu_Xhh_m" + mX + ", gaus_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("expr::CB_frac_Xhh_m" + mX + "('TMath::Min(TMath::Max(CB_frac_p0 + CB_frac_p1 * " + mX + ", 0.0), 1.0)', CB_frac_p0, CB_frac_p1)").c_str());
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
    m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
    m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1));
  }

  // std::pair<double, double> SignalModel::two_sigma_window(const std::string& tag_category, const int& resonance_mass) {
  //   if (tag_category == "0") {
  //     if (mass_category == ""
  //   }
  // }

  void SignalModel::write(const std::string& output_file_name)
  {
    m_wk->var("CB_alpha_p0")->setConstant();
    MSG_INFO("CB_alpha_p0: " << m_wk->var("CB_alpha_p0")->getVal());
    m_wk->var("CB_alpha_p1")->setConstant();
    MSG_INFO("CB_alpha_p1: " << m_wk->var("CB_alpha_p1")->getVal());
    m_wk->var("CB_alpha_p2")->setConstant();
    MSG_INFO("CB_alpha_p2: " << m_wk->var("CB_alpha_p2")->getVal());
    m_wk->var("CB_mu_p0")->setConstant();
    MSG_INFO("CB_mu_p0: " << m_wk->var("CB_mu_p0")->getVal());
    m_wk->var("CB_mu_p1")->setConstant();
    MSG_INFO("CB_mu_p1: " << m_wk->var("CB_mu_p1")->getVal());
    m_wk->var("CB_n_p0")->setConstant();
    MSG_INFO("CB_n_p0: " << m_wk->var("CB_n_p0")->getVal());
    m_wk->var("CB_n_p1")->setConstant();
    MSG_INFO("CB_n_p1: " << m_wk->var("CB_n_p1")->getVal());
    m_wk->var("CB_sigma_p0")->setConstant();
    MSG_INFO("CB_sigma_p0: " << m_wk->var("CB_sigma_p0")->getVal());
    m_wk->var("CB_sigma_p1")->setConstant();
    MSG_INFO("CB_sigma_p1: " << m_wk->var("CB_sigma_p1")->getVal());
    m_wk->var("CB_sigma_p2")->setConstant();
    MSG_INFO("CB_sigma_p2: " << m_wk->var("CB_sigma_p2")->getVal());
    m_wk->var("gaus_mu_p0")->setConstant();
    MSG_INFO("gaus_mu_p0: " << m_wk->var("gaus_mu_p0")->getVal());
    m_wk->var("gaus_mu_p1")->setConstant();
    MSG_INFO("gaus_mu_p1: " << m_wk->var("gaus_mu_p1")->getVal());
    m_wk->var("gaus_sigma_k")->setConstant();
    MSG_INFO("gaus_sigma_k: " << m_wk->var("gaus_sigma_k")->getVal());
    m_wk->var("CB_frac_p0")->setConstant();
    MSG_INFO("CB_frac_p0: " << m_wk->var("CB_frac_p0")->getVal());
    m_wk->var("CB_frac_p1")->setConstant();
    MSG_INFO("CB_frac_p1: " << m_wk->var("CB_frac_p1")->getVal());
    m_wk->factory("mass_resonance[300, 0, 2000]");
    m_wk->factory("expr::CB_mu('CB_mu_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
    m_wk->factory("expr::CB_alpha('CB_alpha_p0 + CB_alpha_p1 * mass_resonance + CB_alpha_p2 * mass_resonance * mass_resonance', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2, mass_resonance)");
    m_wk->factory("expr::CB_sigma('CB_sigma_p0 + CB_sigma_p1 * mass_resonance + CB_sigma_p2 * mass_resonance * mass_resonance', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2, mass_resonance)");
    m_wk->factory("expr::CB_n('CB_n_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
    m_wk->factory("CBShape::CB_PDF(mass, CB_mu, CB_sigma, CB_alpha, CB_n)");
    // m_wk->factory("expr::gaus_mu('gaus_mu_k * CB_mu', gaus_mu_k, CB_mu)");
    m_wk->factory("expr::gaus_mu('gaus_mu_p0 + gaus_mu_p1 * mass_resonance', gaus_mu_p0, gaus_mu_p1, mass_resonance)");
    m_wk->factory("expr::gaus_sigma('gaus_sigma_k * CB_sigma', gaus_sigma_k, CB_sigma)");
    m_wk->factory("Gaussian::gaus_PDF(mass, gaus_mu, gaus_sigma)");
    m_wk->factory("expr::CB_frac('TMath::Min(TMath::Max(CB_frac_p0 + CB_frac_p1 * mass_resonance, 0.0), 1.0)', CB_frac_p0, CB_frac_p1, mass_resonance)");
    m_wk->factory("SUM::signal_PDF(CB_frac * CB_PDF, gaus_PDF)");
    MSG_INFO("Preparing to write workspace to " << output_file_name);
    m_wk->writeToFile(output_file_name.c_str(), false);
  }

  void SignalModel::plot()
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
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) {
      std::string mX(std::to_string(resonance_mass));
      double sigma(m_wk->function(("CB_sigma_Xhh_m" + mX).c_str())->getVal());
      double x_low(std::max(resonance_mass - 4 * sigma, 245.)), x_high(resonance_mass + 4 * sigma);
      RooPlot* frame = m_wk->var("mass")->frame(RooFit::Range(x_low, x_high));
      m_data->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::MarkerColor(kBlack));
      m_wk->pdf("signal_model")->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(kRed));
      // Leo's PDF
      // // m_wk->factory("Leo_CB_alpha_p0[20, -INF, +INF]");
      // // m_wk->factory("Leo_CB_alpha_p1[-0.005, -10, 10]");
      // // m_wk->factory("expr::Leo_CB_alpha('CB_alpha_p0 + CB_alpha_p1 * mass_resonance + CB_alpha_p2 * mass_resonance * mass_resonance', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2, mass_resonance)");
      // m_wk->factory("expr::CB_mu('CB_mu_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
      // m_wk->factory("expr::CB_sigma('CB_sigma_p0 + CB_sigma_p1 * mass_resonance + CB_sigma_p2 * mass_resonance * mass_resonance', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2, mass_resonance)");
      // m_wk->factory("expr::CB_n('CB_n_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
      // m_wk->factory("CBShape::CB_PDF(mass, CB_mu, CB_sigma, CB_alpha, CB_n)");
      // // m_wk->factory("expr::gaus_mu('gaus_mu_k * CB_mu', gaus_mu_k, CB_mu)");
      // m_wk->factory("expr::gaus_mu('gaus_mu_p0 + gaus_mu_p1 * mass_resonance', gaus_mu_p0, gaus_mu_p1, mass_resonance)");
      // m_wk->factory("expr::gaus_sigma('gaus_sigma_k * CB_sigma', gaus_sigma_k, CB_sigma)");
      // m_wk->factory("Gaussian::gaus_PDF(mass, gaus_mu, gaus_sigma)");
      // m_wk->factory("expr::CB_frac('TMath::Min(TMath::Max(CB_frac_p0 + CB_frac_p1 * mass_resonance, 0.0), 1.0)', CB_frac_p0, CB_frac_p1, mass_resonance)");
      // m_wk->factory("SUM::signal_PDF(CB_frac * CB_PDF, gaus_PDF)");
      // Construct a histogram with the pulls of the data w.r.t the curve
      RooPlot* frame_ratio = m_wk->var("mass")->frame(RooFit::Range(x_low, x_high));
      RooHist* pull_hist = frame->pullHist();
      frame_ratio->addPlotable(pull_hist, "P");
      // Setup canvas
      TCanvas canvas("canvas", "canvas", 600, 600);
      gPad->SetBottomMargin(0.35);
      frame->SetMinimum(1e-10);
      frame->GetXaxis()->SetLabelOffset(-100);
      frame->Draw();
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
  }
}
