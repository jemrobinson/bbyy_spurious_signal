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
    // Crystal Ball
    m_wk->factory("CB_mu_p0[2.2, -50, 50]");
    m_wk->factory("CB_mu_p1[1, 0.8, 1.2]");
    m_wk->factory("CB_sigma_p0[30, -1000, 1000]");
    m_wk->factory("CB_sigma_p1[0.06, -100, 100]");
    m_wk->factory("CB_alpha_p0[69, -500, 500]");
    m_wk->factory("CB_alpha_p1[-0.45, -5, 5]");
    m_wk->factory("CB_alpha_p2[0.0007, -0.5, 0.5]");
    m_wk->factory("CB_n[10, 8, 12]");
    // Combination
    m_wk->factory("CB_frac_p0[0.144208, -50, 50]");
    m_wk->factory("CB_frac_p1[-0.00194108, -5, 5]");
    // Gaussian
    m_wk->factory("gaus_mu_p0[0.83, -50, 50]");
    m_wk->factory("gaus_mu_p1[1, 0.5, 1.5]");
    m_wk->factory("gaus_sigma_p0[-40, -1000, 1000]");
    m_wk->factory("gaus_sigma_p1[0.25, -100, 100]");
    m_wk->factory("gaus_sigma_p2[-0.0003, -10, 10]");
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
    // Also build a simple CB+G shape
    m_wk->factory("simple_CB_mu[260.325, 0, +INF]");
    m_wk->factory("simple_CB_sigma[6.5, 0, +INF]");
    m_wk->factory("simple_CB_alpha[2.90303, 0, +INF]");
    m_wk->factory("simple_CB_n[10]");
    m_wk->factory("simple_gaus_mu[273.538, 0, +INF]");
    m_wk->factory("simple_gaus_sigma_k[6.5, 1, +INF]");
    m_wk->factory("simple_CB_frac[0.7, 0.0, 1.0]");
    m_wk->factory("prod::simple_gaus_sigma(simple_gaus_sigma_k, simple_CB_sigma)");
    m_wk->factory("CBShape::simple_CB(mass, simple_CB_mu, simple_CB_sigma, simple_CB_alpha, simple_CB_n)");
    m_wk->factory("Gaussian::simple_gaus(mass, simple_gaus_mu, simple_gaus_sigma)");
    m_wk->factory("SUM::simple_signal_PDF(simple_CB_frac * simple_CB, simple_gaus)");
  }

  void SignalModel::add_mass_point(const int& resonance_mass) {
    std::string mX(std::to_string(resonance_mass));
    m_wk->factory(("expr::CB_mu_Xhh_m" + mX + "('CB_mu_p0 + CB_mu_p1 * " + mX + "', CB_mu_p0, CB_mu_p1)").c_str());
    m_wk->factory(("expr::CB_sigma_Xhh_m" + mX + "('CB_sigma_p0 + CB_sigma_p1 * " + mX + "', CB_sigma_p0, CB_sigma_p1)").c_str());
    m_wk->factory(("expr::CB_alpha_Xhh_m" + mX + "('TMath::Exp(CB_alpha_p0 + CB_alpha_p1 * " + mX + " + CB_alpha_p2 * " + mX + " * " + mX + ")', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2)").c_str());
    m_wk->factory(("expr::CB_frac_Xhh_m" + mX + "('TMath::Exp(CB_frac_p0 + CB_frac_p1 * " + mX + ")', CB_frac_p0, CB_frac_p1)").c_str());
    m_wk->factory(("expr::gaus_mu_Xhh_m" + mX + "('gaus_mu_p0 + gaus_mu_p1 * " + mX + "', gaus_mu_p0, gaus_mu_p1)").c_str());
    m_wk->factory(("expr::gaus_sigma_Xhh_m" + mX + "('gaus_sigma_p0 + gaus_sigma_p1 * " + mX + " + gaus_sigma_p2 * " + mX + " * " + mX + "', gaus_sigma_p0, gaus_sigma_p1, gaus_sigma_p2)").c_str());
    m_wk->factory(("CBShape::CB_Xhh_m" + mX + "(mass, CB_mu_Xhh_m" + mX + ", CB_sigma_Xhh_m" + mX + ", CB_alpha_Xhh_m" + mX + ", CB_n)").c_str());
    m_wk->factory(("Gaussian::gaus_Xhh_m" + mX + "(mass, gaus_mu_Xhh_m" + mX + ", gaus_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("SUM::Xhh_m" + mX + "(CB_frac_Xhh_m" + mX +" * CB_Xhh_m" + mX + ", gaus_Xhh_m" + mX + ")").c_str());
    dynamic_cast<RooSimultaneous*>(m_wk->pdf("signal_model"))->addPdf(*m_wk->pdf(("Xhh_m" + mX).c_str()), mX.c_str());
    // m_wk->factory(("expr::gaus_sigma_Xhh_m" + mX + "('(gaus_sigma_k_p0 + gaus_sigma_k_p1 * " + mX + " + gaus_sigma_k_p2 * " + mX + " * " + mX + ") * CB_sigma_Xhh_m" + mX + "', gaus_sigma_k_p0, gaus_sigma_k_p1, gaus_sigma_k_p2, CB_sigma_Xhh_m" + mX +")").c_str());
  }

  RooCategory* SignalModel::mass_points() {
    return m_wk->cat("mass_points");
  }

  void SignalModel::fit(RooDataSet& data)
  {
    m_data = &data;
    MSG_INFO("Fitting PDFs to data");
    m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(1));
    m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(1));
  }

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
    m_wk->var("CB_n")->setConstant();
    MSG_INFO("CB_n: " << m_wk->var("CB_n")->getVal());
    m_wk->var("CB_sigma_p0")->setConstant();
    MSG_INFO("CB_sigma_p0: " << m_wk->var("CB_sigma_p0")->getVal());
    m_wk->var("CB_sigma_p1")->setConstant();
    MSG_INFO("CB_sigma_p1: " << m_wk->var("CB_sigma_p1")->getVal());
    m_wk->var("gaus_mu_p0")->setConstant();
    MSG_INFO("gaus_mu_p0: " << m_wk->var("gaus_mu_p0")->getVal());
    m_wk->var("gaus_mu_p1")->setConstant();
    MSG_INFO("gaus_mu_p1: " << m_wk->var("gaus_mu_p1")->getVal());
    m_wk->var("gaus_sigma_p0")->setConstant();
    MSG_INFO("gaus_sigma_p0: " << m_wk->var("gaus_sigma_p0")->getVal());
    m_wk->var("gaus_sigma_p1")->setConstant();
    MSG_INFO("gaus_sigma_p1: " << m_wk->var("gaus_sigma_p1")->getVal());
    // m_wk->var("gaus_sigma_k")->setConstant();
    // MSG_INFO("gaus_sigma_k: " << m_wk->var("gaus_sigma_k")->getVal());
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

  // void SignalModel::plot()
  void SignalModel::plot(std::map<std::string, RooDataSet*> dataset_map)
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
      RooPlot* frame = m_wk->var("mass")->frame(RooFit::Range(PlotStyle::mass_low(resonance_mass), PlotStyle::mass_high(resonance_mass)));
      m_data->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::MarkerColor(kBlack));
      m_wk->pdf("signal_model")->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(kRed));
      RooHist* pull_hist = frame->pullHist();
      // Print values
      m_wk->pdf("simple_signal_PDF")->fitTo(*dataset_map[mX], RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
      MSG_INFO(mX);
      MSG_INFO("simple_CB_mu:" << m_wk->var("simple_CB_mu")->getVal() << " vs. " << m_wk->function(("CB_mu_Xhh_m" + mX).c_str())->getVal());
      MSG_INFO("simple_CB_sigma:" << m_wk->var("simple_CB_sigma")->getVal() << " vs. " << m_wk->function(("CB_sigma_Xhh_m" + mX).c_str())->getVal());
      MSG_INFO("simple_CB_alpha:" << m_wk->var("simple_CB_alpha")->getVal() << " vs. " << m_wk->function(("CB_alpha_Xhh_m" + mX).c_str())->getVal());
      MSG_INFO("simple_CB_n:" << m_wk->var("simple_CB_n")->getVal() << " vs. " << m_wk->var("CB_n")->getVal());
      MSG_INFO("simple_gaus_mu:" << m_wk->var("simple_gaus_mu")->getVal() << " vs. " << m_wk->function(("gaus_mu_Xhh_m" + mX).c_str())->getVal());
      MSG_INFO("simple_gaus_sigma:" << m_wk->function("simple_gaus_sigma")->getVal() << " vs. " << m_wk->function(("gaus_sigma_Xhh_m" + mX).c_str())->getVal());
      MSG_INFO("simple_CB_frac:" << m_wk->var("simple_CB_frac")->getVal() << " vs. " << m_wk->function(("CB_frac_Xhh_m" + mX).c_str())->getVal());
      dataset_map[mX]->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(kBlack));
      m_wk->pdf("simple_signal_PDF")->plotOn(frame, RooFit::LineColor(kBlue));
      // Construct a histogram with the pulls of the data w.r.t the curve
      RooPlot* frame_ratio = m_wk->var("mass")->frame(RooFit::Range(PlotStyle::mass_low(resonance_mass), PlotStyle::mass_high(resonance_mass)));
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
