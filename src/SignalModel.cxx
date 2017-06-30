#include "Logger.h"
#include "PlotStyle.h"
#include "RooAbsPdf.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooCBShape.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "SignalModel.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFile.h"


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
    m_wk->factory("CB_alpha_p0[-1, -10, 10]");
    m_wk->factory("CB_alpha_p1[0.0001, -0.1, 0.1]");
    m_wk->factory("CB_alpha_p2[-1e07, -0.0001, 0.0001]");
    m_wk->factory("CB_mu_p0[0.75, -10, 10]");
    m_wk->factory("CB_mu_p1[1, 0.7, 1.3]");
    m_wk->factory("CB_n[5, 0, 20]");
    m_wk->factory("CB_sigma_p0[-10, -100, 100]");
    m_wk->factory("CB_sigma_p1[0.05, -10, 10]");
    m_wk->factory("CB_sigma_p2[-1e-5, -1, 1]");
    // Gaussian tails
    m_wk->factory("gaus_mu_k[1, 0.8, 1.2]");
    m_wk->factory("gaus_sigma_k[4, 1, 20]");
    // Combination
    m_wk->factory("CB_frac[0.95, 0, 1.0]");
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
    m_wk->factory(("expr::CB_alpha_Xhh_m" + mX + "('CB_alpha_p0 + CB_alpha_p1 * " + mX + " + CB_alpha_p2 * " + mX + " * " + mX + "', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2)").c_str());
    m_wk->factory(("expr::CB_sigma_Xhh_m" + mX + "('CB_sigma_p0 + CB_sigma_p1 * " + mX + " + CB_sigma_p2 * " + mX + " * " + mX + "', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2)").c_str());
    m_wk->factory(("CBShape::CB_Xhh_m" + mX + "(mass, CB_mu_Xhh_m" + mX + ", CB_sigma_Xhh_m" + mX + ", CB_alpha_Xhh_m" + mX + ", CB_n)").c_str());
    m_wk->factory(("expr::gaus_mu_Xhh_m" + mX + "('gaus_mu_k * CB_mu_Xhh_m" + mX + "', gaus_mu_k, CB_mu_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("expr::gaus_sigma_Xhh_m" + mX + "('gaus_sigma_k * CB_sigma_Xhh_m" + mX + "', gaus_sigma_k, CB_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("Gaussian::gaus_Xhh_m" + mX + "(mass, gaus_mu_Xhh_m" + mX + ", gaus_sigma_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("SUM::Xhh_m" + mX + "(CB_frac * CB_Xhh_m" + mX + ", gaus_Xhh_m" + mX + ")").c_str());
    dynamic_cast<RooSimultaneous*>(m_wk->pdf("signal_model"))->addPdf(*m_wk->pdf(("Xhh_m" + mX).c_str()), mX.c_str());
  }

  RooCategory* SignalModel::mass_points() {
    return m_wk->cat("mass_points");
  }

  void SignalModel::fit(RooDataSet& data)
  {
    m_data = &data;
    MSG_INFO("Fitting PDFs to data");
    m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(false), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
    m_wk->pdf("signal_model")->fitTo(data, RooFit::SumW2Error(false), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1));
  }

  void SignalModel::write(const std::string& output_file_name)
  {
    m_wk->var("CB_alpha_p0")->setConstant();
    m_wk->var("CB_alpha_p1")->setConstant();
    m_wk->var("CB_alpha_p2")->setConstant();
    m_wk->var("CB_mu_p0")->setConstant();
    m_wk->var("CB_mu_p1")->setConstant();
    m_wk->var("CB_n")->setConstant();
    m_wk->var("CB_sigma_p0")->setConstant();
    m_wk->var("CB_sigma_p1")->setConstant();
    m_wk->var("CB_sigma_p2")->setConstant();
    m_wk->var("gaus_mu_k")->setConstant();
    m_wk->var("gaus_sigma_k")->setConstant();
    m_wk->var("CB_frac")->setConstant();
    m_wk->factory("mass_resonance[300, 0, 2000]");
    m_wk->factory("expr::CB_mu('CB_mu_p0 + CB_mu_p1 * mass_resonance', CB_mu_p0, CB_mu_p1, mass_resonance)");
    m_wk->factory("expr::CB_alpha('CB_alpha_p0 + CB_alpha_p1 * mass_resonance + CB_alpha_p2 * mass_resonance * mass_resonance', CB_alpha_p0, CB_alpha_p1, CB_alpha_p2, mass_resonance)");
    m_wk->factory("expr::CB_sigma('CB_sigma_p0 + CB_sigma_p1 * mass_resonance + CB_sigma_p2 * mass_resonance * mass_resonance', CB_sigma_p0, CB_sigma_p1, CB_sigma_p2, mass_resonance)");
    m_wk->factory("CBShape::CB_PDF(mass, CB_mu, CB_sigma, CB_alpha, CB_n)");
    m_wk->factory("expr::gaus_mu('gaus_mu_k * CB_mu', gaus_mu_k, CB_mu)");
    m_wk->factory("expr::gaus_sigma('gaus_sigma_k * CB_sigma', gaus_sigma_k, CB_sigma)");
    m_wk->factory("Gaussian::gaus_PDF(mass, gaus_mu, gaus_sigma)");
    m_wk->factory("SUM::signal_PDF(CB_frac * CB_PDF, gaus_PDF)");
    MSG_INFO("Preparing to write workspace to " << output_file_name);
    // m_wk->Print("v");
    m_wk->writeToFile(output_file_name.c_str(), false);
  }

  void SignalModel::plot()
  {
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
    canvas.Print(("plots/signal_model_" + m_mass_category + "Mass_" + m_tag_category + "tag_overall.pdf").c_str());
    MSG_INFO("Created \033[1m" << m_mass_category << " mass " << m_tag_category << "-tag\033[0m plot");

    // Plot each fit individually
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) {
      std::string mX(std::to_string(resonance_mass));
      TCanvas canvas("canvas", "canvas", 600, 600);
      RooPlot* frame = m_wk->var("mass")->frame(RooFit::Range(0.85 * resonance_mass, 1.15 * resonance_mass));
      m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::MarkerColor(kBlack));
      m_wk->pdf("signal_model")->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(PlotStyle::colours().at(0)));
      frame->Draw();
      canvas.Print(("plots/" + m_mass_category + "Mass_" + m_tag_category + "tag/signal_model_" + m_mass_category + "Mass_" + m_tag_category + "tag_mX_" + mX + ".pdf").c_str());
      MSG_INFO("Created plots/" << m_mass_category << "Mass_" << m_tag_category << "tag/signal_model_" << m_mass_category << "Mass_" << m_tag_category << "tag_mX_" << mX << ".pdf");
    }
  }
}
