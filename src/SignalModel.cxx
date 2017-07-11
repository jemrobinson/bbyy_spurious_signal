// Local
#include "Logger.h"
#include "PlotStyle.h"
#include "SignalModel.h"
#include "ExpGausExpPDF.h"
// STL
#include <algorithm>
#include <iomanip>
// ROOT and RooFit
#include "RooCategory.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include <typeinfo>


namespace SpuriousSignal {
  /**
   * SignalModel constructor
   */
  SignalModel::SignalModel(const std::string& mass_category, const std::string& tag_category)
    : m_mass_category(mass_category)
    , m_tag_category(tag_category)
    , m_wk(0)
    , m_data(0)
    , m_models({"EGE"})//, "CBGA"})
    , m_model_metaparameters({{"CBGA", {"CBGA_alpha_CB_p0", "CBGA_alpha_CB_p1", "CBGA_mu_CB_p0", "CBGA_mu_CB_p1", "CBGA_mu_GA_p0", "CBGA_mu_GA_p1", "CBGA_n_CB_p0", "CBGA_n_CB_p1", "CBGA_sigma_CB_p0", "CBGA_sigma_CB_p1", "CBGA_sigma_GA_p0", "CBGA_sigma_GA_p1"}},
                          {"EGE",  {"EGE_mu_p0", "EGE_mu_p1", "EGE_sigma_p0", "EGE_sigma_p1", "EGE_kL_p0", "EGE_kL_p1", "EGE_kL_p2", "EGE_kH_p0", "EGE_kH_p1", "EGE_kH_p2"}}})
    , m_model_parameters({{"CBGA", {"CBGA_alpha_CB", "CBGA_f_CB", "CBGA_mu_CB", "CBGA_mu_GA", "CBGA_n_CB", "CBGA_sigma_CB", "CBGA_sigma_GA"}},
                          {"EGE",  {"EGE_mu", "EGE_sigma", "EGE_kL", "EGE_kH"}}})
  {
    m_wk = new RooWorkspace(("signal_model_" + mass_category + "Mass_" + tag_category + "tag").c_str());
    // Set up category to distinguish samples
    RooCategory mass_points("mass_points", "mass_points");
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) {
      mass_points.defineType(std::to_string(resonance_mass).c_str());
    }
    m_wk->import(mass_points);
    m_wk->factory("mass_resonance[300, 0, 2000]");
    // Simultaneous ExpGausExp
    if (std::find(m_models.begin(), m_models.end(), "EGE") != m_models.end()) {
      m_wk->factory(("EGE_mu_p0[" + std::to_string(m_mass_category == "low" ? 3.0 : 1.4) + ", -4, 4]").c_str());
      m_wk->factory(("EGE_mu_p1[" + std::to_string(m_mass_category == "low" ? 1.0 : 0.997) + ", 0.99, 1.01]").c_str());
      m_wk->factory(("EGE_sigma_p0[" + std::to_string(m_mass_category == "low" ? -5.0 : -4.0) + ", -50, 50]").c_str());
      m_wk->factory(("EGE_sigma_p1[" + std::to_string(m_mass_category == "low" ? 0.03 : 0.03) + ", -0.5, 0.5]").c_str());
      m_wk->factory(("EGE_kL_p0[" + std::to_string(m_mass_category == "low" ? 29.1 : 24.8) + ", -100, 100]").c_str());
      m_wk->factory(("EGE_kL_p1[" + std::to_string(m_mass_category == "low" ? -0.15 : -0.11) + ", -1, 1]").c_str());
      m_wk->factory(("EGE_kL_p2[" + std::to_string(m_mass_category == "low" ? 0.0002 : 0.0001) + ", -0.005, 0.005]").c_str());
      m_wk->factory(("EGE_kH_p0[" + std::to_string(m_mass_category == "low" ? 0.27 : -0.28) + ", -10, 10]").c_str());
      m_wk->factory(("EGE_kH_p1[" + std::to_string(m_mass_category == "low" ? 0.002 : 0.004) + ", -0.1, 0.1]").c_str());
      m_wk->factory(("EGE_kH_p2[" + std::to_string(m_mass_category == "low" ? 0.0 : -2.5e-6) + ", -0.001, 0.001]").c_str());
    }
    // Simultaneous Crystal Ball + Gaussian
    if (std::find(m_models.begin(), m_models.end(), "EGE") != m_models.end()) {
      m_wk->factory("CBGA_alpha_CB_p0[5, -10, 10]");
      m_wk->factory("CBGA_alpha_CB_p1[0.0, -0.01, 0.01]");
      // m_wk->factory(("CBGA_alpha_CB_p2[" + parameterised_CBGA(2) + ", -0.001, 0.001]").c_str());
      m_wk->factory("CBGA_f_CB_p0[0.8, -10, 10]");
      m_wk->factory("CBGA_f_CB_p1[0.0, -0.1, 0.1]");
      m_wk->factory("CBGA_mu_CB_p0[0.0, -5, 5]");

      m_wk->factory("CBGA_mu_CB_p1[1.0, 0.99, 1.01]");
      m_wk->factory("CBGA_mu_GA_p0[0.0, -5, 5]");
      m_wk->factory("CBGA_mu_GA_p1[1.0, 0.99, 1.01]");
      m_wk->factory("CBGA_n_CB_p0[8.0, -10, 10]");
      m_wk->factory("CBGA_n_CB_p1[0.0, -0.01, 0.01]");
      m_wk->factory("CBGA_sigma_CB_p0[5, -10, 10]");
      m_wk->factory("CBGA_sigma_CB_p1[0.0, -0.01, 0.01]");
      // m_wk->factory(("CBGA_sigma_CB_p2[" + parameterised_CBGA(13) + ", -0.0001, 0.0001]").c_str());
      m_wk->factory("CBGA_sigma_GA_p0[5, -10, 10]");
      m_wk->factory("CBGA_sigma_GA_p1[0.0, -0.01, 0.01]");
      // m_wk->factory(("CBGA_sigma_GA_p2[" + parameterised_CBGA(16) + ", -0.1, 0.1]").c_str());
    }
    // Leo's parameterised CB+GA
    m_wk->factory(("expr::Leo_CBGA_alpha_CB('" + parameterised_CBGA(0) + " + " + parameterised_CBGA(1) + " * mass_resonance + " + parameterised_CBGA(2) + " * mass_resonance * mass_resonance', mass_resonance)").c_str());
    m_wk->factory(("expr::Leo_CBGA_f_CB('TMath::Exp(" + parameterised_CBGA(3) + " + " + parameterised_CBGA(4) + " * mass_resonance)', mass_resonance)").c_str());
    m_wk->factory("Leo_CBGA_n_CB[10.0]");
    if (m_mass_category == "low") {
      m_wk->factory(("expr::Leo_CBGA_mu_CB('TMath::Exp(" + parameterised_CBGA(5) + " + " + parameterised_CBGA(6) + " * mass_resonance)', mass_resonance)").c_str());
      m_wk->factory(("expr::Leo_CBGA_mu_GA('TMath::Exp(" + parameterised_CBGA(8) + " + " + parameterised_CBGA(9) + " * mass_resonance)', mass_resonance)").c_str());
    } else {
      m_wk->factory(("expr::Leo_CBGA_mu_CB('" + parameterised_CBGA(5) + " + " + parameterised_CBGA(6) + " * mass_resonance + " + parameterised_CBGA(7) + " * mass_resonance * mass_resonance', mass_resonance)").c_str());
      m_wk->factory(("expr::Leo_CBGA_mu_GA('" + parameterised_CBGA(8) + " + " + parameterised_CBGA(9) + " * mass_resonance + " + parameterised_CBGA(10) + " * mass_resonance * mass_resonance', mass_resonance)").c_str());
    }
    m_wk->factory(("expr::Leo_CBGA_sigma_CB('" + parameterised_CBGA(11) + " + " + parameterised_CBGA(12) + " * mass_resonance + " + parameterised_CBGA(13) + " * mass_resonance * mass_resonance', mass_resonance)").c_str());
    m_wk->factory(("expr::Leo_CBGA_sigma_GA('" + parameterised_CBGA(14) + " + " + parameterised_CBGA(15) + " * mass_resonance + " + parameterised_CBGA(16) + " * mass_resonance * mass_resonance', mass_resonance)").c_str());
  }

  void SignalModel::build_simultaneous_PDF(RooRealVar& mass)
  {
    m_wk->import(mass);
    // Individual ExpGausExp (and import the code into the workspace)
    m_wk->factory("individual_EGE_mu[260, 0, 2000]");
    m_wk->factory("individual_EGE_sigma[4, 0.0, 100]");
    m_wk->factory("individual_EGE_kL[2, 0.0, 20]");
    m_wk->factory("individual_EGE_kH[2, 0.0, 20]");
    m_wk->factory("ExpGausExpPDF::individual_EGE(mass, individual_EGE_mu, individual_EGE_sigma, individual_EGE_kL, individual_EGE_kH)");
    m_wk->addClassDeclImportDir("include/");
    m_wk->importClassCode("ExpGausExpPDF*");
    // Individual Crystal Ball + Gaussian
    m_wk->factory("individual_CBGA_alpha_CB[0.1, 0.0, 10.0]");
    m_wk->factory("individual_CBGA_f_CB[0.15, 0.0, 1.0]");
    m_wk->factory("individual_CBGA_mu_CB[260, 0, 2000]");
    m_wk->factory("individual_CBGA_mu_GA[260, 0, 2000]");
    m_wk->factory("individual_CBGA_n_CB[10, 0, 200]");
    m_wk->factory("individual_CBGA_sigma_CB[4, 0.0, 40]");
    m_wk->factory("individual_CBGA_sigma_GA[4, 0.0, 40]");
    m_wk->factory("CBShape::individual_CBGA_CB(mass, individual_CBGA_mu_CB, individual_CBGA_sigma_CB, individual_CBGA_alpha_CB, individual_CBGA_n_CB)");
    m_wk->factory("Gaussian:individual_CBGA_GA(mass,individual_CBGA_mu_GA, individual_CBGA_sigma_GA)");
    m_wk->factory("SUM::individual_CBGA(individual_CBGA_f_CB * individual_CBGA_CB, individual_CBGA_GA)");
    // Leo's CB+GA
    m_wk->factory("CBShape::Leo_CBGA_CB(mass, Leo_CBGA_mu_CB, Leo_CBGA_sigma_CB, Leo_CBGA_alpha_CB, Leo_CBGA_n_CB)");
    m_wk->factory("Gaussian::Leo_CBGA_GA(mass, Leo_CBGA_mu_GA, Leo_CBGA_sigma_GA)");
    m_wk->factory("SUM::Leo_CBGA(Leo_CBGA_f_CB * Leo_CBGA_CB, Leo_CBGA_GA)");
    // Build a simultaneous ExpGausExp and CB+GA
    m_wk->factory("Simultaneous::simultaneous_EGE(mass_points)");
    m_wk->factory("Simultaneous::simultaneous_CBGA(mass_points)");
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category) ) {
      add_mass_point(resonance_mass);
    }
  }

  void SignalModel::add_mass_point(const int& resonance_mass) {
    std::string mX(std::to_string(resonance_mass));
    // Add ExpGausExp
    m_wk->factory(("expr::EGE_mu_Xhh_m" + mX + "('EGE_mu_p0 + EGE_mu_p1 * " + mX + "', EGE_mu_p0, EGE_mu_p1)").c_str());
    m_wk->factory(("expr::EGE_sigma_Xhh_m" + mX + "('TMath::Min(TMath::Max(EGE_sigma_p0 + EGE_sigma_p1 * " + mX + ", " + std::to_string(m_mass_category == "low" ? 3.0 : 6.0 ) + "), " + std::to_string(m_mass_category == "low" ? 10.0 : 18.0 ) + ")', EGE_sigma_p0, EGE_sigma_p1)").c_str());
    m_wk->factory(("expr::EGE_kL_Xhh_m" + mX + "('TMath::Min(TMath::Max(EGE_kL_p0 + EGE_kL_p1 * " + mX + " + EGE_kL_p2 * " + mX + " * " + mX + ", 0.1), " + std::to_string(m_mass_category == "low" ? 3.0 : 1.1) + ")', EGE_kL_p0, EGE_kL_p1, EGE_kL_p2)").c_str());
    m_wk->factory(("expr::EGE_kH_Xhh_m" + mX + "('TMath::Min(TMath::Max(EGE_kH_p0 + EGE_kH_p1 * " + mX + " + EGE_kH_p2 * " + mX + " * " + mX + ", 0.1), " + std::to_string(m_mass_category == "low" ? 3.0 : 0.9) + ")', EGE_kH_p0, EGE_kH_p1, EGE_kH_p2)").c_str());
    m_wk->factory(("ExpGausExpPDF::EGE_Xhh_m" + mX + "(mass, EGE_mu_Xhh_m" + mX +", EGE_sigma_Xhh_m" + mX + ", EGE_kL_Xhh_m" + mX + ", EGE_kH_Xhh_m" + mX + ")").c_str());
    dynamic_cast<RooSimultaneous*>(m_wk->pdf("simultaneous_EGE"))->addPdf(*m_wk->pdf(("EGE_Xhh_m" + mX).c_str()), mX.c_str());
    // Add Crystal Ball + Gaussian
    // 2 params
    m_wk->factory(("expr::CBGA_alpha_CB_Xhh_m" + mX + "('CBGA_alpha_CB_p0 + CBGA_alpha_CB_p1 * " + mX + "', CBGA_alpha_CB_p0, CBGA_alpha_CB_p1)").c_str());
    m_wk->factory(("expr::CBGA_f_CB_Xhh_m" + mX + "('CBGA_f_CB_p0 + CBGA_f_CB_p1 * " + mX + "', CBGA_f_CB_p0, CBGA_f_CB_p1)").c_str());
    m_wk->factory(("expr::CBGA_mu_CB_Xhh_m" + mX + "('CBGA_mu_CB_p0 + CBGA_mu_CB_p1 * " + mX + "', CBGA_mu_CB_p0, CBGA_mu_CB_p1)").c_str());
    m_wk->factory(("expr::CBGA_mu_GA_Xhh_m" + mX + "('CBGA_mu_GA_p0 + CBGA_mu_GA_p1 * " + mX + "', CBGA_mu_GA_p0, CBGA_mu_GA_p1)").c_str());
    m_wk->factory(("expr::CBGA_n_CB_Xhh_m" + mX + "('CBGA_n_CB_p0 + CBGA_n_CB_p1 * " + mX + "', CBGA_n_CB_p0, CBGA_n_CB_p1)").c_str());
    m_wk->factory(("expr::CBGA_sigma_CB_Xhh_m" + mX + "('CBGA_sigma_CB_p0 + CBGA_sigma_CB_p1 * " + mX + "', CBGA_sigma_CB_p0, CBGA_sigma_CB_p1)").c_str());
    m_wk->factory(("expr::CBGA_sigma_GA_Xhh_m" + mX + "('CBGA_sigma_GA_p0 + CBGA_sigma_GA_p1 * " + mX + "', CBGA_sigma_GA_p0, CBGA_sigma_GA_p1)").c_str());
    // 3 params
    // m_wk->factory(("expr::CBGA_sigma_GA_Xhh_m" + mX + "('CBGA_sigma_GA_p0 + CBGA_sigma_GA_p1 * " + mX + "', CBGA_sigma_GA_p0, CBGA_sigma_GA_p1)").c_str());
    // m_wk->factory(("expr::CBGA_alpha_CB_Xhh_m" + mX + "('CBGA_alpha_CB_p0 + CBGA_alpha_CB_p1 * " + mX + " + CBGA_alpha_CB_p2 * " + mX + " * " + mX + "', CBGA_alpha_CB_p0, CBGA_alpha_CB_p1, CBGA_alpha_CB_p2)").c_str());
    // m_wk->factory(("expr::CBGA_f_CB_Xhh_m" + mX + "('TMath::Exp(CBGA_f_CB_p0 + CBGA_f_CB_p1 * " + mX + ")', CBGA_f_CB_p0, CBGA_f_CB_p1)").c_str());
    // m_wk->factory(("expr::CBGA_mu_CB_Xhh_m" + mX + "('CBGA_mu_CB_p0 + CBGA_mu_CB_p1 * " + mX + " + CBGA_mu_CB_p2 * " + mX + " * " + mX + "', CBGA_mu_CB_p0, CBGA_mu_CB_p1, CBGA_mu_CB_p2)").c_str());
    // m_wk->factory(("expr::CBGA_mu_GA_Xhh_m" + mX + "('CBGA_mu_GA_p0 + CBGA_mu_GA_p1 * " + mX + " + CBGA_mu_GA_p2 * " + mX + " * " + mX + "', CBGA_mu_GA_p0, CBGA_mu_GA_p1, CBGA_mu_GA_p2)").c_str());
    // m_wk->factory(("expr::CBGA_n_CB_Xhh_m" + mX + "('TMath::Exp(CBGA_n_CB_p0 + CBGA_n_CB_p1 * " + mX + ")', CBGA_n_CB_p0, CBGA_n_CB_p1)").c_str());
    // m_wk->factory(("expr::CBGA_sigma_CB_Xhh_m" + mX + "('CBGA_sigma_CB_p0 + CBGA_sigma_CB_p1 * " + mX + " + CBGA_sigma_CB_p2 * " + mX + " * " + mX + "', CBGA_sigma_CB_p0, CBGA_sigma_CB_p1, CBGA_sigma_CB_p2)").c_str());
    // m_wk->factory(("expr::CBGA_sigma_GA_Xhh_m" + mX + "('CBGA_sigma_GA_p0 + CBGA_sigma_GA_p1 * " + mX + " + CBGA_sigma_GA_p2 * " + mX + " * " + mX + "', CBGA_sigma_GA_p0, CBGA_sigma_GA_p1, CBGA_sigma_GA_p2)").c_str());
    m_wk->factory(("CBShape::CBGA_CB_Xhh_m" + mX + "(mass, CBGA_mu_CB_Xhh_m" + mX + ", CBGA_sigma_CB_Xhh_m" + mX + ", CBGA_alpha_CB_Xhh_m" + mX + ", CBGA_n_CB_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("Gaussian::CBGA_GA_Xhh_m" + mX + "(mass, CBGA_mu_GA_Xhh_m" + mX + ", CBGA_sigma_GA_Xhh_m" + mX + ")").c_str());
    m_wk->factory(("SUM::CBGA_Xhh_m" + mX + "(CBGA_f_CB_Xhh_m" + mX + " * CBGA_CB_Xhh_m" + mX + ", CBGA_GA_Xhh_m" + mX + ")").c_str());
    dynamic_cast<RooSimultaneous*>(m_wk->pdf("simultaneous_CBGA"))->addPdf(*m_wk->pdf(("CBGA_Xhh_m" + mX).c_str()), mX.c_str());
  }

  RooCategory* SignalModel::mass_points() const {
    return m_wk->cat("mass_points");
  }

  void SignalModel::fit(RooDataSet& data)
  {
    m_data = &data;
    for (const auto& model : m_models) {
      MSG_INFO("Fitting simultaneous " << (model == "EGE" ? "ExpGausExp" : model == "CBGA" ? "CB+GA" : model) << " PDF to input events");
      m_wk->pdf(("simultaneous_" + model).c_str())->fitTo(data, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
      m_wk->pdf(("simultaneous_" + model).c_str())->fitTo(data, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(1));
      // Print the model parameters
      for (auto metaparameter : m_model_metaparameters.at(model)) {
        MSG_INFO("... " << std::left << std::setw(20) << (metaparameter + ":    ") << m_wk->var(metaparameter.c_str())->getVal());
      }
    }
  }

  void SignalModel::write(const std::string& output_file_name)
  {
    // Write EGE to output
    for (auto metaparameter : m_model_metaparameters.at("EGE")) {
      m_wk->var(metaparameter.c_str())->setConstant();
    }
    m_wk->factory("expr::EGE_mu('EGE_mu_p0 + EGE_mu_p1 * mass_resonance', EGE_mu_p0, EGE_mu_p1, mass_resonance)");
    m_wk->factory(("expr::EGE_sigma('TMath::Min(TMath::Max(EGE_sigma_p0 + EGE_sigma_p1 * mass_resonance, " + std::to_string(m_mass_category == "low" ? 3.0 : 6.0 ) + "), " + std::to_string(m_mass_category == "low" ? 10.0 : 18.0 ) + ")', EGE_sigma_p0, EGE_sigma_p1, mass_resonance)").c_str());
    m_wk->factory(("expr::EGE_kL('TMath::Min(TMath::Max(EGE_kL_p0 + EGE_kL_p1 * mass_resonance + EGE_kL_p2 * mass_resonance * mass_resonance, 0.1), " + std::to_string(m_mass_category == "low" ? 3.0 : 1.5) + ")', EGE_kL_p0, EGE_kL_p1, EGE_kL_p2, mass_resonance)").c_str());
    m_wk->factory(("expr::EGE_kH('TMath::Min(TMath::Max(EGE_kH_p0 + EGE_kH_p1 * mass_resonance + EGE_kH_p2 * mass_resonance * mass_resonance, 0.1), " + std::to_string(m_mass_category == "low" ? 3.0 : 1.0) + ")', EGE_kH_p0, EGE_kH_p1, EGE_kH_p2, mass_resonance)").c_str());
    m_wk->factory("ExpGausExpPDF::signal_PDF(mass, EGE_mu, EGE_sigma, EGE_kL, EGE_kH)");
    MSG_INFO("Preparing to write workspace to " << output_file_name);
    m_wk->writeToFile(output_file_name.c_str(), false);
  }

  void SignalModel::plot()
  {
    MSG_INFO("Plotting \033[1m" << m_mass_category << " mass " << m_tag_category << "-tag\033[0m");
    PlotStyle::EnsureAtlasStyle();

    // Loop over signal models
    for (auto model : m_models) {
      // Plot all fits on one canvas
      TCanvas canvas("canvas", "canvas", 600, 600);
      RooPlot* frame = m_wk->var("mass")->frame(RooFit::Range((m_mass_category == "low" ? 245 : 335), (m_mass_category == "low" ? 485 : 1140)));
      for (unsigned int idx = 0; idx < PlotStyle::resonance_masses(m_mass_category).size(); ++idx) {
        std::string mX(std::to_string(PlotStyle::resonance_masses(m_mass_category).at(idx)));
        int colour(PlotStyle::colours().at(idx));
        m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::MarkerColor(colour));
        m_wk->pdf(("simultaneous_" + model).c_str())->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(colour));
      }
      frame->Draw();
      frame->SetMinimum(0);
      canvas.Print(("plots/signal_model_" + model + "_" + m_mass_category + "Mass_" + m_tag_category + "tag_overall.pdf").c_str());
      MSG_INFO("Created plots/signal_model_" + model + "_" << m_mass_category << "Mass_" << m_tag_category << "tag_overall.pdf");

      // Plot each fit individually
      for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) {
        // Construct useful variables
        std::string mX(std::to_string(resonance_mass));
        double mass_low(0.9 * resonance_mass), mass_high(1.1 * resonance_mass);
        int nBins(int(mass_high - mass_low));
        while (nBins > 140) { nBins = int(nBins / 2.0); }
        // Fit with individual ExpGausExp
        RooDataSet* data_slice = dynamic_cast<RooDataSet*>(m_data->reduce(RooFit::Cut(("mass_points==mass_points::" + mX).c_str())));
        for (auto& parameter : m_model_parameters.at(model)) {
          m_wk->var(("individual_" + parameter).c_str())->setVal(m_wk->function((parameter + "_Xhh_m" + mX).c_str())->getVal());
        }
        m_wk->pdf(("individual_" + model).c_str())->fitTo(*data_slice, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::PrintLevel(-1));
        for (auto& parameter : m_model_parameters.at(model)) {
          MSG_INFO("... individual_" << std::left << std::setw(20) << (parameter + ":    ") << m_wk->var(("individual_" + parameter).c_str())->getVal());
          MSG_WARNING(parameter << "[" << m_tag_category << "][\"" << m_mass_category << "\"].append(" << m_wk->var(("individual_" + parameter).c_str())->getVal() << ")");
        }
        // Construct frames
        RooPlot* frame = m_wk->var("mass")->frame(mass_low, mass_high, nBins);
        RooPlot* frame_ratio = m_wk->var("mass")->frame(mass_low, mass_high, nBins);
        m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(kBlack));
        m_wk->pdf(("simultaneous_" + model).c_str())->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(kRed));
        RooHist* pull_hist_full = frame->pullHist(); pull_hist_full->SetLineColor(kRed); pull_hist_full->SetMarkerColor(kRed);
        data_slice->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::Invisible());
        m_wk->pdf(("individual_" + model).c_str())->plotOn(frame, RooFit::LineColor(kBlue));
        RooHist* pull_hist_simple = frame->pullHist(); pull_hist_simple->SetLineColor(kBlue); pull_hist_simple->SetMarkerColor(kBlue);
        if (model == "CBGA") { // Plot Leo's parameterisation if relevant
          m_wk->var("mass_resonance")->setVal(resonance_mass);
          m_wk->pdf("Leo_CBGA")->plotOn(frame, RooFit::LineColor(kGreen + 3));
        }
        // Construct a histogram with the pulls of the data w.r.t the curve
        frame_ratio->addPlotable(pull_hist_full, "P");
        frame_ratio->addPlotable(pull_hist_simple, "P");
        // Setup canvas
        TCanvas canvas_individual("canvas_individual", "canvas_individual", 600, 600);
        gPad->SetBottomMargin(0.35);
        frame->SetMinimum(1e-10);
        frame->GetXaxis()->SetLabelOffset(-100);
        frame->Draw();
        // Draw text
        TLatex textBox; textBox.SetNDC(); textBox.SetTextFont(42); textBox.SetTextSize(0.02);
        if (model == "EGE") {
          textBox.SetTextColor(kBlue);
          textBox.DrawLatex(0.19, 0.90, ("Single " + model + " fit").c_str());
          textBox.DrawLatex(0.19, 0.86, ("#mu_{EGE} = " + PlotStyle::to_string(m_wk->var("individual_EGE_mu")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.82, ("#sigma_{EGE} = " + PlotStyle::to_string(m_wk->var("individual_EGE_sigma")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.78, ("k_{L,EGE} = " + PlotStyle::to_string(m_wk->var("individual_EGE_kL")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.74, ("k_{H,EGE} = " + PlotStyle::to_string(m_wk->var("individual_EGE_kH")->getVal(), 2)).c_str());
          textBox.SetTextColor(kRed);
          textBox.DrawLatex(0.74, 0.90, ("Simultaneous " + model + " fit").c_str());
          textBox.DrawLatex(0.74, 0.86, ("#mu_{EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_mu_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.82, ("#sigma_{EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_sigma_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.78, ("k_{L,EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_kL_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.74, ("k_{H,EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_kH_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
        } else if (model == "CBGA") {
          textBox.SetTextColor(kBlue);
          textBox.DrawLatex(0.19, 0.90, ("Single " + model + " fit").c_str());
          textBox.DrawLatex(0.19, 0.86, ("#alpha_{CB} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_alpha_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.82, ("f_{CB} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_f_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.78, ("#mu_{CB} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_mu_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.74, ("#mu_{GA} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_mu_GA")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.70, ("n_{CB} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_n_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.66, ("#sigma_{CB} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_sigma_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.62, ("#sigma_{GA} = " + PlotStyle::to_string(m_wk->var("individual_CBGA_sigma_GA")->getVal(), 2)).c_str());
          textBox.SetTextColor(kRed);
          textBox.DrawLatex(0.74, 0.90, ("Simultaneous " + model + " fit").c_str());
          textBox.DrawLatex(0.74, 0.86, ("#alpha_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_alpha_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.82, ("f_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_f_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.78, ("#mu_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_mu_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.74, ("#mu_{GA} = " + PlotStyle::to_string(m_wk->function(("CBGA_mu_GA_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.70, ("n_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_n_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.66, ("#sigma_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_sigma_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.62, ("#sigma_{GA} = " + PlotStyle::to_string(m_wk->function(("CBGA_sigma_GA_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
        }
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
        canvas_individual.Print(("plots/" + m_mass_category + "Mass_" + m_tag_category + "tag/signal_model_" + model + "_" + m_mass_category + "Mass_" + m_tag_category + "tag_mX_" + mX + ".pdf").c_str());
        MSG_INFO("Created plots/" << m_mass_category << "Mass_" << m_tag_category << "tag/signal_model_" << model << "_" << m_mass_category << "Mass_" << m_tag_category << "tag_mX_" << mX << ".pdf");
      }
    }
  }

  std::string SignalModel::parameterised_CBGA(const int& parameter) const
  {
    std::vector<double> parameters;
    if (m_mass_category == "low") {
      if (m_tag_category == "1") {
        std::vector<double> _p({95.95, -0.6169, 0.0009952, 0.1473, -0.001951, 4.71, 0.003299, 0.0, 4.801, 0.003105, 0.0, -9.584, 0.06117, -3.98e-05, 248.8, -1.708, 0.003045});
        parameters.insert(parameters.end(), _p.begin(), _p.end());
      } else if (m_tag_category == "2") {
        std::vector<double> _p({36.15, -0.2098, 0.0003138, 0.1913, -0.001764, 4.717, 0.003271, 0.0, 4.7, 0.003387, 0.0, -9.18, 0.05983, -4.638e-05, -105.7, 0.7039, -0.001057});
        parameters.insert(parameters.end(), _p.begin(), _p.end());
      }
    } else if (m_mass_category == "high") {
      if (m_tag_category == "1") {
        std::vector<double> _p({-4.974, 0.01985, -1.36e-05, -0.5062, 6.076e-05, 21.71, 0.935, 3.826e-05, 64.61, 0.8174, 0.0001511, 10.39, -0.01581, 2.328e-05, -21.5, 0.1165, -5.918e-05});
        parameters.insert(parameters.end(), _p.begin(), _p.end());
      } else if (m_tag_category == "2") {
        std::vector<double> _p({2.053, -0.003402, 4.28e-06, -0.4215, -0.000229, 4.241, 0.9892, 2.646e-06, 13.43, 0.9817, 1.394e-05, 3.707, 0.002918, 9.237e-06, -3.17, 0.04017, 3.898e-06});
        parameters.insert(parameters.end(), _p.begin(), _p.end());
      }
    }
    if (parameters.empty()) { return "0.0"; }
    return std::to_string(parameters.at(parameter));
  }
}
