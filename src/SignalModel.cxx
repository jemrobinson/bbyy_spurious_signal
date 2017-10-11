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
  SignalModel::SignalModel(const std::string& mass_category, const std::string& tag_category, const std::vector<std::string>& model_names)
    : m_mass_category(mass_category)
    , m_tag_category(tag_category)
    , m_wk(0)
    , m_data(0)
    , m_models(model_names)
    , m_model_parameters({{"CBGA", {"CBGA_alpha_CB", "CBGA_f_CB", "CBGA_mu_CB", "CBGA_mu_GA", "CBGA_n_CB", "CBGA_sigma_CB", "CBGA_sigma_GA"}},
                          {"EGE",  {"EGE_mu", "EGE_sigma_k", "EGE_kL", "EGE_kH"}}})
    , m_model_metaparameters({{"CBGA", {"CBGA_alpha_CB_p0", "CBGA_alpha_CB_p1", "CBGA_f_CB_p0", "CBGA_f_CB_p1", "CBGA_mu_CB_p0", "CBGA_mu_CB_p1", "CBGA_mu_GA_p0", "CBGA_mu_GA_p1", "CBGA_n_CB", "CBGA_sigma_CB_p0", "CBGA_sigma_CB_p1", "CBGA_sigma_CB_p2", "CBGA_sigma_GA_p0", "CBGA_sigma_GA_p1", "CBGA_sigma_GA_p2"}},
                              {"EGE",  {"EGE_mu_p0", "EGE_mu_p1", "EGE_sigma_p0", "EGE_sigma_p1", "EGE_sigma_p2", "EGE_kL_p0", "EGE_kL_p1", "EGE_kL_p2", "EGE_kH_p0", "EGE_kH_p1", "EGE_kH_p2"}}})
  {
    MSG_INFO("Configuring signal model fits for " << m_models.size() << " models.");
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
      m_wk->factory(("EGE_mu_p0[" + std::to_string(get_initial_value("EGE_mu_p0")) + ", -10, 10]").c_str());
      m_wk->factory(("EGE_mu_p1[" + std::to_string(get_initial_value("EGE_mu_p1")) + ", 0.985, 1.015]").c_str());
      m_wk->factory(("EGE_sigma_p0[" + std::to_string(get_initial_value("EGE_sigma_p0")) + ", -20, 20]").c_str());//-20, 20]").c_str());
      m_wk->factory(("EGE_sigma_p1[" + std::to_string(get_initial_value("EGE_sigma_p1")) + ", -5, 5]").c_str());//-0.1, 0.1]").c_str());
      m_wk->factory(("EGE_sigma_p2[" + std::to_string(get_initial_value("EGE_sigma_p2")) + ", -5, 5]").c_str());//1.0, 8.0]").c_str());
      m_wk->factory(("EGE_kL_p0[" + std::to_string(get_initial_value("EGE_kL_p0")) + ", -100, 100]").c_str());//-100, 300]").c_str());
      m_wk->factory(("EGE_kL_p1[" + std::to_string(get_initial_value("EGE_kL_p1")) + ", -100, 100]").c_str());//-1, 1]").c_str());
      m_wk->factory(("EGE_kL_p2[" + std::to_string(get_initial_value("EGE_kL_p2")) + ", -100, 100]").c_str());//-20000, 20000]").c_str());
      m_wk->factory(("EGE_kH_p0[" + std::to_string(get_initial_value("EGE_kH_p0")) + ", -10, 10]").c_str());//-5, 10]").c_str());
      m_wk->factory(("EGE_kH_p1[" + std::to_string(get_initial_value("EGE_kH_p1")) + ", -10, 10]").c_str());//-0.1, 0.1]").c_str());
      m_wk->factory(("EGE_kH_p2[" + std::to_string(get_initial_value("EGE_kH_p2")) + ", -10, 10]").c_str());//-0.001, 0.001]").c_str());
    }

    // Simultaneous Crystal Ball + Gaussian
    if (std::find(m_models.begin(), m_models.end(), "CBGA") != m_models.end()) {
      m_wk->factory(("CBGA_alpha_CB_p0[" + std::to_string(get_initial_value("CBGA_alpha_CB_p0")) + ", 0, 10]").c_str());
      m_wk->factory(("CBGA_alpha_CB_p1[" + std::to_string(get_initial_value("CBGA_alpha_CB_p1")) + ", -0.1, 0.2]").c_str());
      m_wk->factory(("CBGA_f_CB_p0[" + std::to_string(get_initial_value("CBGA_f_CB_p0")) + ", 0.3, 1.8]").c_str());
      m_wk->factory(("CBGA_f_CB_p1[" + std::to_string(get_initial_value("CBGA_f_CB_p1")) + ", -0.1, 0.1]").c_str());
      m_wk->factory(("CBGA_mu_CB_p0[" + std::to_string(get_initial_value("CBGA_mu_CB_p0")) + ", 0, 30]").c_str());
      m_wk->factory(("CBGA_mu_CB_p1[" + std::to_string(get_initial_value("CBGA_mu_CB_p1")) + ", 0.95, 1.05]").c_str());
      m_wk->factory(("CBGA_mu_GA_p0[" + std::to_string(get_initial_value("CBGA_mu_GA_p0")) + ", 5, 40]").c_str());
      m_wk->factory(("CBGA_mu_GA_p1[" + std::to_string(get_initial_value("CBGA_mu_GA_p1")) + ", 0.9, 1.1]").c_str());
      m_wk->factory(("CBGA_n_CB[" + std::to_string(get_initial_value("CBGA_n_CB")) + ", 5, 200]").c_str());
      m_wk->factory(("CBGA_sigma_CB_p0[" + std::to_string(get_initial_value("CBGA_sigma_CB_p0")) + ", -70, 0]").c_str());
      m_wk->factory(("CBGA_sigma_CB_p1[" + std::to_string(get_initial_value("CBGA_sigma_CB_p1")) + ", 0, 1.0]").c_str());
      m_wk->factory(("CBGA_sigma_CB_p2[" + std::to_string(get_initial_value("CBGA_sigma_CB_p2")) + ", -0.1, 0.1]").c_str());
      m_wk->factory(("CBGA_sigma_GA_p0[" + std::to_string(get_initial_value("CBGA_sigma_GA_p0")) + ", -300, 0]").c_str());
      m_wk->factory(("CBGA_sigma_GA_p1[" + std::to_string(get_initial_value("CBGA_sigma_GA_p1")) + ", 0, 2.5]").c_str());
      m_wk->factory(("CBGA_sigma_GA_p2[" + std::to_string(get_initial_value("CBGA_sigma_GA_p2")) + ", -0.1, 0.2]").c_str());
    }
  }

  void SignalModel::build_simultaneous_PDF(RooRealVar& mass)
  {
    m_wk->import(mass);
    // Individual ExpGausExp (and import the code into the workspace)
    if (std::find(m_models.begin(), m_models.end(), "EGE") != m_models.end()) {
      m_wk->factory("individual_EGE_mu[260, 0, 2000]");
      m_wk->factory("individual_EGE_sigma_k[0.01, 0.01, 0.03]");
      m_wk->factory("prod::individual_EGE_sigma(individual_EGE_sigma_k, individual_EGE_mu)");
      m_wk->factory("individual_EGE_kL[2, 0.0, 20]");
      m_wk->factory("individual_EGE_kH[2, 0.0, 20]");
      m_wk->factory("ExpGausExpPDF::individual_EGE(mass, individual_EGE_mu, individual_EGE_sigma, individual_EGE_kL, individual_EGE_kH)");
      m_wk->addClassDeclImportDir("include/");
      m_wk->importClassCode("ExpGausExpPDF*");
      // Build a simultaneous ExpGausExp
      m_wk->factory("Simultaneous::simultaneous_EGE(mass_points)");
    }
    // Individual Crystal Ball + Gaussian
    if (std::find(m_models.begin(), m_models.end(), "CBGA") != m_models.end()) {
      m_wk->factory("individual_CBGA_alpha_CB[0.1, 0.0, 20.0]");
      m_wk->factory("individual_CBGA_f_CB[0.15, 0.0, 1.0]");
      m_wk->factory("individual_CBGA_mu_CB[260, 0, 2000]");
      m_wk->factory("individual_CBGA_mu_GA[260, 0, 2000]");
      m_wk->factory("individual_CBGA_n_CB[10, 0, 200]");
      m_wk->factory("individual_CBGA_sigma_CB[4, 0.0, 40]");
      m_wk->factory("individual_CBGA_sigma_GA[4, 0.0, 40]");
      m_wk->factory("CBShape::individual_CBGA_CB(mass, individual_CBGA_mu_CB, individual_CBGA_sigma_CB, individual_CBGA_alpha_CB, individual_CBGA_n_CB)");
      m_wk->factory("Gaussian:individual_CBGA_GA(mass,individual_CBGA_mu_GA, individual_CBGA_sigma_GA)");
      m_wk->factory("SUM::individual_CBGA(individual_CBGA_f_CB * individual_CBGA_CB, individual_CBGA_GA)");
      // Build a simultaneous Crystal Ball + Gaussian
      m_wk->factory("Simultaneous::simultaneous_CBGA(mass_points)");
    }
    for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category) ) {
      add_mass_point(resonance_mass);
    }
  }

  void SignalModel::add_mass_point(const int& resonance_mass) {
    std::string mX(std::to_string(resonance_mass));
    // Add ExpGausExp
    if (std::find(m_models.begin(), m_models.end(), "EGE") != m_models.end()) {
      m_wk->factory(("expr::Xhh_m" + mX + "_EGE_mu('EGE_mu_p0 + EGE_mu_p1 * " + mX + "', EGE_mu_p0, EGE_mu_p1)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_EGE_sigma('TMath::Exp(EGE_sigma_p0 + EGE_sigma_p1 * " + mX + " / 100.0 + EGE_sigma_p2 * " + mX + " * " + mX + " / 10000.0)', EGE_sigma_p0, EGE_sigma_p1, EGE_sigma_p2)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_EGE_kL('EGE_kL_p0 + EGE_kL_p1 / ((" + mX + " / 100.0) - 2.5) + EGE_kL_p2 / (((" + mX + " / 100.0) - 2.5) * ((" + mX + " / 100.0) - 2.5))', EGE_kL_p0, EGE_kL_p1, EGE_kL_p2)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_EGE_kH('TMath::Exp(EGE_kH_p0 + EGE_kH_p1 * " + mX + " / 100.0 + EGE_kH_p2 * " + mX + " * " + mX + " / 10000.0)', EGE_kH_p0, EGE_kH_p1, EGE_kH_p2)").c_str());
      m_wk->factory(("ExpGausExpPDF::Xhh_m" + mX + "_EGE(mass, Xhh_m" + mX + "_EGE_mu, Xhh_m" + mX + "_EGE_sigma, Xhh_m" + mX + "_EGE_kL, Xhh_m" + mX + "_EGE_kH)").c_str());
      dynamic_cast<RooSimultaneous*>(m_wk->pdf("simultaneous_EGE"))->addPdf(*m_wk->pdf(("Xhh_m" + mX + "_EGE").c_str()), mX.c_str());
    }
    // Add Crystal Ball + Gaussian
    if (std::find(m_models.begin(), m_models.end(), "CBGA") != m_models.end()) {
      m_wk->factory(("expr::Xhh_m" + mX + "_CBGA_alpha_CB('CBGA_alpha_CB_p0 + CBGA_alpha_CB_p1 * " + mX + "', CBGA_alpha_CB_p0, CBGA_alpha_CB_p1)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_CBGA_f_CB('CBGA_f_CB_p0 + CBGA_f_CB_p1 * " + mX + "', CBGA_f_CB_p0, CBGA_f_CB_p1)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_CBGA_mu_CB('CBGA_mu_CB_p0 + CBGA_mu_CB_p1 * " + mX + "', CBGA_mu_CB_p0, CBGA_mu_CB_p1)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_CBGA_mu_GA('CBGA_mu_GA_p0 + CBGA_mu_GA_p1 * " + mX + "', CBGA_mu_GA_p0, CBGA_mu_GA_p1)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_CBGA_sigma_CB('CBGA_sigma_CB_p0 + CBGA_sigma_CB_p1 * " + mX + " + CBGA_sigma_CB_p2 * " + mX + " * " + mX + "', CBGA_sigma_CB_p0, CBGA_sigma_CB_p1, CBGA_sigma_CB_p2)").c_str());
      m_wk->factory(("expr::Xhh_m" + mX + "_CBGA_sigma_GA('CBGA_sigma_GA_p0 + CBGA_sigma_GA_p1 * " + mX + " + CBGA_sigma_GA_p2 * " + mX + " * " + mX + "', CBGA_sigma_GA_p0, CBGA_sigma_GA_p1, CBGA_sigma_GA_p2)").c_str());
      m_wk->factory(("CBShape::Xhh_m" + mX + "_CBGA_CB(mass, Xhh_m" + mX + "_CBGA_mu_CB, Xhh_m" + mX + "_CBGA_sigma_CB, Xhh_m" + mX + "_CBGA_alpha_CB, CBGA_n_CB)").c_str());
      m_wk->factory(("Gaussian::Xhh_m" + mX + "_CBGA_GA(mass, CBGA_mu_GA_Xhh_m" + mX + ", CBGA_sigma_GA_Xhh_m" + mX + ")").c_str());
      m_wk->factory(("SUM::Xhh_m" + mX + "_CBGA(Xhh_m" + mX + "_CBGA_f_CB * Xhh_m" + mX + "_CBGA_CB, Xhh_m" + mX + "_CBGA_GA)").c_str());
      dynamic_cast<RooSimultaneous*>(m_wk->pdf("simultaneous_CBGA"))->addPdf(*m_wk->pdf(("Xhh_m" + mX + "_CBGA").c_str()), mX.c_str());
    }
  }

  RooCategory* SignalModel::mass_points() const {
    return m_wk->cat("mass_points");
  }

  void SignalModel::fit(RooDataSet& data)
  {
    m_data = &data;
    for (const auto& model : m_models) {
      MSG_INFO("Fitting simultaneous " << (model == "EGE" ? "ExpGausExp" : model == "CBGA" ? "Crystal Ball + Gaussian" : model) << " PDF to input events");
      m_wk->pdf(("simultaneous_" + model).c_str())->fitTo(data, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "minimize"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(-1));
      m_wk->pdf(("simultaneous_" + model).c_str())->fitTo(data, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(-1));
      // Print the model parameters
      for (auto metaparameter : m_model_metaparameters.at(model)) {
        MSG_INFO("... " << std::left << std::setw(24) << (metaparameter + ":    ") << m_wk->var(metaparameter.c_str())->getVal());
      }
    }
  }

  void SignalModel::write(const std::string& output_file_name)
  {
    MSG_INFO("Preparing to write workspace to " << output_file_name);
    for (const auto& model : m_models) {
      // Set meta-parameters constant
      for (auto metaparameter : m_model_metaparameters.at(model)) {
        m_wk->var(metaparameter.c_str())->setConstant();
      }
      // Write EGE to output
      if (model == "EGE") {
        m_wk->factory("expr::EGE_mu('EGE_mu_p0 + EGE_mu_p1 * mass_resonance', EGE_mu_p0, EGE_mu_p1, mass_resonance)");
        m_wk->factory("expr::EGE_sigma('TMath::Exp(EGE_sigma_p0 + EGE_sigma_p1 * mass_resonance / 100.0 + EGE_sigma_p2 * mass_resonance * mass_resonance / 10000.0)', EGE_sigma_p0, EGE_sigma_p1, EGE_sigma_p2, mass_resonance)");
        m_wk->factory("expr::EGE_kL('EGE_kL_p0 + EGE_kL_p1 / ((mass_resonance / 100.0) - 2.5) + EGE_kL_p2 / (((mass_resonance / 100.0) - 2.5) * ((mass_resonance / 100.0) - 2.5))', EGE_kL_p0, EGE_kL_p1, EGE_kL_p2, mass_resonance)");
        m_wk->factory("expr::EGE_kH('TMath::Exp(EGE_kH_p0 + EGE_kH_p1 * mass_resonance / 100.0 + EGE_kH_p2 * mass_resonance * mass_resonance / 10000.0)', EGE_kH_p0, EGE_kH_p1, EGE_kH_p2, mass_resonance)");
        m_wk->factory("ExpGausExpPDF::signal_PDF(mass, EGE_mu, EGE_sigma, EGE_kL, EGE_kH)");
        m_wk->writeToFile(output_file_name.c_str(), false);
      }
    }
  }

  void SignalModel::plot()
  {
    MSG_INFO("Plotting \033[1m" << m_mass_category << " mass " << m_tag_category << "-tag\033[0m");
    PlotStyle::EnsureAtlasStyle();

    // Loop over signal models
    for (auto model : m_models) {
      // Plot all fits on one canvas
      TCanvas canvas("canvas", "canvas", 600, 600);
      RooPlot* frame = m_wk->var("mass")->frame(RooFit::Range(PlotStyle::mass_range(m_mass_category).first, PlotStyle::mass_range(m_mass_category).second));
      for (unsigned int idx = 0; idx < PlotStyle::resonance_masses(m_mass_category).size(); ++idx) {
        std::string mX(std::to_string(PlotStyle::resonance_masses(m_mass_category).at(idx)));
        m_wk->var("mass_resonance")->setVal(PlotStyle::resonance_masses(m_mass_category).at(idx));
        int colour(PlotStyle::colours().at(idx));
        m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::MarkerColor(colour));
        m_wk->pdf(("simultaneous_" + model).c_str())->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(colour));
      }
      frame->Draw();
      frame->SetMinimum(0);
      canvas.Print(("plots/signal_model/overall/overall_signal_model_" + model + "_" + m_mass_category + "Mass_" + m_tag_category + "tag_overall.pdf").c_str());
      MSG_INFO("Created plots/signal_model_" + model + "_" << m_mass_category << "Mass_" << m_tag_category << "tag_overall.pdf");

      // Plot each fit individually
      for (auto resonance_mass : PlotStyle::resonance_masses(m_mass_category)) {
        // Construct useful variables
        std::string mX(std::to_string(resonance_mass));
        m_wk->var("mass_resonance")->setVal(resonance_mass);
        double mass_low(0.9 * resonance_mass), mass_high(1.1 * resonance_mass);
        int nBins(int(mass_high - mass_low));
        while (nBins > 140) { nBins = int(nBins / 2.0); }
        // Construct frames
        RooPlot* frame = m_wk->var("mass")->frame(mass_low, mass_high, nBins);
        RooPlot* frame_ratio = m_wk->var("mass")->frame(mass_low, mass_high, nBins);
        // Fit with individual PDF
        RooDataSet* data_slice = dynamic_cast<RooDataSet*>(m_data->reduce(RooFit::Cut(("mass_points==mass_points::" + mX).c_str())));
        for (auto& parameter : m_model_parameters.at(model)) {
          if (parameter.find("mu") != std::string::npos) { m_wk->var(("individual_" + parameter).c_str())->setVal(resonance_mass); }
        }
        m_wk->pdf(("individual_" + model).c_str())->fitTo(*data_slice, RooFit::SumW2Error(true), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::PrintLevel(-1));
        for (auto& parameter : m_model_parameters.at(model)) {
          MSG_INFO("... individual_" << std::left << std::setw(24) << (parameter + ":    ") << m_wk->var(("individual_" + parameter).c_str())->getVal() << " +/- " << m_wk->var(("individual_" + parameter).c_str())->getError());
          MSG_DEBUG(parameter << "[" << m_tag_category << "][\"" << m_mass_category << "\"].append(" << m_wk->var(("individual_" + parameter).c_str())->getVal() << ")");
        }
        // Plot data and simultaneous EGE
        m_data->plotOn(frame, RooFit::Cut(("mass_points==mass_points::" + mX).c_str()), RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(kBlack));
        TGraph* g_data = dynamic_cast<TGraph*>(frame->getObject(frame->numItems() - 1)); int ndof(g_data->GetN());
        m_wk->pdf(("simultaneous_" + model).c_str())->plotOn(frame, RooFit::Slice(*m_wk->cat("mass_points"), mX.c_str()), RooFit::ProjWData(*m_wk->cat("mass_points"), *m_data), RooFit::LineColor(kRed));
        double chiSquared(frame->chiSquare() * ndof); int nFitParams(m_wk->pdf(("simultaneous_" + model).c_str())->getParameters(m_data)->getSize());
        RooHist* pull_hist_parameterised = frame->pullHist(); pull_hist_parameterised->SetLineColor(kRed); pull_hist_parameterised->SetMarkerColor(kRed);
        // Plot individual
        data_slice->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::Invisible());
        m_wk->pdf(("individual_" + model).c_str())->plotOn(frame, RooFit::LineColor(kBlue));
        double individual_chiSquared(frame->chiSquare() * ndof); int individual_nFitParams(m_wk->pdf(("individual_" + model).c_str())->getParameters(m_data)->getSize());
        RooHist* pull_hist_individual = frame->pullHist(); pull_hist_individual->SetLineColor(kBlue); pull_hist_individual->SetMarkerColor(kBlue);
        // Add pulls to ratio frame
        frame_ratio->addPlotable(pull_hist_parameterised, "P");
        frame_ratio->addPlotable(pull_hist_individual, "P");
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
          textBox.DrawLatex(0.19, 0.82, ("#sigma_{EGE} = " + PlotStyle::to_string(m_wk->function("individual_EGE_sigma")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.78, ("k_{L,EGE} = " + PlotStyle::to_string(m_wk->var("individual_EGE_kL")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.74, ("k_{H,EGE} = " + PlotStyle::to_string(m_wk->var("individual_EGE_kH")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.19, 0.70, ("#chi^{2} / ndof = " + PlotStyle::to_string(individual_chiSquared, 2) + "/" + PlotStyle::to_string(ndof - individual_nFitParams, 0)).c_str());
          textBox.SetTextColor(kRed);
          textBox.DrawLatex(0.74, 0.90, ("Simultaneous " + model + " fit").c_str());
          textBox.DrawLatex(0.74, 0.86, ("#mu_{EGE} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_EGE_mu").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.82, ("#sigma_{EGE} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_EGE_sigma").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.78, ("k_{L,EGE} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_EGE_kL").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.74, ("k_{H,EGE} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_EGE_kH").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.70, ("#chi^{2} / ndof = " + PlotStyle::to_string(chiSquared, 2) + "/" + PlotStyle::to_string(ndof - nFitParams, 0)).c_str());
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
          textBox.DrawLatex(0.19, 0.58, ("#chi^{2} / ndof = " + PlotStyle::to_string(individual_chiSquared, 2) + "/" + PlotStyle::to_string(ndof - individual_nFitParams, 0)).c_str());
          textBox.SetTextColor(kRed);
          textBox.DrawLatex(0.74, 0.90, ("Simultaneous " + model + " fit").c_str());
          textBox.DrawLatex(0.74, 0.86, ("#alpha_{CB} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_CBGA_alpha_CB").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.82, ("f_{CB} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_CBGA_f_CB").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.78, ("#mu_{CB} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_CBGA_mu_CB").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.74, ("#mu_{GA} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_CBGA_mu_GA").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.70, ("n_{CB} = " + PlotStyle::to_string(m_wk->var("CBGA_n_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.66, ("#sigma_{CB} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_CBGA_sigma_CB").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.62, ("#sigma_{GA} = " + PlotStyle::to_string(m_wk->function(("Xhh_m" + mX + "_CBGA_sigma_GA").c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.58, ("#chi^{2} / ndof = " + PlotStyle::to_string(chiSquared, 2) + "/" + PlotStyle::to_string(ndof - nFitParams, 0)).c_str());
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
        canvas_individual.Print(("plots/signal_model/" + model + "/signal_model_" + model + "_" + m_mass_category + "Mass_" + m_tag_category + "tag_mX_" + mX + ".pdf").c_str());
        MSG_INFO("Created plots/" << m_mass_category << "Mass_" << m_tag_category << "tag/signal_model_" << model << "_" << m_mass_category << "Mass_" << m_tag_category << "tag_mX_" << mX << ".pdf");
        // Print individual fit parameters
        for (auto& parameter : m_model_parameters.at(model)) {
          RooRealVar *wk_parameter = m_wk->var(("individual_" + parameter).c_str());
          (void)wk_parameter; // to avoid compiler warnings when running in non-debug mode
          MSG_DEBUG(m_mass_category << " " << m_tag_category << " " << parameter << " " << resonance_mass << " " << wk_parameter->getVal() << " " << wk_parameter->getError());
        }
      }
    }
  }

  void SignalModel::append(std::vector<double>& target, std::vector<double> source) {
    target.insert(target.end(), source.begin(), source.end());
  }

  double SignalModel::get_initial_value(const std::string& name) {
    if (m_mass_category == "low" && m_tag_category == "0") {
      if (name == "EGE_mu_p0") { return -1.33107855; }
      if (name == "EGE_mu_p1") { return 1.004928891; }
      if (name == "EGE_sigma_p0") { return -3.97166592; }
      if (name == "EGE_sigma_p1") { return 2.81727456; }
      if (name == "EGE_sigma_p2") { return -0.2979523; }
      if (name == "EGE_kL_p0") { return 1.17591096; }
      if (name == "EGE_kL_p1") { return -0.30864717; }
      if (name == "EGE_kL_p2") { return 0.03796169; }
      if (name == "EGE_kH_p0") { return -6.14867379; }
      if (name == "EGE_kH_p1") { return 2.98087018; }
      if (name == "EGE_kH_p2") { return -0.35195593;}
    } else if (m_mass_category == "high" && m_tag_category == "0") {
      if (name == "EGE_mu_p0") { return 2.6927747; }
      if (name == "EGE_mu_p1") { return 0.9953319763; }
      if (name == "EGE_sigma_p0") { return 0.35508604; }
      if (name == "EGE_sigma_p1") { return 0.30728317; }
      if (name == "EGE_sigma_p2") { return -0.0111862; }
      if (name == "EGE_kL_p0") { return 0.46598494; }
      if (name == "EGE_kL_p1") { return -0.1558814; }
      if (name == "EGE_kL_p2") { return 0.00359117; }
      if (name == "EGE_kH_p0") { return -1.11877485; }
      if (name == "EGE_kH_p1") { return 9.74486849e-04; }
      if (name == "EGE_kH_p2") { return 1.07415847e-03; }
    } if (m_mass_category == "low" && m_tag_category == "1") {
      if (name == "EGE_mu_p0") { return -2.37223907; }
      if (name == "EGE_mu_p1") { return 1.0096604181; }
      if (name == "EGE_sigma_p0") { return -2.60476481; }
      if (name == "EGE_sigma_p1") { return 2.02458664; }
      if (name == "EGE_sigma_p2") { return -0.19813873; }
      if (name == "EGE_kL_p0") { return 1.08152049; }
      if (name == "EGE_kL_p1") { return -0.33258371; }
      if (name == "EGE_kL_p2") { return 0.0618694; }
      if (name == "EGE_kH_p0") { return -1.36019444; }
      if (name == "EGE_kH_p1") { return 0.38158853; }
      if (name == "EGE_kH_p2") { return -0.02376858; }
    } if (m_mass_category == "high" && m_tag_category == "1") {
      if (name == "EGE_mu_p0") { return 4.47617391; }
      if (name == "EGE_mu_p1") { return 0.992433913; }
      if (name == "EGE_sigma_p0") { return 0.78287279; }
      if (name == "EGE_sigma_p1") { return 0.29029052; }
      if (name == "EGE_sigma_p2") { return -0.00807065; }
      if (name == "EGE_kL_p0") { return 0.97026121; }
      if (name == "EGE_kL_p1") { return -0.23456268; }
      if (name == "EGE_kL_p2") { return -0.66998253; }
      if (name == "EGE_kH_p0") { return -1.03804543; }
      if (name == "EGE_kH_p1") { return 0.12539734; }
      if (name == "EGE_kH_p2") { return -0.00650364; }
    } if (m_mass_category == "low" && m_tag_category == "2") {
      if (name == "EGE_mu_p0") { return 0.1713232; }
      if (name == "EGE_mu_p1") { return 1.00000153; }
      if (name == "EGE_sigma_p0") { return -4.42029036; }
      if (name == "EGE_sigma_p1") { return 3.16102646; }
      if (name == "EGE_sigma_p2") { return  -0.3863813; }
      if (name == "EGE_kL_p0") { return 8.70715388e-01; }
      if (name == "EGE_kL_p1") { return 1.23627296e-01; }
      if (name == "EGE_kL_p2") { return -4.94317474e-04; }
      if (name == "EGE_kH_p0") { return -0.43097996; }
      if (name == "EGE_kH_p1") { return 0.30241119; }
      if (name == "EGE_kH_p2") { return -0.06352152; }
    } if (m_mass_category == "high" && m_tag_category == "2") {
      if (name == "EGE_mu_p0") { return 3.48818577; }
      if (name == "EGE_mu_p1") { return 0.9930448617; }
      if (name == "EGE_sigma_p0") { return 0.8865416; }
      if (name == "EGE_sigma_p1") { return 0.28844091; }
      if (name == "EGE_sigma_p2") { return -0.00945846; }
      if (name == "EGE_kL_p0") { return 1.07873558; }
      if (name == "EGE_kL_p1") { return -0.50070711; }
      if (name == "EGE_kL_p2") { return 0.75016217; }
      if (name == "EGE_kH_p0") { return 0.04114206; }
      if (name == "EGE_kH_p1") { return -0.11020904; }
      if (name == "EGE_kH_p2") { return 0.00727963; }
    }
    // CBGA
    if (m_tag_category == "0") {
      if (name == "CBGA_alpha_CB_p0") { return 0.000724086; }
      if (name == "CBGA_alpha_CB_p1") { return 0.00434892; }
      if (name == "CBGA_f_CB_p0") { return 0.400164; }
      if (name == "CBGA_f_CB_p1") { return 0.000639864; }
      if (name == "CBGA_mu_CB_p0") { return 1.18889; }
      if (name == "CBGA_mu_CB_p1") { return 0.995387; }
      if (name == "CBGA_mu_GA_p0") { return 21.5743; }
      if (name == "CBGA_mu_GA_p1") { return 0.94393; }
      if (name == "CBGA_n_CB") { return 50; }
      if (name == "CBGA_sigma_CB_p0") { return -22.0718; }
      if (name == "CBGA_sigma_CB_p1") { return 0.082453; }
      if (name == "CBGA_sigma_CB_p2") { return 5.61611e-05; }
      if (name == "CBGA_sigma_GA_p0") { return -152.392; }
      if (name == "CBGA_sigma_GA_p1") { return 1.02557; }
      if (name == "CBGA_sigma_GA_p2") { return -0.00155964; }
    } else if (m_tag_category == "1") {
      if (name == "CBGA_alpha_CB_p0") { return 3.123; }
      if (name == "CBGA_alpha_CB_p1") { return -0.0007643; }
      if (name == "CBGA_f_CB_p0") { return 0.9375; }
      if (name == "CBGA_f_CB_p1") { return -0.001607; }
      if (name == "CBGA_mu_CB_p0") { return 12.16; }
      if (name == "CBGA_mu_CB_p1") { return 0.9265; }
      if (name == "CBGA_mu_GA_p0") { return 27.04; }
      if (name == "CBGA_mu_GA_p1") { return 0.8281; }
      if (name == "CBGA_n_CB") { return  10; }
      if (name == "CBGA_sigma_CB_p0") { return -11.68; }
      if (name == "CBGA_sigma_CB_p1") { return 0.07075; }
      if (name == "CBGA_sigma_CB_p2") { return -0.00005255; }
      if (name == "CBGA_sigma_GA_p0") { return -172.8; }
      if (name == "CBGA_sigma_GA_p1") { return 1.105; }
      if (name == "CBGA_sigma_GA_p2") { return -0.001554; }
    } else if (m_tag_category == "2") {
      if (name == "CBGA_alpha_CB_p0") { return 3.6; }
      if (name == "CBGA_alpha_CB_p1") { return -0.002293; }
      if (name == "CBGA_f_CB_p0") { return 1.137; }
      if (name == "CBGA_f_CB_p1") { return -0.001675; }
      if (name == "CBGA_mu_CB_p0") { return 3.155; }
      if (name == "CBGA_mu_CB_p1") { return 0.984; }
      if (name == "CBGA_mu_GA_p0") { return 20.56; }
      if (name == "CBGA_mu_GA_p1") { return 0.8691; }
      if (name == "CBGA_n_CB") { return  10; }
      if (name == "CBGA_sigma_CB_p0") { return -47.54; }
      if (name == "CBGA_sigma_CB_p1") { return 0.3151; }
      if (name == "CBGA_sigma_CB_p2") { return -0.0004661; }
      if (name == "CBGA_sigma_GA_p0") { return -83.45; }
      if (name == "CBGA_sigma_GA_p1") { return 0.5328; }
      if (name == "CBGA_sigma_GA_p2") { return -0.0007261; }
    }
    MSG_WARNING("Failed to retrieve initial value of '" << name << "' for " << m_mass_category << " mass, " << m_tag_category << "-tag");
    return 0.0;
  }
}
