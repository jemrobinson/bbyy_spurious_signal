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
    , m_models({"EGE", "CBGA"})
    , m_model_parameters({{"CBGA", {"CBGA_alpha_CB", "CBGA_f_CB", "CBGA_mu_CB", "CBGA_mu_GA", "CBGA_n_CB", "CBGA_sigma_CB", "CBGA_sigma_GA"}},
                          {"EGE",  {"EGE_mu", "EGE_sigma_k", "EGE_kL", "EGE_kH"}}})
    , m_model_metaparameters({{"CBGA", {"CBGA_alpha_CB_p0", "CBGA_alpha_CB_p1", "CBGA_f_CB_p0", "CBGA_f_CB_p1", "CBGA_mu_CB_p0", "CBGA_mu_CB_p1", "CBGA_mu_GA_p0", "CBGA_mu_GA_p1", "CBGA_n_CB", "CBGA_sigma_CB_p0", "CBGA_sigma_CB_p1", "CBGA_sigma_CB_p2", "CBGA_sigma_GA_p0", "CBGA_sigma_GA_p1", "CBGA_sigma_GA_p2"}},
                              {"EGE",  {"EGE_mu_p0", "EGE_mu_p1", "EGE_sigma_p0", "EGE_sigma_p1", "EGE_kL_p0", "EGE_kL_p1", "EGE_kL_p2", "EGE_kH_p0", "EGE_kH_p1", "EGE_kH_p2"}}})
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
      m_wk->factory(("EGE_mu_p0[" + std::to_string(get_initial_value("EGE_mu_p0")) + ", -10, 10]").c_str());
      m_wk->factory(("EGE_mu_p1[" + std::to_string(get_initial_value("EGE_mu_p1")) + ", 0.985, 1.015]").c_str());
      m_wk->factory(("EGE_sigma_p0[" + std::to_string(get_initial_value("EGE_sigma_p0")) + ", -20, 20]").c_str());
      m_wk->factory(("EGE_sigma_p1[" + std::to_string(get_initial_value("EGE_sigma_p1")) + ", -0.1, 0.1]").c_str());
      m_wk->factory(("EGE_kL_p0[" + std::to_string(get_initial_value("EGE_kL_p0")) + ", -100, 300]").c_str());
      m_wk->factory(("EGE_kL_p1[" + std::to_string(get_initial_value("EGE_kL_p1")) + ", -1, 1]").c_str());
      m_wk->factory(("EGE_kL_p2[" + std::to_string(get_initial_value("EGE_kL_p2")) + ", -20000, 20000]").c_str());
      m_wk->factory(("EGE_kH_p0[" + std::to_string(get_initial_value("EGE_kH_p0")) + ", -5, 10]").c_str());
      m_wk->factory(("EGE_kH_p1[" + std::to_string(get_initial_value("EGE_kH_p1")) + ", -0.1, 0.1]").c_str());
      m_wk->factory(("EGE_kH_p2[" + std::to_string(get_initial_value("EGE_kH_p2")) + ", -0.001, 0.001]").c_str());
    }
    // Simultaneous Crystal Ball + Gaussian
    if (std::find(m_models.begin(), m_models.end(), "CBGA") != m_models.end()) {
      m_wk->factory(("CBGA_alpha_CB_p0[" + std::to_string(get_initial_value("CBGA_alpha_CB_p0")) + ", 0, 10]").c_str());
      m_wk->factory(("CBGA_alpha_CB_p1[" + std::to_string(get_initial_value("CBGA_alpha_CB_p1")) + ", -0.1, 0.2]").c_str());
      m_wk->factory(("CBGA_f_CB_p0[" + std::to_string(get_initial_value("CBGA_f_CB_p0")) + ", 0.4, 1.8]").c_str());
      m_wk->factory(("CBGA_f_CB_p1[" + std::to_string(get_initial_value("CBGA_f_CB_p1")) + ", -0.1, 0.1]").c_str());
      m_wk->factory(("CBGA_mu_CB_p0[" + std::to_string(get_initial_value("CBGA_mu_CB_p0")) + ", 0, 30]").c_str());
      m_wk->factory(("CBGA_mu_CB_p1[" + std::to_string(get_initial_value("CBGA_mu_CB_p1")) + ", 0.95, 1.05]").c_str());
      m_wk->factory(("CBGA_mu_GA_p0[" + std::to_string(get_initial_value("CBGA_mu_GA_p0")) + ", 5, 40]").c_str());
      m_wk->factory(("CBGA_mu_GA_p1[" + std::to_string(get_initial_value("CBGA_mu_GA_p1")) + ", 0.95, 1.05]").c_str());
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
    m_wk->factory("individual_EGE_mu[260, 0, 2000]");
    m_wk->factory("individual_EGE_sigma_k[0.01, 0.01, 0.03]");
    m_wk->factory("prod::individual_EGE_sigma(individual_EGE_sigma_k, individual_EGE_mu)");
    m_wk->factory("individual_EGE_kL[2, 0.0, 20]");
    m_wk->factory("individual_EGE_kH[2, 0.0, 20]");
    m_wk->factory("ExpGausExpPDF::individual_EGE(mass, individual_EGE_mu, individual_EGE_sigma, individual_EGE_kL, individual_EGE_kH)");
    m_wk->addClassDeclImportDir("include/");
    m_wk->importClassCode("ExpGausExpPDF*");
    // Individual Crystal Ball + Gaussian
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
    if (std::find(m_models.begin(), m_models.end(), "EGE") != m_models.end()) {
      m_wk->factory(("expr::EGE_mu_Xhh_m" + mX + "('EGE_mu_p0 + EGE_mu_p1 * " + mX + "', EGE_mu_p0, EGE_mu_p1)").c_str());
      m_wk->factory(("expr::EGE_sigma_Xhh_m" + mX + "('TMath::Min(TMath::Max(EGE_sigma_p0 + EGE_sigma_p1 * " + mX + ", " + std::to_string(m_mass_category == "low" ? 3 : 4) + "), " + std::to_string(m_mass_category == "low" ? 12 : 20) + ")', EGE_sigma_p0, EGE_sigma_p1)").c_str());
      m_wk->factory(("expr::EGE_kL_Xhh_m" + mX + "('TMath::Min(TMath::Max(EGE_kL_p0 + EGE_kL_p1 * " + mX + " + EGE_kL_p2 / " + mX + ", 0.2), 4.0)', EGE_kL_p0, EGE_kL_p1, EGE_kL_p2)").c_str());
      m_wk->factory(("expr::EGE_kH_Xhh_m" + mX + "('TMath::Min(TMath::Max(EGE_kH_p0 + EGE_kH_p1 * " + mX + " + EGE_kH_p2 * " + mX + " * " + mX + ", 0.2), 4.0)', EGE_kH_p0, EGE_kH_p1, EGE_kH_p2)").c_str());
      m_wk->factory(("ExpGausExpPDF::EGE_Xhh_m" + mX + "(mass, EGE_mu_Xhh_m" + mX +", EGE_sigma_Xhh_m" + mX + ", EGE_kL_Xhh_m" + mX + ", EGE_kH_Xhh_m" + mX + ")").c_str());
      dynamic_cast<RooSimultaneous*>(m_wk->pdf("simultaneous_EGE"))->addPdf(*m_wk->pdf(("EGE_Xhh_m" + mX).c_str()), mX.c_str());
    }
    // Add Crystal Ball + Gaussian
    if (std::find(m_models.begin(), m_models.end(), "CBGA") != m_models.end()) {
      m_wk->factory(("expr::CBGA_alpha_CB_Xhh_m" + mX + "('CBGA_alpha_CB_p0 + CBGA_alpha_CB_p1 * " + mX + "', CBGA_alpha_CB_p0, CBGA_alpha_CB_p1)").c_str());
      m_wk->factory(("expr::CBGA_f_CB_Xhh_m" + mX + "('CBGA_f_CB_p0 + CBGA_f_CB_p1 * " + mX + "', CBGA_f_CB_p0, CBGA_f_CB_p1)").c_str());
      m_wk->factory(("expr::CBGA_mu_CB_Xhh_m" + mX + "('CBGA_mu_CB_p0 + CBGA_mu_CB_p1 * " + mX + "', CBGA_mu_CB_p0, CBGA_mu_CB_p1)").c_str());
      m_wk->factory(("expr::CBGA_mu_GA_Xhh_m" + mX + "('CBGA_mu_GA_p0 + CBGA_mu_GA_p1 * " + mX + "', CBGA_mu_GA_p0, CBGA_mu_GA_p1)").c_str());
      m_wk->factory(("expr::CBGA_sigma_CB_Xhh_m" + mX + "('CBGA_sigma_CB_p0 + CBGA_sigma_CB_p1 * " + mX + " + CBGA_sigma_CB_p2 * " + mX + " * " + mX + "', CBGA_sigma_CB_p0, CBGA_sigma_CB_p1, CBGA_sigma_CB_p2)").c_str());
      m_wk->factory(("expr::CBGA_sigma_GA_Xhh_m" + mX + "('CBGA_sigma_GA_p0 + CBGA_sigma_GA_p1 * " + mX + " + CBGA_sigma_GA_p2 * " + mX + " * " + mX + "', CBGA_sigma_GA_p0, CBGA_sigma_GA_p1, CBGA_sigma_GA_p2)").c_str());
      m_wk->factory(("CBShape::CBGA_CB_Xhh_m" + mX + "(mass, CBGA_mu_CB_Xhh_m" + mX + ", CBGA_sigma_CB_Xhh_m" + mX + ", CBGA_alpha_CB_Xhh_m" + mX + ", CBGA_n_CB)").c_str());
      m_wk->factory(("Gaussian::CBGA_GA_Xhh_m" + mX + "(mass, CBGA_mu_GA_Xhh_m" + mX + ", CBGA_sigma_GA_Xhh_m" + mX + ")").c_str());
      m_wk->factory(("SUM::CBGA_Xhh_m" + mX + "(CBGA_f_CB_Xhh_m" + mX + " * CBGA_CB_Xhh_m" + mX + ", CBGA_GA_Xhh_m" + mX + ")").c_str());
      dynamic_cast<RooSimultaneous*>(m_wk->pdf("simultaneous_CBGA"))->addPdf(*m_wk->pdf(("CBGA_Xhh_m" + mX).c_str()), mX.c_str());
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
    // Write EGE to output
    for (auto metaparameter : m_model_metaparameters.at("EGE")) {
      m_wk->var(metaparameter.c_str())->setConstant();
    }
    m_wk->factory("expr::EGE_mu_raw('EGE_mu_p0 + EGE_mu_p1 * mass_resonance', EGE_mu_p0, EGE_mu_p1, mass_resonance)");
    m_wk->factory(("expr::EGE_sigma('TMath::Min(TMath::Max(EGE_sigma_p0 + EGE_sigma_p1 * mass_resonance, " + std::to_string(m_mass_category == "low" ? 3 : 4) + "), " + std::to_string(m_mass_category == "low" ? 12 : 20) + ")', EGE_sigma_p0, EGE_sigma_p1, mass_resonance)").c_str());
    m_wk->factory("expr::EGE_kL('TMath::Min(TMath::Max(EGE_kL_p0 + EGE_kL_p1 * mass_resonance + EGE_kL_p2 / mass_resonance, 0.2), 4.0)', EGE_kL_p0, EGE_kL_p1, EGE_kL_p2, mass_resonance)");
    m_wk->factory("expr::EGE_kH('TMath::Min(TMath::Max(EGE_kH_p0 + EGE_kH_p1 * mass_resonance + EGE_kH_p2 * mass_resonance * mass_resonance, 0.2), 4.0)', EGE_kH_p0, EGE_kH_p1, EGE_kH_p2, mass_resonance)");
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
        m_wk->var("mass_resonance")->setVal(PlotStyle::resonance_masses(m_mass_category).at(idx));
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
          textBox.DrawLatex(0.74, 0.86, ("#mu_{EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_mu_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.82, ("#sigma_{EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_sigma_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.78, ("k_{L,EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_kL_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.74, ("k_{H,EGE} = " + PlotStyle::to_string(m_wk->function(("EGE_kH_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
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
          textBox.DrawLatex(0.74, 0.86, ("#alpha_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_alpha_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.82, ("f_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_f_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.78, ("#mu_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_mu_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.74, ("#mu_{GA} = " + PlotStyle::to_string(m_wk->function(("CBGA_mu_GA_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.70, ("n_{CB} = " + PlotStyle::to_string(m_wk->var("CBGA_n_CB")->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.66, ("#sigma_{CB} = " + PlotStyle::to_string(m_wk->function(("CBGA_sigma_CB_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
          textBox.DrawLatex(0.74, 0.62, ("#sigma_{GA} = " + PlotStyle::to_string(m_wk->function(("CBGA_sigma_GA_Xhh_m" + mX).c_str())->getVal(), 2)).c_str());
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
        canvas_individual.Print(("plots/" + m_mass_category + "Mass_" + m_tag_category + "tag/signal_model_" + model + "_" + m_mass_category + "Mass_" + m_tag_category + "tag_mX_" + mX + ".pdf").c_str());
        MSG_INFO("Created plots/" << m_mass_category << "Mass_" << m_tag_category << "tag/signal_model_" << model << "_" << m_mass_category << "Mass_" << m_tag_category << "tag_mX_" << mX << ".pdf");
      }
    }
  }

  void SignalModel::append(std::vector<double>& target, std::vector<double> source) {
    target.insert(target.end(), source.begin(), source.end());
  }

  double SignalModel::get_initial_value(const std::string& name) {
    if (m_mass_category == "low" && m_tag_category == "0") {
      if (name == "EGE_mu_p0") { return -1.47449437499; }
      if (name == "EGE_mu_p1") { return 1.005375375; }
      if (name == "EGE_sigma_p0") { return -14.7105371; }
      if (name == "EGE_sigma_p1") { return 0.0693068599999; }
      if (name == "EGE_kL_p0") { return -23.5578042939; }
      if (name == "EGE_kL_p1") { return 0.0371280575616; }
      if (name == "EGE_kL_p2") { return 3952.48944757; }
      if (name == "EGE_kH_p0") { return -0.814981007462; }
      if (name == "EGE_kH_p1") { return 0.00491339376534; }
      if (name == "EGE_kH_p2") { return 1.54040320163e-07; }
    } else if (m_mass_category == "high" && m_tag_category == "0") {
      if (name == "EGE_mu_p0") { return 2.10913833994; }
      if (name == "EGE_mu_p1") { return 0.996911067194; }
      if (name == "EGE_sigma_p0") { return -13.9163441936; }
      if (name == "EGE_sigma_p1") { return 0.039597383871; }
      if (name == "EGE_kL_p0") { return 4.9407264299; }
      if (name == "EGE_kL_p1") { return -0.00241285149081; }
      if (name == "EGE_kL_p2") { return -1608.13185168; }
      if (name == "EGE_kH_p0") { return -1.37161926041; }
      if (name == "EGE_kH_p1") { return 0.00494901864028; }
      if (name == "EGE_kH_p2") { return -2.81315062286e-06; }
    } if (m_mass_category == "low" && m_tag_category == "1") {
      if (name == "EGE_mu_p0") { return -2.45754812504; }
      if (name == "EGE_mu_p1") { return 1.010055125; }
      if (name == "EGE_sigma_p0") { return -9.76830174374; }
      if (name == "EGE_sigma_p1") { return 0.05133163375; }
      if (name == "EGE_kL_p0") { return -30.4009357541; }
      if (name == "EGE_kL_p1") { return 0.0458943038413; }
      if (name == "EGE_kL_p2") { return 5236.77268308; }
      if (name == "EGE_kH_p0") { return 0.736566804598; }
      if (name == "EGE_kH_p1") { return -0.00250744455118; }
      if (name == "EGE_kH_p2") { return 7.4759865861e-06; }
    } if (m_mass_category == "high" && m_tag_category == "1") {
      if (name == "EGE_mu_p0") { return 5.34054940692; }
      if (name == "EGE_mu_p1") { return 0.991423952569; }
      if (name == "EGE_sigma_p0") { return 2.50161581028; }
      if (name == "EGE_sigma_p1") { return 0.0144543648221; }
      if (name == "EGE_kL_p0") { return 0.568798720108; }
      if (name == "EGE_kL_p1") { return 0.000286392707423; }
      if (name == "EGE_kL_p2") { return 86.4115958043; }
      if (name == "EGE_kH_p0") { return 1.35761961112; }
      if (name == "EGE_kH_p1") { return -0.00181687865308; }
      if (name == "EGE_kH_p2") { return 1.10214839952e-06; }
    } if (m_mass_category == "low" && m_tag_category == "2") {
      if (name == "EGE_mu_p0") { return -0.0675574997714; }
      if (name == "EGE_mu_p1") { return 1.0008295; }
      if (name == "EGE_sigma_p0") { return -3.6432324; }
      if (name == "EGE_sigma_p1") { return 0.0281058400001; }
      if (name == "EGE_kL_p0") { return -57.1941470189; }
      if (name == "EGE_kL_p1") { return 0.0833342844111; }
      if (name == "EGE_kL_p2") { return 10019.6733521; }
      if (name == "EGE_kH_p0") { return 5.91007635148; }
      if (name == "EGE_kH_p1") { return -0.0295306763403; }
      if (name == "EGE_kH_p2") { return 4.21812194433e-05; }
    } if (m_mass_category == "high" && m_tag_category == "2") {
      if (name == "EGE_mu_p0") { return 3.48818577074; }
      if (name == "EGE_mu_p1") { return 0.99304486166; }
      if (name == "EGE_sigma_p0") { return -0.325412569162; }
      if (name == "EGE_sigma_p1") { return 0.0170567944664; }
      if (name == "EGE_kL_p0") { return 1.07816004197; }
      if (name == "EGE_kL_p1") { return -6.26645623641e-05; }
      if (name == "EGE_kL_p2") { return -7.67419139643; }
      if (name == "EGE_kH_p0") { return 0.860939466908; }
      if (name == "EGE_kH_p1") { return -0.000403346158662; }
      if (name == "EGE_kH_p2") { return 2.5412126773e-07; }
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
