// STL
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
// ROOT
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
// RooFit
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgusBG.h"
#include "RooCBShape.h"
#include "RooCFunction3Binding.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooGamma.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooLandau.h"
#include "RooLognormal.h"
#include "RooMsgService.h"
#include "RooNLLVar.h"
#include "RooNovosibirsk.h"
#include "RooPlot.h"
#include "RooPlot.h"
#include "RooPoisson.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
// Local
#include "ParameterSet.h"
#include "Logger.h"
#include "FitMassPoint.h"

int main(int /*argc*/, char** /*argv*/) {  
  using namespace SpuriousSignal;
  
  int MASS_STEP = 1000; //20;
    
  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  gErrorIgnoreLevel = kBreak;
  gStyle->SetOptTitle(0);
  
  // Setup colours and mass categories
  std::vector<int> colours = {kViolet, kGreen + 3, kBlue, kRed, kMagenta, kCyan};
  std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)},
                                                                  {"high", std::make_pair<int, int>(335, 1140)}};
  // std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)}};

  TFile f_inputWS("/afs/cern.ch/user/a/andari/public/Leo/SigParam_low_cbga_cut5/Parameterized/SM/res_SM_CBGA_Parameterized_workspace.root");
  RooWorkspace wk(*(RooWorkspace*)f_inputWS.Get("signalWS"));
  // Gaussian
  RooRealVar a_muGANom_SM_c2(*(RooRealVar*)wk.obj("a_muGANom_SM_c2")); a_muGANom_SM_c2.setConstant(true);
  RooRealVar b_muGANom_SM_c2(*(RooRealVar*)wk.obj("b_muGANom_SM_c2")); b_muGANom_SM_c2.setConstant(true);
  RooRealVar c_muGANom_SM_c2(*(RooRealVar*)wk.obj("c_muGANom_SM_c2")); c_muGANom_SM_c2.setConstant(true);
  RooRealVar a_sigmaGANom_SM_c2(*(RooRealVar*)wk.obj("a_sigmaGANom_SM_c2")); a_sigmaGANom_SM_c2.setConstant(true);
  RooRealVar b_sigmaGANom_SM_c2(*(RooRealVar*)wk.obj("b_sigmaGANom_SM_c2")); b_sigmaGANom_SM_c2.setConstant(true);
  // Crystal Ball
  RooRealVar a_muCBNom_SM_c2(*(RooRealVar*)wk.obj("a_muCBNom_SM_c2")); a_muCBNom_SM_c2.setConstant(true);
  RooRealVar b_muCBNom_SM_c2(*(RooRealVar*)wk.obj("b_muCBNom_SM_c2")); b_muCBNom_SM_c2.setConstant(true);
  RooRealVar c_muCBNom_SM_c2(*(RooRealVar*)wk.obj("c_muCBNom_SM_c2")); c_muCBNom_SM_c2.setConstant(true);
  RooRealVar a_sigmaCBNom_SM_c2(*(RooRealVar*)wk.obj("a_sigmaCBNom_SM_c2")); a_sigmaCBNom_SM_c2.setConstant(true);
  RooRealVar b_sigmaCBNom_SM_c2(*(RooRealVar*)wk.obj("b_sigmaCBNom_SM_c2")); b_sigmaCBNom_SM_c2.setConstant(true);
  RooRealVar nCB_SM_c2(*(RooRealVar*)wk.obj("nCB_SM_c2")); nCB_SM_c2.setConstant(true);
  RooRealVar a_alphaCB_SM_c2(*(RooRealVar*)wk.obj("a_alphaCB_SM_c2")); a_alphaCB_SM_c2.setConstant(true);
  RooRealVar b_alphaCB_SM_c2(*(RooRealVar*)wk.obj("b_alphaCB_SM_c2")); b_alphaCB_SM_c2.setConstant(true);
  // Combination
  RooRealVar fracCB_SM_c2(*(RooRealVar*)wk.obj("fracCB_SM_c2")); fracCB_SM_c2.setConstant(true);

  // std::vector<double> 

  for (auto mass_kv : mass_categories) {
    // Define parameter ranges
    double peak_min((mass_kv.first == "low") ? 260 : 0);
    double peak_max((mass_kv.first == "low") ? 280 : 1000);
    
    // Define data parameters
    RooRealVar mass("mass", "m_{yyjj}", mass_kv.second.first, mass_kv.second.second);
    // RooFormulaVar mass_scaled("mass_scaled", ("(mass / " + std::to_string(mass_kv.second.first) + ")").c_str(), RooArgList(mass));
    RooRealVar weight("weight", "event weight", 1e-10, 1e10);

    // Get the data
    RooDataSet* ptr_raw_data = RooDataSet::read(("input/m_yyjj_SM_bkg_" + mass_kv.first + "Mass_0tag_tightIsolated.csv").c_str(), RooArgList(mass, weight));
    RooDataSet data("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));
    MSG_INFO("Loaded " << data.numEntries() << " events");

    // Construct vector of mass points
    int mass_point = mass_kv.second.first - MASS_STEP;
    std::vector<double> mass_points((mass_kv.second.second - mass_kv.second.first) / MASS_STEP + 1); 
    std::generate(mass_points.begin(), mass_points.end(), [&mass_point, &MASS_STEP]{ return mass_point += MASS_STEP; }); 

    // Signal: CB+Gaus
    RooRealVar mass_resonance("mass_resonance", "mass of resonance", 305);
    RooFormulaVar signal_gaus_mean("signal_gaus_mean", "a_muGANom_SM_c2 + b_muGANom_SM_c2 * (mass_resonance / 100. - 1) + c_muGANom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muGANom_SM_c2, b_muGANom_SM_c2, c_muGANom_SM_c2, mass_resonance));
    RooFormulaVar signal_gaus_sigma("signal_gaus_sigma", "a_sigmaGANom_SM_c2 + b_sigmaGANom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaGANom_SM_c2, b_sigmaGANom_SM_c2, mass_resonance));
    RooGaussian signal_gaus("signal_gaus", "signal gaus", mass, signal_gaus_mean, signal_gaus_sigma);
    RooFormulaVar signal_CB_mean("signal_CB_mean", "a_muCBNom_SM_c2 + b_muCBNom_SM_c2 * (mass_resonance / 100. - 1) + c_muCBNom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muCBNom_SM_c2, b_muCBNom_SM_c2, c_muCBNom_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_sigma("signal_CB_sigma", "a_sigmaCBNom_SM_c2 + b_sigmaCBNom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaCBNom_SM_c2, b_sigmaCBNom_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_alpha("signal_CB_alpha", "a_alphaCB_SM_c2 + b_alphaCB_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_alphaCB_SM_c2, b_alphaCB_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_n("signal_CB_n", "nCB_SM_c2", RooArgList(nCB_SM_c2));
    RooCBShape signal_CB("signal_CB", "signal_CB", mass, signal_CB_mean, signal_CB_sigma, signal_CB_alpha, signal_CB_n);
    RooFormulaVar signal_frac_CB("signal_frac_CB", "fracCB_SM_c2", RooArgList(fracCB_SM_c2));
    RooAddPdf signal_PDF("signal_PDFa", "CB + Gaussian", RooArgList(signal_CB, signal_gaus), signal_frac_CB);
    
    // Novosibirsk
    RooRealVar novosibirsk_peak("novosibirsk_peak", "peak of Novosibirsk", 270, peak_min, peak_max);
    RooRealVar novosibirsk_tail("novosibirsk_tail", "tail of Novosibirsk", -1, -20, 20);
    RooRealVar novosibirsk_width("novosibirsk_width", "width of Novosibirsk", 30, 0, 500);
    RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", mass, novosibirsk_peak, novosibirsk_width, novosibirsk_tail);
    
    // Modified Gamma
    // RooRealVar gamma_alpha("gamma_alpha", "alpha of Gamma", (mass_kv.first == "low" ? 1.5 : 0.08), 0, 100);
    RooRealVar gamma_alpha0("gamma_alpha0", "alpha0 of Gamma", (mass_kv.first == "low" ? 1.5 : 0.08), 0, 1000);
    RooRealVar gamma_alpha1("gamma_alpha1", "alpha1 of Gamma", 0.003, -0.01, 0.01);
    RooFormulaVar gamma_alpha("gamma_alpha", "gamma_alpha0 + gamma_alpha1 * mass", RooArgList(mass, gamma_alpha0, gamma_alpha1));
    // RooRealVar gamma_theta("gamma_theta", "theta of Gamma", (mass_kv.first == "low" ? 55 : 220), 0, 1000);
    RooRealVar gamma_theta0("gamma_theta0", "theta0 of Gamma", (mass_kv.first == "low" ? 0.07 : 190), 0, 1000);
    RooRealVar gamma_theta1("gamma_theta1", "theta1 of Gamma", 0.1, -1, 1);
    RooFormulaVar gamma_theta("gamma_theta", "gamma_theta0 + gamma_theta1 * mass", RooArgList(mass, gamma_theta0, gamma_theta1));
    RooRealVar gamma_mu("gamma_mu", "minimum of Gamma", mass_kv.second.first, 0, mass_kv.second.second);
    // RooGamma modified_gamma("modified_gamma", "modified_gamma", mass, gamma_alpha, gamma_theta, gamma_mu);
    RooGenericPdf modified_gamma("modified_gamma", "modified_gamma","TMath::GammaDist(mass, gamma_alpha, gamma_mu, gamma_theta)", RooArgList(mass, gamma_alpha, gamma_theta, gamma_mu));

    // Inverse Gaussian
    // RooRealVar invgaus_mu("invgaus_mu", "mean of inverse Gaussian", 100, 0, 1000);
    RooRealVar invgaus_mu0("invgaus_mu0", "mu0 of inverse Gaussian", 100, 0, 10000);
    RooRealVar invgaus_mu1("invgaus_mu1", "mu1 of inverse Gaussian", 0, -1, 1);
    RooFormulaVar invgaus_mu("invgaus_mu", "invgaus_mu0 + invgaus_mu1 * mass", RooArgList(mass, invgaus_mu0, invgaus_mu1));
    // RooRealVar invgaus_lambda("invgaus_lambda", "lambda of inverse Gaussian", 0.25, 0, 100);
    RooRealVar invgaus_lambda0("invgaus_lambda0", "lambda0 of inverse Gaussian", 0.25, 0, 100);
    RooRealVar invgaus_lambda1("invgaus_lambda1", "lambda1 of inverse Gaussian", 0, -1, 1);
    RooFormulaVar invgaus_lambda("invgaus_lambda", "invgaus_lambda0 + invgaus_lambda1 * mass", RooArgList(mass, invgaus_lambda0, invgaus_lambda1));
    RooGenericPdf invgaus("invgaus", "invgaus", "TMath::Sqrt(invgaus_lambda / (2 * TMath::Pi() * mass * mass * mass)) * TMath::Exp((-invgaus_lambda * (mass - invgaus_mu) * (mass - invgaus_mu)) / (2 * invgaus_mu * invgaus_mu * mass))", RooArgList(mass, invgaus_mu, invgaus_lambda));

    // Modified Landau
    RooRealVar landau_mean("landau_mean", "mean of Landau", 270, peak_min, peak_max);
    RooRealVar landau_sigma0("landau_sigma0", "sigma0 of Landau", 100, 0, 1000);
    RooRealVar landau_sigma1("landau_sigma1", "sigma0 of Landau", 0, -1, 1);
    RooFormulaVar landau_sigma("landau_sigma", "landau_sigma0 + landau_sigma1 * mass", RooArgList(mass, landau_sigma0, landau_sigma1));
    // RooLandau modified_landau("modified_landau","modified_landau", mass, landau_mean, landau_sigma);
    RooGenericPdf modified_landau("modified_landau", "modified_landau", "TMath::Landau(mass, landau_mean, landau_sigma)", RooArgList(mass, landau_mean, landau_sigma));

    // // Gamma + Gaussian
    // RooRealVar gammagaus_alpha("gammagaus_alpha", "gamma of Gamma", (mass_kv.first == "low" ? 1.4 : 0.0), 0, 100);
    // RooRealVar gammagaus_theta("gammagaus_theta", "theta of Gamma", (mass_kv.first == "low" ? 65 : 190), 0, 1000);
    // RooRealVar gammagaus_mu("gammagaus_mu", "minimum of Gamma", mass_kv.second.first, 0, mass_kv.second.second);
    // RooGamma gammagaus_Gamma("gammagaus_Gamma", "gammagaus_Gamma", mass, gammagaus_alpha, gammagaus_theta, gammagaus_mu);
    // RooRealVar gammagaus_peak("gammagaus_peak", "peak of Gaussian", 270, peak_min - 20, peak_max + 20);
    // RooRealVar gammagaus_sigma("gammagaus_sigma", "sigma of Gaussian", 18, 0, 1000);
    // RooGaussian gammagaus_gauss("gammagaus_gauss", "gammagaus_gauss", mass, gammagaus_peak, gammagaus_sigma);
    // RooRealVar gammagaus_fracgamma("gammagaus_fracgamma", "gammagaus_fracgamma", 0.8, 0.0, 1.0);
    // RooAddPdf gammagaus("gammagaus", "gammagaus", RooArgList(gammagaus_Gamma, gammagaus_gauss), gammagaus_fracgamma);
        
    // // Modified Cauchy
    // // RooRealVar cauchy_s1("cauchy_s1", "s1 of Cauchy", 0.01, -0.1, 0.1);
    // // RooFormulaVar cauchy_s("cauchy_s", "@1 * exp(@2 * @0)", RooArgList(mass, cauchy_s0, cauchy_s1)); // exponential as we require s > 0 everywhere
    // RooRealVar cauchy_mu("cauchy_mu", "location of Cauchy", 270, peak_min, peak_max);
    // RooRealVar cauchy_s0("cauchy_s0", "s0 of Cauchy", (mass_kv.first == "low" ? 200 : 150), 0, 10000);
    // RooRealVar cauchy_s1("cauchy_s1", "s1 of Cauchy",  (mass_kv.first == "low" ? 30.0 : 2.0), -1000, 1000);
    // RooFormulaVar cauchy_s("cauchy_s", ("cauchy_s0 * exp(cauchy_s1 * ((mass / " + std::to_string(mass_kv.second.first) + ") - 1.0))").c_str(), RooArgList(mass, cauchy_s0, cauchy_s1)); // exponential as we require beta > 0 everywhere
    // RooGenericPdf modified_cauchy("modified_cauchy", "modified_cauchy", "TMath::CauchyDist(mass, cauchy_mu, cauchy_s)", RooArgList(mass, cauchy_mu, cauchy_s));

    // // Log Normal
    // RooRealVar log_normal_m0("log_normal_m0", "median of log normal", 270, peak_min, peak_max);
    // RooRealVar log_normal_k("log_normal_k", "exp(shape) of log normal", 5, -100, 100);
    // // RooRealVar log_normal_k0("log_normal_k0", "k0 of log normal", 50, -1000, 1000);
    // // RooRealVar log_normal_k1("log_normal_k1", "k1 of log normal", 0, -0.1, 0.1);
    // // RooFormulaVar log_normal_k("log_normal_k", ("log_normal_k0 + log_normal_k1 * ((mass / " + std::to_string(mass_kv.second.first) + ") - 1.0)").c_str(), RooArgList(mass, log_normal_k0, log_normal_k1));
    // RooLognormal log_normal("log_normal", "log_normal", mass, log_normal_m0, log_normal_k);
            
    // // Rayleigh
    // RooRealVar rayleigh_shift("rayleigh_shift", "location of Rayleigh", 270, peak_min - 50, peak_max + 50);
    // RooRealVar rayleigh_sigma("rayleigh_sigma", "sigma of Cauchy", 0, 10000);
    // RooGenericPdf rayleigh("rayleigh", "rayleigh", "((mass - rayleigh_shift) / (rayleigh_sigma * rayleigh_sigma)) * TMath::Exp((-(mass - rayleigh_shift) * (mass - rayleigh_shift)) / (2 * rayleigh_sigma * rayleigh_sigma))", RooArgList(mass, rayleigh_shift, rayleigh_sigma));


    // // Johnson
    // RooRealVar johnson_xi("johnson_xi", "location of Johnson", 270, peak_min - 50, peak_max + 50);
    // RooRealVar johnson_lambda("johnson_lambda", "location of Johnson", -100, 100);
    // RooRealVar johnson_gamma("johnson_gamma", "location of Johnson", -100, 100);
    // RooRealVar johnson_delta("johnson_delta", "location of Johnson", -100, 100);
    // RooJohnsonSU johnson("johnson","johnson", mass, johnson_xi, johnson_lambda, johnson_gamma, johnson_delta);

    // // Inverted Argus
    // RooRealVar argus_m0("argus_m0", "m0 of Argus", 270, 0, 1000);
    // RooRealVar argus_c("argus_c", "range of Argus", mass_kv.second.second, 200, 2000);
    // RooFormulaVar argus_mass_inv("argus_mass_inv", "argus_c - mass", RooArgList(argus_c, mass));
    // RooArgusBG inverted_argus("inverted_argus","inverted_argus", argus_mass_inv, argus_m0, argus_c);
    // 
    // // Poisson
    // RooRealVar poisson_shift("poisson_shift", "shift of Poisson", 270, peak_min - 50, peak_max + 50);
    // RooRealVar poisson_lambda("poisson_lambda", "lambda of Poisson", 0, 1000);
    // RooFormulaVar poisson_mass("poisson_mass", "mass - poisson_shift", RooArgList(mass, poisson_shift));
    // RooPoisson poisson("inverted_argus","inverted_argus", poisson_mass, poisson_lambda);

    // // Convolved Gaussian + exponential
    // RooRealVar gaussexp_mg("gaussexp_mg", "gaussexp_mg", 270, peak_min - 50, peak_max + 50);
    // RooRealVar gaussexp_sg("gaussexp_sg", "gaussexp_sg", 20, 0, 1000);
    // RooGaussian gaussexp_gauss("gaussexp_gauss", "gaussexp_gauss", mass, gaussexp_mg, gaussexp_sg);
    // RooRealVar gaussexp_c("gaussexp_c","gaussexp_c", -0.006, -20, 20);
    // RooExponential gaussexp_exp("gaussexp_exp", "gaussexp_exp", mass, gaussexp_c);
    // // // Set #bins to be used for FFT sampling to 10000
    // // t.setBins(10000,"cache") ; 
    // // Construct landau (x) gauss
    // RooFFTConvPdf gaussexp("gaussexp", "gaussexp_gauss (X) gaussexp_exp", mass, gaussexp_gauss, gaussexp_exp);
      
    
    // Number of signal and background events
    RooRealVar nSig("nSig", "number of signal events", 0, 100);
    RooRealVar nBkg("nBkg", "number of background events", 1000, 10000);

    // Setup fit functions
    std::vector<RooAbsPdf*> bkg_functions;
    bkg_functions.push_back(&novosibirsk);
    bkg_functions.push_back(&modified_gamma);
    // bkg_functions.push_back(&gammagaus);
    // bkg_functions.push_back(&invgaus);
    bkg_functions.push_back(&modified_landau);
    // bkg_functions.push_back(&modified_cauchy);
    // bkg_functions.push_back(&log_normal);
    // bkg_functions.push_back(&rayleigh);
    // bkg_functions.push_back(&inverted_argus);
    // bkg_functions.push_back(&poisson);
    // bkg_functions.push_back(&gaussexp);
    

    // Construct parameter sets that need to be remembered
    // std::vector<ParameterSet> parameter_sets;
    // parameter_sets.push_back(ParameterSet("Novosibirsk", {&novosibirsk_peak, &novosibirsk_width, &novosibirsk_tail, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Modified Gamma", {&gamma_gamma0, &gamma_gamma1, &gamma_mu, &gamma_beta0, &gamma_beta1, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Gammagaus", {&gammagaus_gamma, &gammagaus_mu, &gammagaus_beta, &gammagaus_peak, &gammagaus_sigma, &gammagaus_fracgamma, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Modified Cauchy", {&cauchy_mu, &cauchy_s0, &cauchy_s1, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Log Normal", {&log_normal_m0, &log_normal_k0, &log_normal_k1, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Rayleigh", {&rayleigh_shift, &rayleigh_sigma, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Inverted Argus", {&argus_m0, &argus_c, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Poisson", {&poisson_shift, &poisson_lambda, &nSig, &nBkg, &mass_resonance}));
    // parameter_sets.push_back(ParameterSet("Gaussian * exponential", {&gaussexp_mg, &gaussexp_sg, &gaussexp_c, &nSig, &nBkg, &mass_resonance}));

    // Recreate output graph file
    TFile f_root_output((std::string("output/fit_functions_") + mass_kv.first + ".root").c_str(), "RECREATE");

    // Do background-only fits
    MSG_INFO("Performing background-only fits for " << bkg_functions.size() << " fit functions.");
    FitMassPoint fits_bkg_only(data, bkg_functions, mass_kv.first);
    fits_bkg_only.fit(true);
    fits_bkg_only.plot(mass.frame(), -1, f_root_output);
    
    double m(invgaus_mu.getVal()), l(invgaus_lambda.getVal());
    std::cout << "Inv gaus, mu: " << m << ", lambda: " << l << ", mode: " << m * (sqrt(1 + ((9 * m * m) / (4 * l *l))) - (3 * m) / (2 * l)) << std::endl;
    
    // Close output file
    f_root_output.Close();

    // // Do S+B fits for different backgrounds
    // std::vector<RooAbsPdf*> splusb_functions;
    // for (auto bkg_function : bkg_functions) {
    //   splusb_functions.push_back(new RooAddPdf("signal_plus_" + bkg_function->getTitle(), "signal + " + bkg_function->getTitle(), RooArgList(signal_PDF, *bkg_function), RooArgList(nSig, nBkg)));
    // }
    // FitMassPoint fits_splusb(data, splusb_functions, mass_kv.first);
    // MSG_INFO("Will perform fits for \033[1m" << mass_points.size() << "\033[0m mass points in the range \033[1m" << mass_kv.second.first << " - " << mass_kv.second.second << "\033[0m GeV");
    // for (auto mass_point : mass_points) {
    //   mass_resonance.setVal(mass_point);
    //   mass_resonance.setConstant(true);
    //   MSG_INFO("... working on mX = " << mass_point);
    //   fits_splusb.fit();
    //   fits_splusb.plot(mass.frame(), mass_point);
    // }

    // // Recreate output text file
    // std::ofstream f_output;
    // f_output.open("output/spurious_signal_" + mass_kv.first + ".csv", std::ios::trunc);
    // f_output.close();
    //     
    // // Selection criterion
    // for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
    //   // Selection criterion
    //   double Z_spurious(0), Z_spurious_max(-99);
    //   MSG_INFO("Using \033[1m" << mass_kv.first << " mass \033[0m fit function \033[1m" << splusb_functions.at(idx)->getTitle() << "\033[0m");
    // 
    //   // Iterate over mass points
    //   for (int iMass = mass_kv.second.first; iMass < mass_kv.second.second; iMass += MASS_STEP) {
    //     
    //     fit_status = splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"), RooFit::Hesse(false), RooFit::Minos(false), RooFit::PrintLevel(PRINT_LEVEL), RooFit::Save(true))->status();
    //     if (fit_status != 0) { MSG_ERROR("... simplex fit did not converge: status = " << fit_status); }
    //     MSG_DEBUG("... simplex: " << parameter_sets.at(idx));
    //     fit_status = splusb_functions.at(idx)->fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(false), RooFit::PrintLevel(PRINT_LEVEL), RooFit::Save(true))->status();
    //     if (fit_status != 0) { MSG_ERROR("... migrad + improve + hesse fit did not converge: status = " << fit_status); }
    //     MSG_INFO("... migrad + improve + hesse: " << parameter_sets.at(idx));
    //     // Record these parameters if Z > max(Z) so far
    //     Z_spurious = nSig.getVal() / nSig.getError(); // (N_S +/- errS_SumW2_true) / (errS_SumW2_false)
    //     MSG_INFO("... resonance mass = " << mass_resonance.getVal() << " => nSig: " << nSig.getVal() << " +/- " << nSig.getError() << ", nBkg: " << nBkg.getVal() << " +/- " << nBkg.getError() << " [" << fit_status << "] => Z = " << Z_spurious);
    //     if (Z_spurious > Z_spurious_max) {
    //       Z_spurious_max = Z_spurious;
    //       MSG_INFO("Recording new max(z) = " << Z_spurious_max << " at " << mass_resonance.getVal());
    //       parameter_sets.at(idx).record_values();
    //     }
    //     // Write parameters to output file
    //     f_output.open("output/spurious_signal_" + mass_kv.first + ".csv", std::ios::app);
    //     f_output << bkg_functions.at(idx)->getTitle() << " " << mass_resonance.getVal() << " " << nSig.getVal() << " " << nSig.getError() << " " << nBkg.getVal() << " " << nBkg.getError() << " " << fit_status << " " << Z_spurious << std::endl;
    //     f_output.close();
    //   }
    // }
    // 
    // // Recreate output graph file and open a new canvas for plotting
    // TFile f_root_output((std::string("output/fit_functions_") + mass_kv.first + ".root").c_str(), "RECREATE");
    // TCanvas canvas("canvas", "canvas", 800, 600);
    // 
    // // Plot PDF fits to the data (MC in this case)
    // RooPlot* frame = mass.frame();
    // data.plotOn(frame);
    // f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "data");
    // for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
    //   parameter_sets.at(idx).restore_values();
    //   MSG_INFO("Restored " << bkg_functions.at(idx)->getTitle() << " parameters at resonant mass = " << mass_resonance.getVal() << ": nSig = " << nSig.getVal() << ", nBkg = " << nBkg.getVal());
    //   // Background
    //   splusb_functions.at(idx)->plotOn(frame, RooFit::Components(bkg_functions.at(idx)->getTitle()), RooFit::LineColor(colours.at(idx)), RooFit::LineStyle(kDashed));
    //   f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "bkg_" + bkg_functions.at(idx)->getTitle());
    //   // Signal + background
    //   splusb_functions.at(idx)->plotOn(frame, RooFit::LineColor(colours.at(idx)));
    //   f_root_output.WriteObject((TGraph*)frame->getObject(frame->numItems() - 1), "splusb_" + bkg_functions.at(idx)->getTitle());
    // }
    // // for (int i = 0; i < frame->numItems(); ++i) { MSG_INFO(frame->getObject(i)->GetName()); }
    // frame->Draw();
    // frame->Draw();
    // TLegend legend(0.6, 0.7, 0.9, 0.9);
    // legend.AddEntry(frame->findObject("signal_plus_novosibirsk_Norm[mass]"), "Novosibirsk", "L");
    // legend.AddEntry(frame->findObject("signal_plus_modified_gamma_Norm[mass]"), "Modified Gamma", "L");
    // legend.AddEntry(frame->findObject("signal_plus_modified_cauchy_Norm[mass]"), "Modified Cauchy", "L");
    // legend.Draw();    
    // canvas.Print((std::string("output/m_yyjj_") + mass_kv.first + "Mass.pdf").c_str());
    // MSG_INFO("Saved " + mass_kv.first + " mass plot");
    // f_root_output.Close();
  }
  return 0;
}

