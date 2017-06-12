// // Local
// #include "D3PDReader.h"
// #include "AnalysisFactory.h"
// #include "EventLoop.h"
// #include "FileHandler.h"
// #include "Interpreter.h"
// #include "TrigDecisionToolD3PDWrapper.h"
// STL
#include <string>
#include <vector>
#include <iostream>
#include <utility>
// #include <fstream>
// // BOOST
// #include <boost/foreach.hpp>
// #include <boost/lexical_cast.hpp>
// ROOT
#include "TCanvas.h"
#include "TError.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
// RooFit
// #include "RooFit.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooGenericPdf.h"
#include "RooNovosibirsk.h"
#include "RooGamma.h"
#include "RooPlot.h"
#include "RooMsgService.h"
#include "RooCFunction3Binding.h"
#include "RooNLLVar.h"
#include "RooWorkspace.h"
#include "RooCBShape.h"
// #include "RooMinuit.h"
// Google
#ifdef DEBUG
#include <google/profiler.h>
#endif

int main(int /*argc*/, char** /*argv*/) {
  #ifdef DEBUG
  ProfilerStart("out.prof");
  #endif
  
  // Disable RooFit and ROOT messages
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR); //WARNING FATAL
  gErrorIgnoreLevel = kBreak;

  // Setup colours and mass categories
  std::vector<int> colours = {kViolet, kGreen + 3, kBlue, kRed};
  // std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)},
  //                                                                 {"high", std::make_pair<int, int>(335, 1140)}};
  std::map< std::string, std::pair<int, int> > mass_categories = {{"low", std::make_pair<int, int>(245, 485)}};

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

// (a_alphaCB_SM_c0,a_alphaCB_SM_c1,a_alphaCB_SM_c2,a_muCBNom_SM_c0,a_muCBNom_SM_c1,a_muCBNom_SM_c2,a_muGANom_SM_c0,a_muGANom_SM_c1,a_muGANom_SM_c2,a_sigmaCBNom_SM_c0,a_sigmaCBNom_SM_c1,a_sigmaCBNom_SM_c2,a_sigmaGANom_SM_c0,a_sigmaGANom_SM_c1,a_sigmaGANom_SM_c2,b_alphaCB_SM_c0,b_alphaCB_SM_c1,b_alphaCB_SM_c2,b_muCBNom_SM_c0,b_muCBNom_SM_c1,b_muCBNom_SM_c2,b_muGANom_SM_c0,b_muGANom_SM_c1,b_muGANom_SM_c2,b_sigmaCBNom_SM_c0,b_sigmaCBNom_SM_c1,b_sigmaCBNom_SM_c2,b_sigmaGANom_SM_c0,b_sigmaGANom_SM_c1,b_sigmaGANom_SM_c2,c_muCBNom_SM_c0,c_muCBNom_SM_c1,c_muCBNom_SM_c2,c_muGANom_SM_c0,c_muGANom_SM_c1,c_muGANom_SM_c2,fracCB_SM_c0,fracCB_SM_c1,fracCB_SM_c2,functionName,mResonance,m_yy,m_yy_m1000000_c0,m_yy_m1000000_c1,m_yy_m1000000_c2,m_yy_m260000_c0,m_yy_m260000_c1,m_yy_m260000_c2,m_yy_m275000_c0,m_yy_m275000_c1,m_yy_m275000_c2,m_yy_m300000_c0,m_yy_m300000_c1,m_yy_m300000_c2,m_yy_m325000_c0,m_yy_m325000_c1,m_yy_m325000_c2,m_yy_m350000_c0,m_yy_m350000_c1,m_yy_m350000_c2,m_yy_m400000_c0,m_yy_m400000_c1,m_yy_m400000_c2,m_yy_m450000_c0,m_yy_m450000_c1,m_yy_m450000_c2,m_yy_m500000_c0,m_yy_m500000_c1,m_yy_m500000_c2,m_yy_m750000_c0,m_yy_m750000_c1,m_yy_m750000_c2,nCB_SM_c0,nCB_SM_c1,nCB_SM_c2,sigYield_SM_m1000000_c0,sigYield_SM_m1000000_c1,sigYield_SM_m1000000_c2,sigYield_SM_m260000_c0,sigYield_SM_m260000_c1,sigYield_SM_m260000_c2,sigYield_SM_m275000_c0,sigYield_SM_m275000_c1,sigYield_SM_m275000_c2,sigYield_SM_m300000_c0,sigYield_SM_m300000_c1,sigYield_SM_m300000_c2,sigYield_SM_m325000_c0,sigYield_SM_m325000_c1,sigYield_SM_m325000_c2,sigYield_SM_m350000_c0,sigYield_SM_m350000_c1,sigYield_SM_m350000_c2,sigYield_SM_m400000_c0,sigYield_SM_m400000_c1,sigYield_SM_m400000_c2,sigYield_SM_m450000_c0,sigYield_SM_m450000_c1,sigYield_SM_m450000_c2,sigYield_SM_m500000_c0,sigYield_SM_m500000_c1,sigYield_SM_m500000_c2,sigYield_SM_m750000_c0,sigYield_SM_m750000_c1,sigYield_SM_m750000_c2,signalCates0,signalCates1,signalCates2,wt,yieldVar_a_SM_c0,yieldVar_a_SM_c1,yieldVar_a_SM_c2,yieldVar_b_SM_c0,yieldVar_b_SM_c1,yieldVar_b_SM_c2,yieldVar_c_SM_c0,yieldVar_c_SM_c1,yieldVar_c_SM_c2,yieldVar_d_SM_c0,yieldVar_d_SM_c1,yieldVar_d_SM_c2)
  // wk.Print();
  // RooFormulaVar* alphaCB_SM_c2 = (RooFormulaVar*)wk.function("alphaCB_SM_c2");
  // std::cout << ((RooRealVar*)wk.obj("a_alphaCB_SM_c2"))->getVal() << std::endl;
  

  // alphaCB_SM_c2->Print();
  // alphaCB_SM_c2->Ge
  // a_alphaCB_SM_c2,b_alphaCB_SM_c2,mRegularized,mResonance
  // std::cout << "alphaCB_SM_c2->getParameter(0): " << (RooRealVar*)(alphaCB_SM_c2->getParameter(0))->getTitle() << ", " << (RooRealVar*)(alphaCB_SM_c2->getParameter(0))->getVal() << std::endl;
  // RooFormulaVar muCBNom_SM_m260000_c1

    // RooCBShape RooAbsReal &_m, RooAbsReal &_m0, RooAbsReal &_sigma, RooAbsReal &_alpha, RooAbsReal &_n)
  // RooCBShape::pdfCB_SM_c0

  // RooCBShape::pdfCB_SM_m750000_c1[ m=m_yy_m750000_c1 m0=muCB_SM_m750000_c1 sigma=sigmaCB_SM_m750000_c1 alpha=alphaCB_SM_m750000_c1 n=nCB_SM_c1 ] = 0.5784
  // RooProduct::muCB_SM_m260000_c1[ muCBNom_SM_m260000_c1 ] = 262.832
  // RooFormulaVar::muCBNom_SM_m260000_c1[ actualVars=(a_muCBNom_SM_c1,b_muCBNom_SM_c1,c_muCBNom_SM_c1) formula="@0+@1*1.600000+@2*1.600000*1.600000+260.000000" ] = 262.832
// RooAddPdf::sigPdf_SM_c1[ fracCB_SM_c1 * pdfCB_SM_c1 + [%] * pdfGA_SM_c1 ] = 0.169803


// RooFormulaVar::alphaCB_SM_c2[ actualVars=(a_alphaCB_SM_c2,b_alphaCB_SM_c2,mRegularized,mResonance) formula="@0+@1*@2" ] = 1.72284
// RooFormulaVar::muCBNom_SM_c2[ actualVars=(a_muCBNom_SM_c2,b_muCBNom_SM_c2,c_muCBNom_SM_c2,mRegularized,mResonance) formula="@0+@1*@3+@2*@3*@3+@4" ] = 1012.89
// // RooProduct::muCB_SM_c2[ muCBNom_SM_c2 ] = 1012.89
// RooFormulaVar::muGANom_SM_c2[ actualVars=(a_muGANom_SM_c2,b_muGANom_SM_c2,c_muGANom_SM_c2,mRegularized,mResonance) formula="@0+@1*@3+@2*@3*@3+@4" ] = 996.2
// // RooProduct::muGA_SM_c2[ muGANom_SM_c2 ] = 996.2
// RooFormulaVar::sigmaCBNom_SM_c2[ actualVars=(a_sigmaCBNom_SM_c2,b_sigmaCBNom_SM_c2,mRegularized,mResonance) formula="@0+@1*@2" ] = 44.9995
// // RooProduct::sigmaCB_SM_c2[ sigmaCBNom_SM_c2 ] = 44.9995
// RooFormulaVar::sigmaGANom_SM_c2[ actualVars=(a_sigmaGANom_SM_c2,b_sigmaGANom_SM_c2,mRegularized,mResonance) formula="@0+@1*@2" ] = 17.1496
// // RooProduct::sigmaGA_SM_c2[ sigmaGANom_SM_c2 ] = 17.1496


  for (auto kv : mass_categories) {
    // Define data parameters
    RooRealVar mass("mass", "m_{yyjj}", kv.second.first, kv.second.second);
    RooRealVar weight("weight", "event weight", 1e-10, 1e10);

    // Get the data
    std::string input_directory("/afs/cern.ch/user/j/jrobinso/HGamma/PyAnalysis/text/");
    RooDataSet* ptr_raw_data = RooDataSet::read((input_directory + "m_yyjj_SM_bkg_" + kv.first + "Mass_0tag_tightIsolated.txt").c_str(), RooArgList(mass, weight));
    RooDataSet data("data", "data", RooArgSet(mass, weight), RooFit::Import(*ptr_raw_data), RooFit::WeightVar(weight));

    // Signal: CB+Gaus
    RooRealVar mass_resonance("mass_resonance", "mass of resonance", 305);
    mass_resonance.setConstant(true);
    // RooRealVar signal_gaus_mean("signal_gaus_mean", "mean of signal gaussian", 30, 0, 100);
    // RooRealVar signal_gaus_sigma("signal_gaus_sigma", "width of signal gaussian", -1, -10, 10);
    RooFormulaVar signal_gaus_mean("signal_gaus_mean", "a_muGANom_SM_c2 + b_muGANom_SM_c2 * (mass_resonance / 100. - 1) + c_muGANom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muGANom_SM_c2, b_muGANom_SM_c2, c_muGANom_SM_c2, mass_resonance));
    RooFormulaVar signal_gaus_sigma("signal_gaus_sigma", "a_sigmaGANom_SM_c2 + b_sigmaGANom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaGANom_SM_c2, b_sigmaGANom_SM_c2, mass_resonance));
    RooGaussian signal_gaus("signal_gaus", "signal gaus", mass, signal_gaus_mean, signal_gaus_sigma);
    // RooRealVar signal_CB_mean("signal_CB_mean", "mean of Crystal Ball", 260, 200, 500);
    // RooRealVar signal_CB_sigma("signal_CB_sigma", "width of Crystal Ball", 30, 0, 100);
    // RooRealVar signal_CB_alpha("signal_CB_alpha", "alpha of Crystal Ball", -1, -10, 10);
    // RooRealVar signal_CB_n("signal_CB_n", "n of Crystal Ball", -1, -10, 10);
    RooFormulaVar signal_CB_mean("signal_CB_mean", "a_muCBNom_SM_c2 + b_muCBNom_SM_c2 * (mass_resonance / 100. - 1) + c_muCBNom_SM_c2 * (mass_resonance / 100. - 1) * (mass_resonance / 100. - 1) + mass_resonance", RooArgList(a_muCBNom_SM_c2, b_muCBNom_SM_c2, c_muCBNom_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_sigma("signal_CB_sigma", "a_sigmaCBNom_SM_c2 + b_sigmaCBNom_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_sigmaCBNom_SM_c2, b_sigmaCBNom_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_alpha("signal_CB_alpha", "a_alphaCB_SM_c2 + b_alphaCB_SM_c2 * (mass_resonance / 100. - 1)", RooArgList(a_alphaCB_SM_c2, b_alphaCB_SM_c2, mass_resonance));
    RooFormulaVar signal_CB_n("signal_CB_n", "nCB_SM_c2", RooArgList(nCB_SM_c2));
    RooCBShape signal_CB("signal_CB", "signal_CB", mass, signal_CB_mean, signal_CB_sigma, signal_CB_alpha, signal_CB_n);
    // RooRealVar signal_frac_CB("signal_frac_CB", "fraction of Crystal Ball", 0.5, 0, 1.0);
    RooFormulaVar signal_frac_CB("signal_frac_CB", "fracCB_SM_c2", RooArgList(fracCB_SM_c2));
    // RooFormulaVar signal_frac_gaus("signal_frac_gaus", "1.0 - signal_frac_CB", RooArgList(signal_frac_CB)); 
    // std::cout << "signal_frac_CB: " << signal_frac_CB.getVal() << std::endl;   
    // RooAddPdf signal_PDF("signal_PDFa", "CB + Gaussian", RooArgList(signal_CB, signal_gaus), RooArgList(signal_frac_CB, signal_frac_gaus));
    RooAddPdf signal_PDF("signal_PDFa", "CB + Gaussian", RooArgList(signal_CB, signal_gaus), signal_frac_CB);

    // for( int i = 340; i < 360; ++i) {
    //   mass.setVal(i); //change value
    //   std::cout << "signal_gaus_mean: " << mass.getVal() << " -> " << signal_gaus_mean.getVal() << std::endl;
    //   std::cout << "signal_gaus_sigma: " << mass.getVal() << " -> " << signal_gaus_sigma.getVal() << std::endl;
    //   std::cout << "signal_gaus: " << mass.getVal() << " -> " << signal_gaus.getVal() << std::endl;
    //   std::cout << "signal_CB_mean: " << mass.getVal() << " -> " << signal_CB_mean.getVal() << std::endl;
    //   std::cout << "signal_CB_sigma: " << mass.getVal() << " -> " << signal_CB_sigma.getVal() << std::endl;
    //   std::cout << "signal_CB_alpha: " << mass.getVal() << " -> " << signal_CB_alpha.getVal() << std::endl;
    //   std::cout << "signal_CB_n: " << mass.getVal() << " -> " << signal_CB_n.getVal() << std::endl;
    //   std::cout << "signal_CB: " << mass.getVal() << " -> " << signal_CB.getVal() << std::endl;
    //   std::cout << "signal_frac_CB: " << mass.getVal() << " -> " << signal_frac_CB.getVal() << std::endl;
    // }
    // signal_PDF.Print();
    

    // Novosibirsk
    RooRealVar novosibirsk_peak("novosibirsk_peak", "peak of Novosibirsk", 260, 200, 500);
    RooRealVar novosibirsk_width("novosibirsk_width", "width of Novosibirsk", 30, 0, 100);
    RooRealVar novosibirsk_tail("novosibirsk_tail", "tail of Novosibirsk", -1, -10, 10);
    RooNovosibirsk novosibirsk("novosibirsk", "novosibirsk", mass, novosibirsk_peak, novosibirsk_width, novosibirsk_tail);
    
    
    // Modified Cauchy
    RooRealVar cauchy_mu("cauchy_mu", "mean of Cauchy", -500, -1e4, 1e4);
    RooRealVar cauchy_t0("cauchy_t0", "t0 of Cauchy", 1000, -1e5, 1e5);
    RooRealVar cauchy_t1("cauchy_t1", "t1 of Cauchy", 0, -10, 10);
    RooFormulaVar cauchy_t("cauchy_t", "cauchy_t0 + cauchy_t1 * mass", RooArgList(cauchy_t0, cauchy_t1, mass));
    RooRealVar cauchy_s0("cauchy_s0", "s0 of Cauchy", -150, -1e5, 1e5);
    RooRealVar cauchy_s1("cauchy_s1", "s1 of Cauchy", 1, -100, 100);
    RooFormulaVar cauchy_s("cauchy_s", "cauchy_s0 + cauchy_s1 * mass", RooArgList(cauchy_s0, cauchy_s1, mass));
    RooGenericPdf modified_cauchy("modified_cauchy", "modified_cauchy", "TMath::CauchyDist(mass - cauchy_mu, cauchy_t, cauchy_s)", RooArgList(mass, cauchy_mu, cauchy_t, cauchy_s));
    
    // Modified Gamma
    RooRealVar gamma_Gamma("gamma_Gamma", "gamma of Gamma", 5, 0, 20);
    RooRealVar gamma_beta0("gamma_beta0", "beta0 of Gamma", -20, -20, 20);
    RooRealVar gamma_beta1("gamma_beta1", "beta1 of Gamma", 0.1, -1, 1);
    RooFormulaVar gamma_beta("gamma_beta", "TMath::Abs(gamma_beta0 + gamma_beta1 * mass)", RooArgList(gamma_beta0, gamma_beta1, mass));
    RooRealVar gamma_mu("gamma_mu", "mean of Gamma", 245, 220, 260);
    RooGenericPdf modified_gamma("modified_gamma", "modified_gamma", "TMath::GammaDist(mass, gamma_Gamma, gamma_mu, TMath::Abs(gamma_beta0 + gamma_beta1 * mass))", RooArgList(mass, gamma_Gamma, gamma_mu, gamma_beta0, gamma_beta1));
    
    // Setup fit functions
    std::vector<RooAbsPdf*> bkg_functions;
    bkg_functions.push_back(&novosibirsk);
    // bkg_functions.push_back(&modified_gamma);
    // bkg_functions.push_back(&modified_cauchy);
    RooRealVar nSig("nSig", "number of signal events", 10, 0, 1e3);
    RooRealVar nBkg("nBkg", "number of background events", 10, 0, 1e6);
    // RooAddPdf SplusB("SplusB", "Signal + background", RooArgList(signal_PDF, *bkg_function), RooArgList(nSig, nBkg));
    RooAddPdf SplusB("SplusB", "Signal + background", RooArgList(signal_PDF, novosibirsk), RooArgList(nSig, nBkg));
    
    // Fit PDFs to data
    std::vector<RooFitResult*> fit_results;
    for (auto bkg_function : bkg_functions) {
      SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "simplex"));
      SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit", "migradimproved"), RooFit::Hesse(true), RooFit::Minos(true));
      // SplusB.fitTo(data, RooFit::SumW2Error(true), RooFit::NumCPU(4), RooFit::Minimizer("Minuit2", "migrad"));
      // m.simplex();
      // m.migrad();
      // m.improve();
      // m.hesse();
      // m.minos();      
      std::cout << "nSig: " << nSig.getVal() << std::endl;
      std::cout << "nBkg: " << nBkg.getVal() << std::endl;
    }
    
    // Make a TCanvas
    TCanvas canvas("canvas", "canvas", 800, 600);

    // RooDataSet* data1 = signal_PDF.generate(mass, 10000);

    // Plot PDF
    RooPlot* frame = mass.frame();
    // data1->plotOn(frame);
    // signal_CB.plotOn(frame, RooFit::LineColor(kRed));
    // signal_gaus.plotOn(frame, RooFit::LineColor(kBlue));
    // signal_PDF.plotOn(frame, RooFit::LineColor(kBlack));
    data.plotOn(frame);
    for (unsigned idx = 0; idx < bkg_functions.size(); ++idx) {
      bkg_functions.at(idx)->plotOn(frame, RooFit::LineColor(colours.at(idx)), RooFit::LineStyle(kDashed));
    }
    SplusB.plotOn(frame, RooFit::LineColor(kRed));
    frame->Draw();
    // xframe->GetYaxis()->SetLimits()
    
    // for (int i = 0; i < frame->numItems(); ++i) {
    //   std::cout << frame->getObject(i)->GetName() << std::endl;
    // }
          
    // TLegend legend(0.6, 0.7, 0.9, 0.9);
    // legend.AddEntry(frame->findObject("novosibirsk_Norm[mass]"), "Novosibirsk", "L");
    // legend.AddEntry(frame->findObject("modified_cauchy_Norm[mass]"), "Cauchy", "L");
    // legend.AddEntry(frame->findObject("modified_gamma_Norm[mass]"), "Gamma", "L");
    // legend.Draw();
    
    canvas.Print((std::string("output/m_yyjj_") + kv.first + "Mass.pdf").c_str());
  }
  
  #ifdef DEBUG
  ProfilerStop();
  #endif
  return 0;
}
    
    
    // gErrorIgnoreLevel = kBreak;
    
    //  m.simplex();
    //  m.migrad();
    //  m.improve();
    //  m.hesse();
    //  m.minos();
    
    // // RooFitResult *r = 0;
    // // "nll", "nll", "-2*log(@0)", RooArgSet(*pdf));
    // RooNLLVar *nll = new RooNLLVar("nll", "-2*log(L)", *sig, wdata, DataError(RooAbsData::SumW2));
    // RooMinimizer m(*nllt) ;
    // mini.optimizeConst(kTRUE) ;
    // mini.setEps(100);
    // mini.setMaxIterations(10000);
    // mini.minimize("Minuit2","Migrad") ;
    // mini.hesse() ;
    // r = mini.save() ;
  
   
    // // fit_results.push_back(bkg_function->fitTo(data, RooFit::SumW2Error(true), RooFit::Strategy(2), RooFit::NumCPU(4))); //RooFit::Extended(true),
    // std::cout << "vvvvvvvvvv" << std::endl;
    // fit_results.back()->printArgs(std::cout);
    // std::cout << "^^^^^^^^^^" << std::endl;
    // std::cout << "Result: " << fit_results.back()->printArgs() << std::endl;
    // gErrorIgnoreLevel = kInfo;
    // bkg_function->Print();
  
  


  // std::vector<std::string> params(argv, argv+argc);
  // D3PDReader::main::process_args(params);
  //
  // // If we are asking for plotting
  // if( D3PDReader::main::type == "plot" ) {
  //
  //   D3PDReader::AnalysisFactory analysisFactory;
  //   std::vector<D3PDReader::AnalysisPtr> analyses = analysisFactory.loadAnalysis( D3PDReader::main::optionSets[0] );
  //   if( analyses.size() == 0 ) { std::cerr << "D3PDReader::main::main: No analysis loaded." << std::endl; exit(1); }
  //
  //   D3PDReader::Plotter plotter;
  //   BOOST_FOREACH( D3PDReader::AnalysisPtr analysis, analyses ) { analysis->initialise( D3PDReader::main::optionSets[0] ); }
  //   analyses[0]->plot( analyses, plotter );
  //
  // // If we are asking for cutflow
  // } else if( D3PDReader::main::type == "cutflow" ) {
  //
  //   D3PDReader::AnalysisFactory analysisFactory;
  //   analysisFactory.parseOptions( D3PDReader::main::optionSets[0] );
  //   std::string analysisName( D3PDReader::main::optionSets.back().at(0).second );
  //   D3PDReader::AnalysisPtr cutflow( analysisFactory.loadAnalysis( analysisName, D3PDReader::main::inputFile ) );
  //   cutflow->printCutFlow();
  //
  // // We must be asking to run an analysis
  // } else {
  //
  //   // Set up and fill TChains
  //   D3PDReader::main::TChainPtr qcdChain( new TChain( "qcd" ) );
  //   D3PDReader::main::TChainPtr trigChain( new TChain( "qcdMeta/TrigConfTree" ) );
  //   D3PDReader::main::fillChains( D3PDReader::main::inputFile, qcdChain, trigChain );
  //
  //   // Set up Interpreter and TrigDecision objects
  //   const D3PDReader::InterpreterPtr ntuple( new D3PDReader::Interpreter( qcdChain.get(), D3PDReader::main::isMC ) );
  //   const D3PDReader::TrigDecisionPtr trigTool( new D3PD::TrigDecisionToolD3PD( ntuple->fChain, trigChain.get() ) );
  //   D3PDReader::EventLoop eventLoop( ntuple, trigTool );
  //
  //   // Iterate over analyses and add analysis to event loop
  //   for( std::vector<std::vector<std::pair<std::string, std::string> > >::const_iterator analysisOptionSet = D3PDReader::main::optionSets.begin();
  //        analysisOptionSet != D3PDReader::main::optionSets.end(); ++analysisOptionSet ) {
  //     eventLoop.addAnalysis( *analysisOptionSet, D3PDReader::main::isMC );
  //   }
  //
  //   // Run event loop
  //   eventLoop.initialise( D3PDReader::main::optionSets );
  //   eventLoop.execute( D3PDReader::main::nEvents );
  //   eventLoop.finalise();
  // }



// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //  Fill TChains from filelist
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void D3PDReader::main::fillChains( std::string inputFileList, TChainPtr inputChain, TChainPtr inputTrigChain ) {
//   int nFiles(0);
//   if( inputFileList == "" ) { return; }
//   std::ifstream inputFiles( inputFileList.c_str() );
//   if( inputFiles.is_open() ) { std::string inputFile;
//     while( inputFiles.good() ) {
//       getline( inputFiles, inputFile );
//       if( inputFile != "" && inputFile.at(0) != '#' ) {
//         int success_chain = inputChain->Add( inputFile.c_str() );
//         int success_trigChain = inputTrigChain->Add( inputFile.c_str() );
//         if( success_chain != 1 || success_trigChain != 1 ) { std::cerr << "Error adding file " << inputFile << "!" << std::endl; }
//         ++nFiles;
//       }
//     }
//     inputFiles.close();
//   }
//   std::cout << "Added " << nFiles << " files" << std::endl;
//   return;
// }

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// //  Process command line arguments
// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// void D3PDReader::main::process_args(std::vector<std::string> params) {
//
//   if( params.size() < 2) {
//     std::cout << params.at(0) << ": not enough arguments" << std::endl;
//     std::cout << " -a analysis type\n -in dataset list\n -out outputfile\n -nEvents nEvents" << std::endl;
//     exit(1);
//   }
//
//   for( std::vector<std::string>::const_iterator itParam = params.begin(); itParam != params.end(); ++itParam ) {
//     if( itParam->at(0) != '-' ) { continue; }
//
//     std::string delimiter( *itParam );
//     std::string argument( *(itParam+1) );
//     delimiter.erase(0, 1);
//
//     // If name is "a" then setup new analysis
//     if( !delimiter.compare(0, 1, "a") ) {
//       std::vector< std::pair<std::string, std::string> > optionSet;
//       optionSets.push_back( optionSet );
//     }
//
//     // If "in", "nEvents", "type" or "isMC" we only want one
//     if( !delimiter.compare(0, 2, "in") ) {
//       if( inputFile != "" ) { std::cout << "Overriding existing input file " << inputFile << " with " << argument << "!" << std::endl; }
//       inputFile = argument;
//     } else if( !delimiter.compare(0, 7, "nEvents") ) {
//       nEvents = boost::lexical_cast<int>(argument);
//     } else if( !delimiter.compare(0, 4, "type") ) {
//       if( type != "" ) { std::cout << "Overriding existing type " << type << " with " << argument << "!" << std::endl; }
//       type = argument;
//     } else if( !delimiter.compare(0, 4, "isMC") ) {
//       isMC = true;
//     // Otherwise we add to the per analysis list
//     } else {
//       if( !delimiter.compare(0, 3, "out") ) { FileHandler output(argument,FileHandler::RECREATE); }
//       std::pair<std::string, std::string> optionPair( delimiter, argument );
//       optionSets.back().push_back( optionPair );
//     }
//   }
// }
