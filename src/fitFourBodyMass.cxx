// // Local
// #include "D3PDReader.h"
// #include "AnalysisFactory.h"
// #include "EventLoop.h"
// #include "FileHandler.h"
// #include "Interpreter.h"
// #include "TrigDecisionToolD3PDWrapper.h"
// // STL
// #include <fstream>
// // BOOST
// #include <boost/foreach.hpp>
// #include <boost/lexical_cast.hpp>
// RooFit
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPlot.h"
// Google
#ifdef DEBUG
#include <google/profiler.h>
#endif

int main(int argc, char** argv) {
  #ifdef DEBUG
  ProfilerStart("out.prof");
  #endif

  RooRealVar x("x", "x", -10, 10);
  RooRealVar mean("mean", "mean of Gaussian", 0, -10, 10);
  RooRealVar sigma("sigma", "width of Gaussian", 3);

  RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);

  // Plot PDF
  RooPlot* xframe = x.frame();
  gauss.plotOn(xframe);
  xframe->Draw();

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

  #ifdef DEBUG
  ProfilerStop();
  #endif
  return 0;
}

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
