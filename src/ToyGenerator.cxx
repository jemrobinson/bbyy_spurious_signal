#include "ToyGenerator.h"
#include "RooRandom.h"
#include "RooDataSet.h"

namespace SpuriousSignal {

  ToyGenerator *ToyGenerator::s_instance = 0;

  ToyGenerator::ToyGenerator() {
    RooRandom::randomGenerator()->SetSeed(20171011);
  }

  ToyGenerator* ToyGenerator::instance() {
    if (!s_instance) {
      s_instance = new ToyGenerator;
    }
    return s_instance;
  }

  RooDataSet* ToyGenerator::generate(RooDataSet* input_dataset) {
    // Construct the output dataset
    if (input_dataset->numEntries() == 0) { return 0; }
    RooDataSet* output_dataset = new RooDataSet(input_dataset->GetName(), input_dataset->GetTitle(), *input_dataset->get(0));
    const RooArgSet *event(0);

    // Iterate over input events, selected if they pass
    for (int idx = 0; idx < input_dataset->numEntries(); ++idx) {
      event = input_dataset->get(idx);
      // Keep event if its weight is more than a random number
      if (RooRandom::uniform() < input_dataset->weight() ) {
        output_dataset->add(*event);
      }
    }
    return output_dataset;
  }
}
