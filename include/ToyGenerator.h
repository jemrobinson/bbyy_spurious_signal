#pragma once

class RooDataSet;

namespace SpuriousSignal {
  class ToyGenerator {
  public:
    /**
     * Singleton accessor
     */
    static ToyGenerator* instance();

    /**
     * Generate a dataset given an input dataset and weight to select on
     */
    RooDataSet* generate(RooDataSet* input_dataset);

  private:
    /**
     * Default constructor
     */
    ToyGenerator();

    /**
     * Static instance
     */
    static ToyGenerator *s_instance;

  };
}