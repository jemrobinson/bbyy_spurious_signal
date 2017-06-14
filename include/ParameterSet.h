#pragma once
#include "RooRealVar.h"
#include <string>

namespace SpuriousSignal {
  class ParameterSet {
  public:
    /**
     * ParameterSet constructor
     */
    ParameterSet(const std::string &name, std::vector<RooRealVar*> variables);
    
    void record_values();
    void restore_values();
    
  private:
    const std::string m_name;
    std::vector<RooRealVar*> m_variables;
    std::vector<double> m_values;
  };
}
