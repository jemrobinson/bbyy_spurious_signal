#pragma once
// STL
#include <iostream>
#include <string>
#include <vector>

class RooRealVar;

namespace SpuriousSignal {
  class ParameterSet {
  public:
    /**
     * ParameterSet constructor
     */
    ParameterSet(const std::string& name, std::vector<RooRealVar*> variables);

    void record_values();
    void restore_values();

    void read_from_file(const std::string& file_name);
    void write_to_file(const std::string& file_name);

    friend std::ostream& operator<<(std::ostream& os, const ParameterSet& p);

  private:
    const std::string m_name;
    std::vector<RooRealVar*> m_variables;
    std::vector<double> m_values;
  };
}
