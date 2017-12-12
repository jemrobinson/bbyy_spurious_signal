// Local
#include "ParameterSet.h"
#include "Logger.h"
// STL
#include <fstream>
#include <sstream>
// ROOT and RooFit
#include "RooRealVar.h"

namespace SpuriousSignal {
  /**
   * ParameterSet constructor
   */
  ParameterSet::ParameterSet(const std::string& name, std::vector<RooRealVar*> variables)
    : m_name(name)
    , m_variables(variables)
  {
    this->record_values();
  }

  void ParameterSet::record_values()
  {
    m_values.clear();

    for (auto variable : m_variables) {
      MSG_VERBOSE(m_name << ": recording values from: " << variable->getTitle() << ": " << variable->getVal());
      m_values.push_back(variable->getVal());
    }
  }

  void ParameterSet::restore_values()
  {
    for (unsigned int idx = 0; idx < m_variables.size(); ++idx) {
      MSG_VERBOSE(m_name << ": restoring variable: " << m_variables.at(idx)->getTitle() << " to " << m_values.at(idx));
      m_variables.at(idx)->setVal(m_values.at(idx));
    }
  }

  std::ostream& operator<<(std::ostream& os, const ParameterSet& p)
  {
    for (auto variable : p.m_variables) {
      os << variable->getTitle() << " [" << variable->getVal() << " +/- " << variable->getError() << "], ";
    }

    return os;
  }

  void ParameterSet::write_to_file(const std::string& file_name) {
    std::ofstream f_output;
    f_output.open(file_name, std::ios::app);
    for (auto variable : m_variables) {
      f_output << variable->GetName() << "\t" << variable->getVal() << "\t" << variable->getError() << std::endl;
    }
    f_output.close();
  }

  void ParameterSet::read_from_file(const std::string& file_name) {
    std::string _line;
    std::ifstream f_input(file_name);
    if (f_input.is_open()) {
      while (std::getline(f_input, _line)) {
        std::stringstream _ssline(_line);
        std::vector<std::string> _tokens;
        std::string _token;
        while (std::getline(_ssline, _token, '\t')) {
          _tokens.push_back(_token);
        }
        for (auto variable : m_variables) {
          if (std::string(variable->GetName()) == _tokens.at(0)) {
            variable->setVal(std::stod(_tokens.at(1)));
            variable->setError(std::stod(_tokens.at(2)));
            break;
          }
        }
      }
      f_input.close();
    }
  }

}