#include "ParameterSet.h"
#include "Logger.h"

namespace SpuriousSignal {
  /**
   * ParameterSet constructor
   */
  ParameterSet::ParameterSet(const std::string &name, std::vector<RooRealVar*> variables)
    : m_name(name)
    , m_variables(variables) {
    this->record_values();
  }

  void ParameterSet::record_values() {
    m_values.clear();
    for (auto variable : m_variables) {
      MSG_VERBOSE(m_name << ": recording values from: " << variable->getTitle() << ": " << variable->getVal());
      m_values.push_back(variable->getVal());
    }
  }
  
  void ParameterSet::restore_values() {
    for (unsigned int idx = 0; idx < m_variables.size(); ++idx) {
      MSG_VERBOSE(m_name << ": restoring variable: " << m_variables.at(idx)->getTitle() << " to " << m_values.at(idx));
      m_variables.at(idx)->setVal(m_values.at(idx));
    }
  }
  
  std::ostream& operator<<(std::ostream& os, const ParameterSet& m) {
    for (auto variable : m.m_variables) {
      os << variable->getTitle() << " [" << variable->getVal() << " +/- " << variable->getError() << "], ";
    }
    return os;
  }
}