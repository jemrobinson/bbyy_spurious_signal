#pragma once
#include <string>
#include <vector>

namespace SpuriousSignal {
  struct PlotStyle {
    static void EnsureAtlasStyle();
    static std::vector<int> colours();
    static int colour(const std::string& fn_name);
    static std::string label(const std::string& fn_name);
    static std::vector<int> resonance_masses(const std::string& mass_category);
  };
}