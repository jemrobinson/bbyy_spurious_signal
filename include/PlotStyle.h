#pragma once
// STL
#include <string>
#include <vector>

namespace SpuriousSignal {
  struct PlotStyle {
    /**
     * Apply ATLAS plot style if not already done
     */
    static void EnsureAtlasStyle();

    /**
     * Diverging colours
     */
    static std::vector<int> colours(const int& nColours = 10);

    /**
     * Colour corresponding to a named function
     */
    static int colour(const std::string& fn_name);

    /**
     * Label corresponding to a named function
     */
    static std::string label(const std::string& fn_name);

    /**
     * Pair of edges giving mass range
     */
    static std::pair<int, int> mass_range(const std::string& mass_category);

    /**
     * List of masses in a given category
     */
    static std::vector<int> resonance_masses(const std::string& mass_category);

    template <typename T>
    static std::string to_string(const T& input, const int& ndp) {
      std::string s_out = std::to_string(input);
      int nCharacters(ndp < 1 ? ndp : 1 + ndp);
      return s_out.substr(0, s_out.find(".") + nCharacters);
    }

  };
}