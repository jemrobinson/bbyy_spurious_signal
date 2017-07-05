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
    static std::vector<int> colours();

    /**
     * Colour corresponding to a named function
     */
    static int colour(const std::string& fn_name);

    /**
     * Label corresponding to a named function
     */
    static std::string label(const std::string& fn_name);

    /**
     * List of masses in a given category
     */
    static std::vector<int> resonance_masses(const std::string& mass_category);
  };
}