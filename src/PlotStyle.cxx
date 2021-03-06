#include "PlotStyle.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include <string>
#include <map>
#include <stdexcept>

namespace SpuriousSignal {
  void PlotStyle::EnsureAtlasStyle()
  {
    static TStyle* atlasStyle = 0;

    if (atlasStyle == 0) {
      atlasStyle = new TStyle("ATLAS", "ATLAS style");
      // use plain black on white colors
      Int_t icol = 0; // WHITE
      atlasStyle->SetFrameBorderMode(icol);
      atlasStyle->SetFrameFillColor(icol);
      atlasStyle->SetCanvasBorderMode(icol);
      atlasStyle->SetCanvasColor(icol);
      atlasStyle->SetPadBorderMode(icol);
      atlasStyle->SetPadColor(icol);
      atlasStyle->SetStatColor(icol);
      // set the paper & margin sizes
      atlasStyle->SetPaperSize(20, 26);
      // set margin sizes
      atlasStyle->SetPadTopMargin(0.05);
      atlasStyle->SetPadRightMargin(0.05);
      atlasStyle->SetPadBottomMargin(0.16);
      atlasStyle->SetPadLeftMargin(0.16);
      // set title offsets (for axis label)
      atlasStyle->SetTitleXOffset(1.6); // 1.4
      atlasStyle->SetTitleYOffset(1.6); // 1.4
      // use large fonts
      Int_t font = 42; // Helvetica scalable
      Double_t tsize = 0.04; //0.05;
      // Int_t font = 43; // Helvetica scalable
      // Double_t tsize = 25; //30;
      atlasStyle->SetTextFont(font);
      atlasStyle->SetTextSize(tsize);
      atlasStyle->SetLabelFont(font, "x");
      atlasStyle->SetTitleFont(font, "x");
      atlasStyle->SetLabelFont(font, "y");
      atlasStyle->SetTitleFont(font, "y");
      atlasStyle->SetLabelFont(font, "z");
      atlasStyle->SetTitleFont(font, "z");
      atlasStyle->SetLabelSize(tsize, "x");
      atlasStyle->SetTitleSize(tsize, "x");
      atlasStyle->SetLabelSize(tsize, "y");
      atlasStyle->SetTitleSize(tsize, "y");
      atlasStyle->SetLabelSize(tsize, "z");
      atlasStyle->SetTitleSize(tsize, "z");
      // use bold lines and markers
      atlasStyle->SetMarkerStyle(20);
      atlasStyle->SetMarkerSize(0.8);
      atlasStyle->SetHistLineWidth(2.);
      atlasStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes
      // get rid of X error bars (as recommended in ATLAS figure guidelines)
      atlasStyle->SetErrorX(0.0001);
      // get rid of error bar caps
      atlasStyle->SetEndErrorSize(0.);
      // do not display any of the standard histogram decorations
      atlasStyle->SetOptTitle(0);
      atlasStyle->SetOptStat(0);
      atlasStyle->SetOptFit(0);
      // put tick marks on top and RHS of plots
      atlasStyle->SetPadTickX(1);
      atlasStyle->SetPadTickY(1);
      // now set the style and end
      gROOT->SetStyle("ATLAS");
      gROOT->ForceStyle();
    }
  }

  std::vector<int> PlotStyle::colours(const int& nColours)
  {
    if (nColours == 6) {
      // Six qualitative colours
      return std::vector<int>({TColor::GetColor("#1b9e77"), TColor::GetColor("#d95f02"), TColor::GetColor("#7570b3"), TColor::GetColor("#e7298a"), TColor::GetColor("#66a61e"), TColor::GetColor("#e6ab02")});
    } else if (nColours == 10) {
      // Ten qualitative colours
      return std::vector<int>({TColor::GetColor("#c00026"), TColor::GetColor("#a00026"), TColor::GetColor("#800026"), TColor::GetColor("#bd0026"), TColor::GetColor("#e31a1c"), TColor::GetColor("#fc4e2a"), TColor::GetColor("#fd8d3c"), TColor::GetColor("#feb24c"), TColor::GetColor("#fed976"), TColor::GetColor("#ffeda0")});
    }
    return std::vector<int>();
  }

  int PlotStyle::colour(const std::string& fn_name)
  {
    std::map<std::string, int> _map = {{"novosibirsk", TColor::GetColor("#e7298a")}, {"modified_gamma", TColor::GetColor("#1b9e77")}, {"modified_landau", TColor::GetColor("#7570b3")},
                                       {"exppoly1", TColor::GetColor("#66c2a5")}, {"exppoly2", TColor::GetColor("#fc8d62")}, {"invpoly2", TColor::GetColor("#8da0cb")}, {"invpoly3", TColor::GetColor("#e78ac3")}, {"powerlaw", TColor::GetColor("#a6d854")}};
    return _map[fn_name];
  }

  std::string PlotStyle::label(const std::string& fn_name)
  {
    std::map<std::string, std::string> _map = {{"novosibirsk", "Novosibirsk"}, {"modified_gamma", "Modified Gamma"}, {"modified_landau", "Modified Landau"},
                                               {"exppoly1", "Exp. (x)"}, {"exppoly2", "Exp. (x^{2})"}, {"invpoly2", "Inv. poly. (x^{-2})"}, {"invpoly3", "Inv. poly. (x^{-3})"}, {"powerlaw", "Power-law"}};
    return _map[fn_name];
  }

  std::pair<int, int> PlotStyle::mass_range(const std::string& mass_category) {
    if (mass_category == "low") {
      return std::make_pair<int, int>(245, 610);
    } else if (mass_category == "high") {
      return std::make_pair<int, int>(335, 1140);
    }
    throw std::invalid_argument("Did not recognise mass category '" + mass_category + "'");
  }

  std::vector<int> PlotStyle::resonance_masses(const std::string& mass_category)
  {
    if (mass_category == "low") {
      return std::vector<int>({260, 275, 300, 325, 350, 400, 450, 500});
    } else if (mass_category == "high") {
      return std::vector<int>({400, 450, 500, 750, 1000});
    }
    throw std::invalid_argument("Did not recognise mass category '" + mass_category + "'");
  }
}
