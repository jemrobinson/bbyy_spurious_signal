#pragma once

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "TObject.h"

class RooExpGausExp : public RooAbsPdf {

public:

  RooExpGausExp() {}
  RooExpGausExp(const char* name, const char* title, RooAbsReal& _m,
                RooAbsReal& _m0, RooAbsReal& _sigma, RooAbsReal& _kL,
                RooAbsReal& _kH);

  RooExpGausExp(const RooExpGausExp& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooExpGausExp(*this, newname); }
  inline virtual ~RooExpGausExp() {}

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  double gaussianIntegral(double tmin, double tmax) const;
  double exponentialIntegral(double tmin, double tmax, double k) const;

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy kLo;
  RooRealProxy kHi;

  Double_t evaluate() const;

private:
  ClassDef(RooExpGausExp, 2);
};
