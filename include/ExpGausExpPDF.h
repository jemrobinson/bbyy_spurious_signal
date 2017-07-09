#pragma once

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "TObject.h"

class ExpGausExpPDF : public RooAbsPdf {

public:

  ExpGausExpPDF() {}
  ExpGausExpPDF(const char* name, const char* title, RooAbsReal& _m,
                RooAbsReal& _m0, RooAbsReal& _sigma, RooAbsReal& _kL,
                RooAbsReal& _kH);

  ExpGausExpPDF(const ExpGausExpPDF& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new ExpGausExpPDF(*this, newname); }
  inline virtual ~ExpGausExpPDF() {}

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
  ClassDef(ExpGausExpPDF, 2);
};
