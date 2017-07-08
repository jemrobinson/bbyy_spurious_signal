#pragma once

#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "TObject.h"

class RooDSCB : public RooAbsPdf {

public:

  RooDSCB() {}
  RooDSCB(const char* name, const char* title, RooAbsReal& _m,
          RooAbsReal& _m0, RooAbsReal& _sigma, RooAbsReal& _alphaLo,
          RooAbsReal& _nLo, RooAbsReal& _alphaHi, RooAbsReal& _nHi);

  RooDSCB(const RooDSCB& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const
  {
    return new RooDSCB(*this, newname);
  }
  inline virtual ~RooDSCB() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

  double gaussianIntegral(double tmin, double tmax) const;
  double powerLawIntegral(double tmin, double tmax, double alpha, double n) const;

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigma;
  RooRealProxy alphaLo;
  RooRealProxy nLo;
  RooRealProxy alphaHi;
  RooRealProxy nHi;

  Double_t evaluate() const;

private:
  ClassDef(RooDSCB, 1);
};
