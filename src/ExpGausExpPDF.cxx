#include "ExpGausExpPDF.h"
#include "Math/ProbFuncMathCore.h"
#include "TMath.h"

ClassImp(ExpGausExpPDF)

//_____________________________________________________________________________
ExpGausExpPDF::ExpGausExpPDF(const char* name, const char* title,
                 RooAbsReal& _m, RooAbsReal& _m0,
                 RooAbsReal& _sigma, RooAbsReal& _kLo,
                 RooAbsReal& _kHi) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  kLo("kLo", "Low-side k", this, _kLo),
  kHi("kHi", "High-side k", this, _kHi)
{
}


//_____________________________________________________________________________
ExpGausExpPDF::ExpGausExpPDF(const ExpGausExpPDF& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma),
  kLo("kLo", this, other.kLo),
  kHi("kHi", this, other.kHi)
{
}


//_____________________________________________________________________________
Double_t ExpGausExpPDF::evaluate() const
{
  Double_t t = (m - m0) / sigma;

  if (t < -kLo) {
    return TMath::Exp(0.5 * kLo * kLo + kLo * t);
  } else if (t > kHi) {
    return TMath::Exp(0.5 * kHi * kHi - kHi * t);
  }

  return TMath::Exp(-0.5 * t * t);
}


//_____________________________________________________________________________
Int_t ExpGausExpPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, m)) {
    return 1;
  }
  return 0;
}


//_____________________________________________________________________________
Double_t ExpGausExpPDF::analyticalIntegral(Int_t /*code*/, const char* rangeName) const
{
  double result = 0;
  double sig = fabs((Double_t)sigma);
  double tmin = (m.min(rangeName) - m0) / sig;
  double tmax = (m.max(rangeName) - m0) / sig;

  if (tmin < -kLo) {
    result += exponentialIntegral(tmin, TMath::Min(tmax, -kLo), kLo);
  }

  if (tmin < kHi && tmax > -kLo) {
    result += gaussianIntegral(TMath::Max(tmin, -kLo), TMath::Min(tmax, kHi));
  }

  if (tmax > kHi) {
    result += exponentialIntegral(-tmax, TMath::Min(-tmin, -kHi), kHi);
  }

  return sig * result;
}

//_____________________________________________________________________________
double ExpGausExpPDF::gaussianIntegral(double tmin, double tmax) const
{
  return sqrt(TMath::TwoPi()) * (ROOT::Math::gaussian_cdf(tmax) - ROOT::Math::gaussian_cdf(tmin));
}

//_____________________________________________________________________________
double ExpGausExpPDF::exponentialIntegral(double tmin, double tmax, double k) const
{
  return (TMath::Exp(0.5 * k * k + k * tmax) - TMath::Exp(0.5 * k * k + k * tmin)) / k;
}
