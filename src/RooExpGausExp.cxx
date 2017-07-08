#include "RooExpGausExp.h"
#include "Math/ProbFuncMathCore.h"
#include "TMath.h"

ClassImp(RooExpGausExp)

//_____________________________________________________________________________
RooExpGausExp::RooExpGausExp(const char* name, const char* title,
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
RooExpGausExp::RooExpGausExp(const RooExpGausExp& other, const char* name) :
  RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma),
  kLo("kLo", this, other.kLo),
  kHi("kHi", this, other.kHi)
{
}


//_____________________________________________________________________________
Double_t RooExpGausExp::evaluate() const
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
Int_t RooExpGausExp::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
  if (matchArgs(allVars, analVars, m))
    return 1;

  return 0;
}


//_____________________________________________________________________________
Double_t RooExpGausExp::analyticalIntegral(Int_t code, const char* rangeName) const
{
  assert(code == 1);
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
double RooExpGausExp::gaussianIntegral(double tmin, double tmax) const
{
  return sqrt(TMath::TwoPi()) * (ROOT::Math::gaussian_cdf(tmax) - ROOT::Math::gaussian_cdf(tmin));
}

//_____________________________________________________________________________
double RooExpGausExp::exponentialIntegral(double tmin, double tmax, double k) const
{
  return (TMath::Exp(0.5 * k * k + k * tmax) - TMath::Exp(0.5 * k * k + k * tmin)) / k;
}
