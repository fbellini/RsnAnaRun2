// Flatte' Formula - fit f0(980)

#include "Riostream.h"
#include "RooFit.h"
#include "RooFlatte.h"
#include "RooAbsCategory.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooComplex.h"
#include "TMath.h"

using namespace std;

ClassImp(RooFlatte)

    RooFlatte::RooFlatte(const char* name, const char* title,
        RooAbsReal& _x, RooAbsReal& _mean,
        RooAbsReal& _g0, RooAbsReal& _m0a, RooAbsReal& _m0b,
        RooAbsReal& _g1, RooAbsReal& _m1a, RooAbsReal& _m1b)
    : RooAbsPdf(name, title)
    , x("x", "Dependent", this, _x)
    , mean("mean", "Mean", this, _mean)
    , g0("g0", "Channel 1 coupling", this, _g0)
    , m0a("m0a", "M1 in channel 1", this, _m0a)
    , m0b("m0b", "M2 in channel 1", this, _m0b)
    , g1("g1", "Channel 2 coupling", this, _g1)
    , m1a("m1a", "M1 in channel 2", this, _m1a)
    , m1b("m1b", "M2 in channel 2", this, _m1b)
{
}

RooFlatte::RooFlatte(const RooFlatte& other,
    const char* name)
    : RooAbsPdf(other, name)
    , x("x", this, other.x)
    , mean("mean", this, other.mean)
    , g0("g0", this, other.g0)
    , m0a("m0a", this, other.m0a)
    , m0b("m0b", this, other.m0b)
    , g1("g1", this, other.g1)
    , m1a("m1a", this, other.m1a)
    , m1b("m1b", this, other.m1b)
{
}

Double_t RooFlatte::evaluate() const
{
  //Flatté amplitude Flatté-like fit (assuming that the couplings are to KK and ππ)

  if (g0 < 0 || g1 < 0) {
    return (0);
  }

  Double_t s = x * x;

  // Energy, centre of mass p^2 of first channel
  Double_t E0a = 0.5 * (s + m0a * m0a - m0b * m0b) / x;
  Double_t qSq0 = E0a * E0a - m0a * m0a;

  // Energy, centre of mass p^2 of second channel
  Double_t E1a = 0.5 * (s + m1a * m1a - m1b * m1b) / x;
  Double_t qSq1 = E1a * E1a - m1a * m1a;

  RooComplex gamma0 = (qSq0 > 0) ? RooComplex(g0 * sqrt(qSq0), 0) : RooComplex(0, g0 * sqrt(-qSq0));

  RooComplex gamma1 = (qSq1 > 0) ? RooComplex(g1 * sqrt(qSq1), 0) : RooComplex(0, g1 * sqrt(-qSq1));

  RooComplex gamma = gamma0 + gamma1;

  RooComplex partB = RooComplex(0.0, 2 * mean / x) * gamma;
  RooComplex partA(mean * mean - s, 0);

  RooComplex denom = partA - partB;

  //RooComplex T(mean*sqrt(g0*g1),0);
  RooComplex T(1, 0);
  T = T / denom;

  return (T.abs2());
}
