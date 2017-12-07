/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 
#include "RooFit.h"
#include "Riostream.h" 
#include "RooRelBW.h" 
#include "RooAbsReal.h" 
#include "RooRealVar.h"
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 
using namespace std;

ClassImp(RooRelBW) 

 RooRelBW::RooRelBW(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _mass,
                        RooAbsReal& _width) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   mass("mass","mass",this,_mass),
   width("width","width",this,_width)
 { 
 } 


 RooRelBW::RooRelBW(const RooRelBW& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   mass("mass",this,other.mass),
   width("width",this,other.width)
 { 
 } 



 Double_t RooRelBW::evaluate() const 
 { 
   //expression for function BW with gamma=f(width,m_pi,m_K)
   const Double_t mkaon2 = 0.4937*0.4937; //GeV
   const Double_t mpion2 = 0.1396*0.1396; //GeV
   Double_t mass2 = (mass*mass);
   Double_t mass4 = (mass2*mass2);
   Double_t factx4 = 1./(x*x*x*x);
   Double_t diffnum = (x*x)-mpion2-mkaon2;
   Double_t prod = 4*(mpion2*mkaon2);
   Double_t diffden = (mass2-mpion2-mkaon2);   
   Double_t ratio = (diffnum*diffnum-prod)/(diffden*diffden-prod);
   Double_t gamma = (width*mass4*factx4)*(TMath::Power(ratio, 1.5));
   Double_t arg = (x*x)-mass2;
   Double_t arg2 = (arg*arg);
   Double_t mw = gamma*mass;
   return (x*mw)/(arg2+(mw*mw));
 } 