#ifndef LinearInterpolation_h
#define LinearInterpolation_h 1

#include <iostream>
#include <fstream>
#include <vector>
#include "TROOT.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"

class LinearInterpolation
{
   public:
   LinearInterpolation(Double_t [], Double_t [], Int_t);
   LinearInterpolation(std::vector<Double_t>, std::vector<Double_t>);
   static LinearInterpolation* Initialization(TString, Double_t, Double_t);
   LinearInterpolation(FILE*, Double_t, Double_t);
  ~LinearInterpolation();
   Double_t GetValue(Double_t);
   Double_t GetXmin();
   Double_t GetXmax();
   Double_t GetYmin();
   Double_t GetYmax();
   Int_t    GetNBins();
   Double_t* GetXMass();
   Double_t* GetYMass();

   private:
   Double_t *MASSXX;
   Double_t *MASSYY;
   Int_t    NUMBINS;
   Double_t YMIN, YMAX;
};

#endif