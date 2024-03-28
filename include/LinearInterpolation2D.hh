#ifndef LinearInterpolation2D_h
#define LinearInterpolation2D_h 1

#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"

#include "LinearInterpolation.hh"

struct projectionX
{
   Double_t X = 0;
   LinearInterpolation *slice;

   projectionX(Double_t data_1, LinearInterpolation* data_2)
   {
      X = data_1;
      slice = data_2;
   }

   ~projectionX()
   {
      delete slice;
   }
};

class LinearInterpolation2D
{
   public:
   static LinearInterpolation2D* Initialization(TString, Double_t, Double_t, Double_t);
   LinearInterpolation2D(FILE*, Double_t, Double_t, Double_t);
   LinearInterpolation2D(std::vector<projectionX*>);
  ~LinearInterpolation2D();
   Double_t GetValue(Double_t, Double_t);
   LinearInterpolation* GetXProjection(Double_t);
   LinearInterpolation* GetYProjection(Double_t);
   Double_t GetXmin();
   Double_t GetXmax();
   Double_t GetYmin();
   Double_t GetYmax();
   Double_t GetZmin();
   Double_t GetZmax();

   private:
   std::vector<projectionX*> LinInter1DList;
};

#endif