#ifndef LinearInterpolation3D_h
#define LinearInterpolation3D_h 1

#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TMath.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"

#include "LinearInterpolation2D.hh"

struct sliceYZ
{
   Double_t X = 0;
   LinearInterpolation2D *surface;

   sliceYZ(Double_t data_1, LinearInterpolation2D* data_2)
   {
      X = data_1;
      surface = data_2;
   }

   ~sliceYZ()
   {
      delete surface;
   }
};

class LinearInterpolation3D
{
   public:
   static LinearInterpolation3D* Initialization(TString, Double_t, Double_t, Double_t, Double_t);
   LinearInterpolation3D(FILE*, Double_t, Double_t, Double_t, Double_t);
  ~LinearInterpolation3D();
   Double_t GetValue(Double_t, Double_t, Double_t);
   
   Double_t GetXmin();
   Double_t GetXmax();
   /*Double_t GetYmin();
   Double_t GetYmax();
   Double_t GetZmin();
   Double_t GetZmax();
   Double_t GetPmin();
   Double_t GetPmax();*/

   private:
   std::vector<sliceYZ*> LinInter2DList;
};

#endif