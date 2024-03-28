#ifndef RunAction_h
#define RunAction_h 1

#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TH3I.h"
#include "TProfile.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"

class RunAction
{
  public:
  RunAction();
 ~RunAction();
  void Run();
  void EndOfCurrentStepOfRun();
  void EndOfRun();

  private:

};

#endif