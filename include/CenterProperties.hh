#ifndef CenterProperties_h
#define CenterProperties_h 1

#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

#include "LinearInterpolation.hh"

struct shell
{
   TString shell_name = "";
   Double_t shell_energy = 0;
   std::vector<std::pair<TString,Double_t>> radiation_times,
                                            nonradiation_times;

   shell(TString nam, Double_t energy, std::pair<TString,Double_t> rad_times[] = nullptr, Int_t num1 = 0, std::pair<TString,Double_t> nonrad_times[] = nullptr, Int_t num2 = 0) 
   {
      shell_name = nam;
      shell_energy = energy;

      for (int i = 0; i < num1; i++)
      {
         radiation_times.push_back(rad_times[i]);
      }

      for (int i = 0; i < num2; i++)
      {
         nonradiation_times.push_back(nonrad_times[i]);
      }
   }

   ~shell()
   {
      if(radiation_times.size() > 0) radiation_times.clear();
      if(nonradiation_times.size() > 0) nonradiation_times.clear();
   }
};

struct Container_ccr
{
   TString full_particle_name,
           center_name, center_type;
   LinearInterpolation* rate = nullptr;

   Container_ccr(TString pnam, TString cnam, TString ctype, LinearInterpolation* crate)
   {
      full_particle_name = pnam;
      center_name = cnam;
      center_type = ctype;
      rate = crate;
   }

   ~Container_ccr()
   {
      if(rate != nullptr) delete rate;
   }
};

class CenterProperties
{
  public:
  CenterProperties(TString);
 ~CenterProperties();

  void SetName(TString);
  void SetGroundState(TString);
  void SetLevel(shell*);
  void SetCaptureRate(Container_ccr*);

  TString GetName();
  TString GetGroundState();
  shell* GetLevel(TString);
  std::vector<Container_ccr*> GetCaptureRates(TString);

  private:
  TString Name = "", ground = "";
  std::vector<shell*> Level_List;
  std::vector<Container_ccr*> Capture_Rates;
};

#endif
