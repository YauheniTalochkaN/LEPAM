#ifndef MaterialProperties_h
#define MaterialProperties_h 1

#include <iostream>
#include <vector>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"

#include "LinearInterpolation.hh"
#include "LinearInterpolation2D.hh"
#include "CenterProperties.hh"

struct Container
{
   LinearInterpolation* func = nullptr;
   TString full_particle_name = "";
   
   Container(TString name, LinearInterpolation* data)
   {
      func = data;
      full_particle_name = name;
   }

   ~Container()
   {
      delete func;
   }
};

struct Container_imfp
{
   LinearInterpolation* imfp = nullptr;
   TString full_particle_name = "";
   LinearInterpolation* els = nullptr;

   Container_imfp(TString name, LinearInterpolation* data_1, LinearInterpolation* data_2)
   {
      imfp = data_1;
      full_particle_name = name;
      els = data_2;
   }

   ~Container_imfp()
   {
      delete imfp;
      delete els;
   }
};

struct Container_emfp
{
   LinearInterpolation* ait = nullptr;
   TString full_particle_name = "";
   LinearInterpolation2D* mls = nullptr;

   Container_emfp(TString name, LinearInterpolation* data_1, LinearInterpolation2D* data_2)
   {
      ait = data_1;
      full_particle_name = name;
      mls = data_2;
   }

   ~Container_emfp()
   {
      delete ait;
      delete mls;    
   }
};

struct Container_ait
{
   LinearInterpolation* ait = nullptr;
   TString full_particle_name = "";
   LinearInterpolation2D* mls = nullptr;
   LinearInterpolation* displow = nullptr;

   Container_ait(TString name, LinearInterpolation* data_1, LinearInterpolation2D* data_2, LinearInterpolation* data_3)
   {
      ait = data_1;
      full_particle_name = name;
      mls = data_2;
      displow = data_3;
   }

   ~Container_ait()
   {
      delete ait;
      delete mls;
      delete displow;
   }
};

struct Container_exadt
{
   LinearInterpolation* ait = nullptr;
   LinearInterpolation2D* des = nullptr;

   Container_exadt(LinearInterpolation* data_1, LinearInterpolation2D* data_2)
   {
      ait = data_1;
      des = data_2;
   }

   ~Container_exadt()
   {
      delete ait;
      delete des;
   }
};


class MaterialProperties
{
  public:
  MaterialProperties();
 ~MaterialProperties();

  void SetBandGap(Double_t);
  void SetDielectricPermittivity(Double_t);
  void SetEffectiveMass(Container*);
  void SetUnitCellVolume(Double_t);
  void SetDensityOfStates(LinearInterpolation*);
  void SetGroupVelocity(Container*);
  void SetInelasticMeanFreePath(Container_imfp*);
  void SetElasticMeanFreePath(Container_emfp*);
  void SetThermalAverageInteractionTime(Container_ait*);
  void SetRadiationDecayTimeOfExciton(LinearInterpolation*);
  void SetElectronHoleRadiationRecombinationTime(LinearInterpolation2D*);
  void SetExcitonDissociationTime(Container_exadt*);
  void SetCenterProperties(CenterProperties*);

  Double_t GetBandGap();
  Double_t GetDielectricPermittivity();
  Double_t GetUnitCellVolume();
  LinearInterpolation* GetEffectiveMass(TString);
  LinearInterpolation* GetDensityOfStates();
  LinearInterpolation* GetGroupVelocity(TString);
  std::vector<Container_imfp*> GetListOfInelasticMeanFreePath(TString);
  std::vector<Container_emfp*> GetListOfElasticAverageInteractionTime(TString);
  std::vector<Container_ait*> GetListOfThermalAverageInteractionTime(TString);
  LinearInterpolation* GetRadiationDecayTimeOfExciton();
  LinearInterpolation2D* GetElectronHoleRadiationRecombinationTime();
  Container_exadt* GetExcitonDissociationTime();
  std::vector<CenterProperties*> GetCenterList();
  CenterProperties* GetCenter(TString);

  private:
  Double_t Band_gap = 0;
  Double_t Unit_Cell_Volume = 0;
  Double_t Dielectric_Permittivity = 1;
  std::vector<Container*> Effective_Mass;
  LinearInterpolation* Density_Of_States = nullptr;
  std::vector<Container*> Group_Velocity;
  std::vector<Container_imfp*> Inelastic_Mean_Free_Paths;
  std::vector<Container_emfp*> Elastic_Mean_Free_Path;
  std::vector<Container_ait*> Thermal_Average_Interaction_Time; 
  LinearInterpolation* Radiation_Decay_Time_Of_Exciton = nullptr;
  LinearInterpolation2D* Electron_Hole_Radiation_Recombination_Time = nullptr;
  Container_exadt* Exciton_Dissociation_Time = nullptr;
  std::vector<CenterProperties*> List_Of_Centers;
};

#endif