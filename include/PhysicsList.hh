#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <iostream>
#include <thread>
#include <omp.h>
#include <list>
#include <utility>
#include "TROOT.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVirtualGeoTrack.h"

#include "LinearInterpolation.hh"
#include "MaterialProperties.hh"
#include "Track.hh"

#define const_c 299792458.0
#define const_hr 1.054571817E-34 
#define const_e 1.60217662E-19
#define const_me 9.10938356E-31
#define const_k 8.9875517873681764E9
#define const_kB 1.380649E-23
#define const_Ry 2.179872361103542E-18 
#define unit_J 1.0
#define unit_K 1.0
#define unit_meV 1.60217662E-22
#define unit_eV 1.60217662E-19
#define unit_keV 1.60217662E-16
#define unit_MeV 1.60217662E-13
#define unit_m 1.0
#define unit_cm 1.0E-2
#define unit_mm 1.0E-3
#define unit_mkm 1.0E-6
#define unit_nm 1.0E-9
#define unit_A 1.0E-10
#define unit_s 1.0
#define unit_mks 1.0E-6
#define unit_ns 1.0E-9
#define unit_ps 1.0E-12
#define unit_fs 1.0E-15
#define unit_V 1.0
#define unit_kV 1.0E3
#define unit_MV 1.0E6
#define unit_T 1.0
#define unit_mT 1.0E-3
#define unit_mkT 1.0E-6

struct PeriodicBoundaries
{
  TVector3 center = TVector3(0,0,0);
  TVector3 edges = TVector3(1, 1, 1);

  PeriodicBoundaries(TVector3 v1, TVector3 v2)
  {
     center = v1;
     edges = v2;
  }
};

class PhysicsList
{
  public:
  PhysicsList();
 ~PhysicsList();
  void UpdateTimeCut();
  void SetEHGlobalTime(Double_t);
  Double_t GetEHGlobalTime();
  void SetTGeoTrack(Track*);
  void SetElectricFieldStrength(TVector3);
  TVector3 GetElectricFieldStrength();
  void SetMagneticFieldInduction(TVector3);
  TVector3 GetMagneticFieldInduction();
  void SetPeriodicBoundaries(TVector3, TVector3);
  PeriodicBoundaries* GetPeriodicBoundaries();
  bool ApplyPeriodicity(TVector3&);
  void DrawTrack(Track*);
  TVector3 RandomTVector3();
  Double_t RandomValue(LinearInterpolation*, Double_t, Double_t);
  Double_t GetDistanceBetweenTwoTracks(TVector3, TVector3, TVector3, TVector3);
  Int_t DiscreteRandomValue(std::vector<Double_t>);
  TString GetTGeoVolumeAtPoint(TVector3);
  void CutTrackLine(Track*);
  Double_t GenerateInelasticScatteringEnergy(Double_t, Double_t, TString);
  void GenerateVectorInElasticScattering(Double_t, Double_t, TVector3, TVector3&);
  bool GenerateVectorInInelasticScattering(Double_t, Double_t, Double_t, TVector3, TVector3&, TVector3&);
  void GenerateVectorInPhononScattering(Double_t, Double_t, Double_t, TVector3, TVector3&);
  void ElasticScattering(Track*, Int_t);
  void ThermalScattering(Track*, Int_t);
  void InelasticScattering(Track*, Int_t);
  void InelasticHoleScattering(Track*);
  void ElectronHoleInteraction(std::vector<Track*>);
  void RadiationDecayOfExciton(Track*);
  void DecayOfCenter(Track*, TString, TString, Double_t, bool);
  void ExcitonDissociation(Track*);
  void CaptureParticleByCenter(Track*, TString, TString, Double_t);
  void UpdateVelocity(Track*);
  void UpdateMass(Track*);
  Double_t CalculateMomentum(Double_t, Double_t);
  Double_t CalculateVelocity(Double_t, Double_t);
  void ApplySingleTrackProcesses(std::vector<Track*>);
  void WhichProcesse(Track*);
  void ApplyMultiTrackProcesses(std::vector<Track*>);
  void ParticleStep(std::vector<Track*>);
  Track* CreateElectron(TrackPoint, TString type = "", TString volume_name = "", MaterialProperties* mp = nullptr);
  Track* CreateHole(TrackPoint, TString type = "", TString volume_name = "", MaterialProperties* mp = nullptr);
  Track* CreateExciton(TrackPoint, TString type = "", TString volume_name = "", MaterialProperties* mp = nullptr, Double_t mex = 2.0*const_me);
  Track* CreatePhoton(TrackPoint, TString volume_name = "", MaterialProperties* mp = nullptr);
  Track* CreateCenter(TrackPoint, TString, TString, TString volume_name = "", MaterialProperties* mp = nullptr);

  private:
  Double_t TIME_CUT_eh, EH_GLOBAL_TIME;
  TVector3 Electric_Field_Strength, Magnetic_Field_Induction;
  PeriodicBoundaries* PeriodicCube = nullptr;
};

#endif