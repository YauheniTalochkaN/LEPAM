#ifndef RunManager_h
#define RunManager_h 1

#include <iostream>
#include <ctime>
#include <chrono>
#include "TString.h"
#include "TProfile.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "TrackPropagation.hh"
#include "PrimaryGeneratorAction.hh"
#include "PhysicsList.hh"

#define geo_size_factor 1E12
#define geo_time_factor 1E18

struct InitialParticleEnergy
{
  TString name;
  Double_t energy;

  InitialParticleEnergy(TString nam, Double_t ek)
   {
      name = nam;
      energy = ek;
   }
};

class RunManager
{
  public:
  ~RunManager();
   static RunManager* getInstance();
   void ReadParam(TString);
   void Initialization();
   void Launch();
   void SaveGeom(TString);
   void SetInitialParticleEnergy();

   void SetDetectorConstruction(DetectorConstruction*);
   void SetRunAction(RunAction*);
   void SetPrimaryGeneratorAction(PrimaryGeneratorAction*);
   void SetTrackPropagation(TrackPropagation*);
   void SetPhysicsList(PhysicsList*);
   void SetTGeoManager(TGeoManager*);

   DetectorConstruction* GetDetectorConstruction();
   RunAction* GetRunAction();
   PrimaryGeneratorAction* GetPrimaryGeneratorAction();
   TrackPropagation* GetTrackPropagation();
   PhysicsList* GetPhysicsList();
   TGeoManager* GetTGeoManager();

   void SetVisible(bool);
   void SetNumBeams(Int_t);
   bool GetVisible();
   Int_t GetNumBeams();
   Int_t GetGlobalIter();

  private:
   RunManager();
   static RunManager* instance;

   RunAction* Run_Action;
   DetectorConstruction* Detector_Construction;
   TrackPropagation* Track_Propagation;
   PrimaryGeneratorAction* Primary_GeneratorAction;
   PhysicsList* Physics_List;
   TGeoManager* geometry;

   bool VISIBLE;
   Int_t NUM_BEAMs;
   Int_t Global_Iter;
   std::vector<InitialParticleEnergy*> List_Of_Initial_Particle_Energy;
};

#endif