#ifndef Track_h
#define Track_h 1

#include <iostream>
#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TString.h"
#include "TVirtualGeoTrack.h"

#include "MaterialProperties.hh"
#include "PrimaryGeneratorAction.hh"

struct TrackPoint
{
  TVector3  momentum_direction = TVector3(0,0,1);
  TVector3  point = TVector3(0,0,0);
  Double_t  energy = 0;
  Double_t  kinetic_energy = 0;
  Double_t  Ecoulomb_energy = 0;
  Double_t  velocity = 0;
  Double_t  time = 0;
  Double_t  time_step = 1.0E-15; 
  TString   volume = "";
  MaterialProperties *matprop = nullptr;

  TrackPoint& operator = (const TrackPoint& tp) 
  {
      momentum_direction = tp.momentum_direction;
      point = tp.point;
      kinetic_energy = tp.kinetic_energy;
      Ecoulomb_energy = tp.Ecoulomb_energy;
      energy = tp.energy;
      velocity = tp.velocity;
      time = tp.time;
      time_step = tp.time_step;
      volume = tp.volume;
      matprop = tp.matprop;
      
      return *this;
  }
};

class PhysicsList;

class Track
{
  friend PhysicsList;

  public:
    Track();
    ~Track();
    TString GetParticleName();
    TString GetParticleType();
    Int_t GetTrackID();
    TVector3 GetInitialPoint();
    Double_t GetCreationTime();
    TString GetCreatorProcess();
    Double_t GetMass();
    Double_t GetCharge();
    TrackPoint GetPreStepPoint();
    TrackPoint GetPostStepPoint();
    bool GetKilled();

  private:
    TString particle_name;
    TString particle_type;
    Int_t track_id;
    TVector3 initial_point;
    Double_t creation_time; 
    TString creator_process;
    Double_t mass;
    Double_t charge;
    TrackPoint pre_step_point;
    TrackPoint post_step_point;
    bool killed;
    bool should_be_killed;
    bool visible;
    TVirtualGeoTrack* geo_track = nullptr;

    friend void PrimaryGeneratorAction::SetParticleEnergy(TString, Double_t);
};

#endif

