#ifndef TrackPropagation_h
#define TrackPropagation_h 1

#include <iostream>
#include <vector>
#include <list>
#include <thread>
#include <omp.h>
#include "TROOT.h"
#include "TGeoManager.h"
#include "TMath.h"
#include "TVirtualGeoTrack.h"

#include "PhysicsList.hh"

class TrackPropagation
{
  public:
  TrackPropagation();
 ~TrackPropagation();
  Int_t MoveTracks();
  Track* GetTrack(Int_t);
  Int_t SizeofTrackList();
  void ClearTrackList();
  void LockMutex();
  void UnLockMutex();
  void SetNumThreads(Int_t);
  Int_t GetNumThreads();
  void ThreadInitialization();
  void ClearNewTrackList();
  void CopyTracksFromNewTrackList();
  void AddNewTrack(Track*);

  private:
  std::mutex MTX;
  Int_t NUM_THREADs;
  std::vector<Track*> ListOfTracks;
  std::vector<Track*> *NewTrackList;
};

#endif