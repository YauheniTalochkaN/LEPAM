#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <iostream>
#include "TROOT.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoCompositeShape.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TProfile.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TVector3.h"

#include "MaterialProperties.hh"

struct PhysicalVolume
{
  TGeoVolume* volume = nullptr;
  MaterialProperties* properties = nullptr;

  PhysicalVolume(TGeoVolume* tgeo, MaterialProperties* prop)
  {
     volume = tgeo;
     properties = prop;
  }

  ~PhysicalVolume()
  {
     delete properties;
     volume = nullptr;
  }
};

class DetectorConstruction
{
  public:
  DetectorConstruction();
 ~DetectorConstruction();
  void BuildGeometry();
  void AddPhysicalVolume(TGeoVolume*, MaterialProperties*);
  MaterialProperties* GetMaterialProperties(TString);
  MaterialProperties* BuildProperties(TString);
  void ClearListOfPhysicalVolumes();
   
  private:
  std::vector <PhysicalVolume*> ListOfPhysicalVolumes;
};

#endif