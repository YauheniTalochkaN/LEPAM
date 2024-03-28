#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TGeoManager.h" 

void DrawGeo(TString name = "Geometry.root")
{
    TFile* f = new TFile(name);
    TGeoManager* geom = (TGeoManager*)f->Get("Geometry");

    geom->SetTopVisible(1);

    geom->GetTopVolume()->Draw();
    geom->DrawTracks();

    f->Close();
    f = nullptr;
    geom = nullptr;
}