#include <iostream>
#include <filesystem>
#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TGeoManager.h"
#include "TVirtualGeoTrack.h"
#include "TGeoNode.h"

void GetTracks(TString name = "Geometry.root")
{
    TFile* f = new TFile(name);
    TGeoManager* geom = (TGeoManager*)f->Get("Geometry");

    std::filesystem::create_directories("./Track_Data");

    for (int i = 0; i < geom->GetNtracks(); i++)
    {
        TVirtualGeoTrack* geo_track = geom->GetTrack(i);
        
        if((geo_track->GetLineColor() == kRed) || (geo_track->GetLineColor() == kBlue) || (geo_track->GetLineColor() == kMagenta) || (geo_track->GetLineColor() == kGreen))
        {
            TString file_name = "./Track_Data/";
            if(geo_track->GetLineColor() == kRed) file_name += "Red";
            if(geo_track->GetLineColor() == kBlue) file_name += "Blue";
            if(geo_track->GetLineColor() == kMagenta) file_name += "Magenta";
            if(geo_track->GetLineColor() == kGreen) file_name += "Green";
            file_name += "_" + std::to_string(i) + ".txt";

            std::cout << file_name << "\n";

            FILE *ff = fopen(file_name, "w");

            for (int j = 0; j < geo_track->GetNpoints(); j++)
            {
               Double_t xx, yy, zz, tt;

               geo_track->GetPoint(j, xx, yy, zz, tt);

               fprintf(ff, "%lg %lg %lg %lg \n", xx, yy, zz, tt);
	        }

            fclose(ff);	
        }
        geo_track = nullptr;
    }

    f->Close();
    f = nullptr;
    geom = nullptr;	
}