#include "RunManager.hh"

RunManager* RunManager::instance = nullptr;

RunManager::RunManager()
{
    VISIBLE = true;
    NUM_BEAMs = 1;
    Global_Iter = 1;

    std::cout << "Welcome to the LEPAM!\n";
}

RunManager::~RunManager()
{
  delete Detector_Construction;
  delete Run_Action;
  delete Track_Propagation;
  delete Primary_GeneratorAction;
  delete Physics_List;
  delete geometry;

  std::cout << "Goodbye :)\n";
}

RunManager* RunManager::getInstance()
{
    if (instance == nullptr)
    {
        instance = new RunManager();
    }

    return instance;
}

void RunManager::ReadParam(TString innam)
{
  FILE *f = fopen(innam, "r");
  if (!f)
  {
    std::cout << "RunManager: There isn't a input file! \n";
  }
  else
  {
    std::cout << "Config file: " << innam << "\n";

    char str[256], STR[256];
    while (str == fgets(str, 255, f))
    {
      memset(STR, 0, sizeof(STR));
      for (int i = 0; i < 256; ++i)
      {
        if ('/' == str[i] && '/' == str[i + 1])
          break;
        else if (0 == str[i])
          break;
        else
          STR[i] = str[i];
      }
  
      char what[128] = "", cn[128] = "";
      int nit = sscanf(STR, "%s %s", &what[0], &cn[0]);
      if (nit > 1)
      {
        double n = 0;
        int nit1 = sscanf(cn, "%lf", &n);
  
        if (0 == strcmp(what, "NUM_THREADs"))
        {
          if (nit > 1 && nit1 > 0)
          {
            Track_Propagation->SetNumThreads(n);
          }
        }
        if (0 == strcmp(what, "NUM_BEAMs"))
        {
          if (nit > 1 && nit1 > 0)
          {
            SetNumBeams(n);
          }
        }
        if (0 == strcmp(what, "VISIBLE"))
        {
          if (nit > 1 && nit1 > 0)
          {
            SetVisible(n);
          }
        }
      }
  
      memset(what, 0, 128);
      memset(cn, 0, 128);
      char cn2[128] = "";
      double n = 0;
      int nit2 = sscanf(STR, "%s %s %lf %s", &what[0], &cn[0], &n, &cn2[0]);
  
      if (nit2 > 3)
      {
        if (0 == strcmp(what, "SET_ENERGY"))
        {
          double units = unit_eV;
          if(0 == strcmp(cn2, "eV")) units = unit_eV;
          if(0 == strcmp(cn2, "keV")) units = unit_keV;
          if(0 == strcmp(cn2, "MeV")) units = unit_MeV;
          if(0 == strcmp(cn2, "meV")) units = unit_meV;
  
          List_Of_Initial_Particle_Energy.push_back(new InitialParticleEnergy(cn, n * units));
          std::cout << "Kinetic energy of " << cn << ": " << n << " " << cn2 << "\n";
        }
      }
    }
  fclose(f);
  }
}

void RunManager::SetVisible(bool vis)
{
    VISIBLE = vis;
    std::cout << "Track visualization: " << VISIBLE << "\n";
}

bool RunManager::GetVisible()
{
    return VISIBLE;
}

void RunManager::SetNumBeams(Int_t num)
{
    NUM_BEAMs = num;
    std::cout << "Total number of beams: " << NUM_BEAMs << "\n";
}

Int_t RunManager::GetNumBeams()
{
    return NUM_BEAMs;
}

Int_t RunManager::GetGlobalIter()
{
    return Global_Iter;
}

void RunManager::SetInitialParticleEnergy()
{   
    if(List_Of_Initial_Particle_Energy.size() == 0) return;
    
    for (auto it = List_Of_Initial_Particle_Energy.begin() ; it != List_Of_Initial_Particle_Energy.end(); ++it)
    {
      Primary_GeneratorAction->SetParticleEnergy((*it)->name, (*it)->energy);
    }
}

void RunManager::Initialization()
{
    Detector_Construction->BuildGeometry();
    Track_Propagation->ThreadInitialization();
}

void RunManager::Launch()
{ 
    std::cout << "\n";

    while(Global_Iter <= NUM_BEAMs)
    {
      Run_Action->Run();
      Track_Propagation->ClearTrackList();
      Global_Iter++;
    }
}

void RunManager::SaveGeom(TString outname)
{
    TFile* outfile = new TFile(outname,"recreate");
    TString topdir = gDirectory->GetPath();
    
    if(VISIBLE) geometry->DrawTracks();
    
    gDirectory->cd(topdir);
    outfile->WriteTObject(geometry);
    outfile->Close();

    std::cout << outname << " is saved! \n";

    delete outfile;    
}

void RunManager::SetDetectorConstruction(DetectorConstruction* local_obj)
{
    Detector_Construction = local_obj;
}

void RunManager::SetRunAction(RunAction* local_obj) 
{
    Run_Action = local_obj;
}

void RunManager::SetPrimaryGeneratorAction(PrimaryGeneratorAction* local_obj)
{
    Primary_GeneratorAction = local_obj;
}

void RunManager::SetTrackPropagation(TrackPropagation* local_obj)
{
    Track_Propagation = local_obj;
}

void RunManager::SetPhysicsList(PhysicsList* local_obj)
{
    Physics_List = local_obj;
}

void RunManager::SetTGeoManager(TGeoManager* local_obj)
{
    geometry = local_obj;
}

DetectorConstruction* RunManager::GetDetectorConstruction()
{
    return Detector_Construction;
}

RunAction* RunManager::GetRunAction() 
{
    return Run_Action;
}

PrimaryGeneratorAction* RunManager::GetPrimaryGeneratorAction()
{
    return Primary_GeneratorAction;
}

TrackPropagation* RunManager::GetTrackPropagation()
{
    return Track_Propagation;
}

PhysicsList* RunManager::GetPhysicsList()
{
    return Physics_List;
}

TGeoManager* RunManager::GetTGeoManager()
{
    return geometry;
}
