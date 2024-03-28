#include "RunAction.hh"
#include "RunManager.hh"

RunAction::RunAction()
{

}

RunAction::~RunAction()
{
  
}

void RunAction::Run()
{   
    std::cout << "Beam number: " << RunManager::getInstance()->GetGlobalIter() << "\r" << std::flush;

    RunManager::getInstance()->GetPrimaryGeneratorAction()->GeneratePrimaries();
    RunManager::getInstance()->GetTrackPropagation()->CopyTracksFromNewTrackList();
    RunManager::getInstance()->SetInitialParticleEnergy();
    RunManager::getInstance()->GetPhysicsList()->SetEHGlobalTime(0);

    Int_t live_tracks = 1;
    while (live_tracks > 0)
    {
        live_tracks = RunManager::getInstance()->GetTrackPropagation()->MoveTracks();

        EndOfCurrentStepOfRun();
    }

    EndOfRun();
}

void RunAction::EndOfCurrentStepOfRun()
{

}

void RunAction::EndOfRun()
{
    TFile *fout = 0;
    fout = new TFile("Track_Data.root", "update");

    char topdir[256];
    memset(topdir,0,sizeof(topdir));
	  strncpy(topdir, gDirectory->GetPath(), 240);

    Double_t xx, yy, zz, cxx, cyy, czz, ctime, En;
    Int_t id, numbeam;
    Int_t particle;
    
    TTree *tout = (TTree*)fout->Get("Track_Data");
    if(tout == nullptr) 
    {
      tout = new TTree("Track_Data", "Track_Data");
      tout->Branch("beam",&numbeam,"beam/I");
      tout->Branch("particle",&particle,"particle/I");
      tout->Branch("id",&id,"id/I");
      tout->Branch("Energy",&En,"Energy/D");
      tout->Branch("x",&xx,"x/D");
      tout->Branch("y",&yy,"y/D");
      tout->Branch("z",&zz,"z/D");
      tout->Branch("cx",&cxx,"cx/D");
      tout->Branch("cy",&cyy,"cy/D");
      tout->Branch("cz",&czz,"cz/D");
      tout->Branch("ctime",&ctime,"ctime/D");
    }
    else
    {
      tout->SetBranchAddress("beam",&numbeam);
      tout->SetBranchAddress("particle",&particle);
      tout->SetBranchAddress("id",&id);
      tout->SetBranchAddress("Energy",&En);
      tout->SetBranchAddress("x",&xx);
      tout->SetBranchAddress("y",&yy);
      tout->SetBranchAddress("z",&zz);
      tout->SetBranchAddress("cx",&cxx);
      tout->SetBranchAddress("cy",&cyy);
      tout->SetBranchAddress("cz",&czz);
      tout->SetBranchAddress("ctime",&ctime);
    }

    int numtracks = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    for (int i = 0; i < numtracks; i++)
    {
        Track *track = RunManager::getInstance()->GetTrackPropagation()->GetTrack(i);

        numbeam = RunManager::getInstance()->GetGlobalIter();
        particle = (track->GetParticleName() == "h") ? 0 : ((track->GetParticleName() == "e-") ? 1 : ((track->GetParticleName() == "ex") ? 2 : ((track->GetParticleName() == "photon") ? 3 : 4)));
        id = track->GetTrackID();
        En = track->GetPostStepPoint().energy / unit_eV; 
        xx = (track->GetPostStepPoint().point.X() - track->GetInitialPoint().X()) / unit_nm;
        yy = (track->GetPostStepPoint().point.Y() - track->GetInitialPoint().Y()) / unit_nm;
        zz = (track->GetPostStepPoint().point.Z() - track->GetInitialPoint().Z()) / unit_nm;
        cxx = track->GetInitialPoint().X() / unit_nm;
        cyy = track->GetInitialPoint().Y() / unit_nm;
        czz = track->GetInitialPoint().Z() / unit_nm; 
        ctime = track->GetCreationTime() / unit_fs;       

        tout->Fill();
    }

    gDirectory->cd(topdir);
    tout->Write(0,TObject::kWriteDelete,0);
    tout = nullptr;

    fout->Close();
    fout = nullptr;
}
