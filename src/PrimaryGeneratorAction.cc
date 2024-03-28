#include "PrimaryGeneratorAction.hh"
#include "RunManager.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
{

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  
}

void PrimaryGeneratorAction::SetParticleEnergy(TString name, Double_t energy)
{
    int track_num = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    if(track_num == 0) return;

    for (int i = 0; i < track_num; i++)
    {
        Track* track = RunManager::getInstance()->GetTrackPropagation()->GetTrack(i);
        if(track->particle_name == name)
        {
            track->pre_step_point.kinetic_energy = energy; 
            track->post_step_point.kinetic_energy = energy;

            track->pre_step_point.energy = energy;
            track->post_step_point.energy = energy;
            
            if((track->particle_name == "e-") || (track->particle_name == "ex"))
            {
                TString volume_name = RunManager::getInstance()->GetPhysicsList()->GetTGeoVolumeAtPoint(track->post_step_point.point);
                if(volume_name != "")
                {
                    MaterialProperties *mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
                    if(mp != nullptr)
                    {
                        Double_t Eg = mp->GetBandGap();
                        if(Eg >= 0) {track->pre_step_point.energy += Eg; track->post_step_point.energy += Eg;}
                    }
                }
            }

            if(track->particle_name == "h")
            {
                track->pre_step_point.kinetic_energy *= -1.0;
                track->post_step_point.kinetic_energy *= -1.0;

                TString volume_name = RunManager::getInstance()->GetPhysicsList()->GetTGeoVolumeAtPoint(track->post_step_point.point);
                if(volume_name != "")
                {
                    MaterialProperties *mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
                    if(mp != nullptr)
                    {
                        LinearInterpolation *dos = mp->GetDensityOfStates();
                        if(dos != nullptr) 
                        {
                            if(energy < dos->GetXmin()) 
                            {
                                track->pre_step_point.kinetic_energy = 0; 
                                track->post_step_point.kinetic_energy = 0;
                            }
                        }
                    }
                }
            }

            RunManager::getInstance()->GetPhysicsList()->UpdateMass(track);
            RunManager::getInstance()->GetPhysicsList()->UpdateVelocity(track);
        }
    }
}

void PrimaryGeneratorAction::GeneratePrimaries()
{    
    TrackPoint primary_particle;
    primary_particle.kinetic_energy = 1.0 * unit_eV;
    primary_particle.momentum_direction = TVector3(1, 0, 0); // RunManager::getInstance()->GetPhysicsList()->RandomTVector3();
    primary_particle.point = TVector3(0, 0, 0);
    primary_particle.time = 0;
    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(RunManager::getInstance()->GetPhysicsList()->CreateElectron(primary_particle));

/*  TrackPoint primary_particle;
    primary_particle.energy = -1.0 * unit_eV;
    primary_particle.momentum_direction = TVector3(1, 0, 0); // RunManager::getInstance()->GetPhysicsList()->RandomTVector3();
    primary_particle.point = TVector3(0, 0, 0);
    primary_particle.time = 0;
    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(RunManager::getInstance()->GetPhysicsList()->CreateHole(primary_particle));*/

/*  TrackPoint primary_particle;
    primary_particle.kinetic_energy = 100*unit_meV;
    primary_particle.momentum_direction = TVector3(1, 0, 0); // RunManager::getInstance()->GetPhysicsList()->RandomTVector3();
    primary_particle.point = TVector3(0,0,0);
    primary_particle.time = 0;
    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(RunManager::getInstance()->GetPhysicsList()->CreateExciton(primary_particle));*/
/*
    for (int i = 0; i < 1; i++)
    {
        Double_t ksi = gRandom->Gaus(0.2, 0.01) * unit_eV;
        while ((ksi < 0.1 * unit_eV) || (ksi > 0.7 * unit_eV))
        {
            ksi = gRandom->Gaus(0.2, 0.01) * unit_eV;
        }
    
        TVector3 position = TVector3(-48.0 * unit_nm + 96.0 * unit_nm * (1.0 - gRandom->Rndm()), -48.0 * unit_nm + 96.0 * unit_nm * (1.0 - gRandom->Rndm()), -48.0 * unit_nm + 96.0 * unit_nm * (1.0 - gRandom->Rndm()));
            
        TrackPoint primary_particle;
        primary_particle.kinetic_energy = ksi;
        primary_particle.momentum_direction = RunManager::getInstance()->GetPhysicsList()->RandomTVector3();
        primary_particle.point = position;
        primary_particle.time = 0;
        RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(RunManager::getInstance()->GetPhysicsList()->CreateElectron(primary_particle));

        ksi = gRandom->Gaus(0.2, 0.01) * unit_eV;
        while ((ksi < 0.1 * unit_eV) || (ksi > 0.7 * unit_eV))
        {
            ksi = gRandom->Gaus(0.2, 0.01) * unit_eV;
        }

        TrackPoint primary_particle2;
        primary_particle2.energy = -ksi;
        primary_particle2.momentum_direction = RunManager::getInstance()->GetPhysicsList()->RandomTVector3();
        primary_particle2.point = position + 4.0 * unit_nm * RunManager::getInstance()->GetPhysicsList()->RandomTVector3();
        primary_particle2.time = 0;
        RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(RunManager::getInstance()->GetPhysicsList()->CreateHole(primary_particle2));
    }*/
}
