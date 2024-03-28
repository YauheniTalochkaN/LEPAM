#include "Track.hh"

Track::Track()
{

}

Track::~Track()
{

}

TString Track::GetParticleName()
{
    return particle_name;
}

TString Track::GetParticleType()
{
    return particle_type;
}

Int_t Track::GetTrackID()
{
    return track_id;
}

TVector3 Track::GetInitialPoint()
{
    return initial_point;
}

Double_t Track::GetCreationTime()
{
    return creation_time;
}

TString Track::GetCreatorProcess()
{
    return creator_process;
}

Double_t Track::GetMass()
{
    return mass;
}

Double_t Track::GetCharge()
{
    return charge;
}

TrackPoint Track::GetPreStepPoint()
{
    return pre_step_point;
}

TrackPoint Track::GetPostStepPoint()
{
    return post_step_point;
}

bool Track::GetKilled()
{
    return killed;
}
