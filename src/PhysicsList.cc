#include "PhysicsList.hh"
#include "RunManager.hh"

PhysicsList::PhysicsList()
{
    TIME_CUT_eh = 1.0 * unit_fs;
    EH_GLOBAL_TIME = 0;
    Electric_Field_Strength = TVector3(0,0,0);
    Magnetic_Field_Induction = TVector3(0,0,0);
}

PhysicsList::~PhysicsList() 
{
    if(PeriodicCube != nullptr) delete PeriodicCube;
}

void PhysicsList::UpdateTimeCut()
{

}

void PhysicsList::SetEHGlobalTime(Double_t val)
{
    EH_GLOBAL_TIME = val;
}

Double_t PhysicsList::GetEHGlobalTime()
{
    return EH_GLOBAL_TIME;
}

void PhysicsList::SetElectricFieldStrength(TVector3 Evec)
{
    Electric_Field_Strength = Evec;
}

TVector3 PhysicsList::GetElectricFieldStrength()
{
    return Electric_Field_Strength;
}

void PhysicsList::SetMagneticFieldInduction(TVector3 Bvec)
{
    Magnetic_Field_Induction = Bvec;
}

TVector3 PhysicsList::GetMagneticFieldInduction()
{
    return Magnetic_Field_Induction;
}

void PhysicsList::SetPeriodicBoundaries(TVector3 v1, TVector3 v2)
{
    PeriodicCube = new PeriodicBoundaries(v1, v2);
}

PeriodicBoundaries* PhysicsList::GetPeriodicBoundaries()
{
    return PeriodicCube;
}

bool PhysicsList::ApplyPeriodicity(TVector3& vec)
{
    if(PeriodicCube != nullptr)
    {
        Int_t n1 = (Int_t)(2.0 * (vec-PeriodicCube->center).X() / PeriodicCube->edges.X()),
              n2 = (Int_t)(2.0 * (vec-PeriodicCube->center).Y() / PeriodicCube->edges.Y()),
              n3 = (Int_t)(2.0 * (vec-PeriodicCube->center).Z() / PeriodicCube->edges.Z());
    
        Double_t xx = vec.X() - n1 * PeriodicCube->edges.X() / 2.0, yy = vec.Y() - n2 * PeriodicCube->edges.Y() / 2.0, zz = vec.Z() - n3 * PeriodicCube->edges.Z() / 2.0;

        if((abs(n1) > 0) && (abs(n1) % 2 > 0)) {if(xx >= 0) {xx -= 0.5 * PeriodicCube->edges.X();} else xx += 0.5 * PeriodicCube->edges.X();}
        if((abs(n2) > 0) && (abs(n2) % 2 > 0)) {if(yy >= 0) {yy -= 0.5 * PeriodicCube->edges.Y();} else yy += 0.5 * PeriodicCube->edges.Y();}
        if((abs(n3) > 0) && (abs(n3) % 2 > 0)) {if(zz >= 0) {zz -= 0.5 * PeriodicCube->edges.Z();} else zz += 0.5 * PeriodicCube->edges.Z();}

        vec = TVector3(xx, yy, zz);

        if((abs(n1) > 0) || (abs(n2) > 0) || (abs(n3) > 0)) return true;
        else return false;
    }
    else return false;
}

TVector3 PhysicsList::RandomTVector3()
{
    Double_t rx, ry, rz;
    gRandom->Sphere(rx, ry, rz, 1);
    return TVector3(rx, ry, rz);
}

Double_t PhysicsList::RandomValue(LinearInterpolation* func, Double_t xmin, Double_t xmax)
{
    if(xmin == xmax) return xmin;
    
    if(xmin < func->GetXmin()) xmin = func->GetXmin();
    if(xmax > func->GetXmax()) xmax = func->GetXmax();

    Double_t ymin = func->GetYmin(), 
             ymax = func->GetYmax();

    if((xmin != func->GetXmin()) || (xmax != func->GetXmax()))
    {
    Double_t *xx_mass = func->GetXMass();
    Double_t *yy_mass = func->GetYMass();
    Int_t nbins = func->GetNBins();

    ymax = func->GetYmin();
    ymin = func->GetYmax();

    for (int i = 0; i < nbins; i++)
    {
        if((xx_mass[i] >= xmin) && (xx_mass[i] <= xmax))
        {
            if(yy_mass[i] > ymax) ymax = yy_mass[i];
            if(yy_mass[i] < ymin) ymin = yy_mass[i];
        }
    }
    }

    if(ymin == ymax) return xmin + gRandom->Rndm() * (xmax - xmin);

    while(true)
    {
        Double_t x = xmin + gRandom->Rndm() * (xmax - xmin);
        Double_t F = (func->GetValue(x) - ymin) / (ymax - ymin);
        Double_t q = gRandom->Rndm();
	  
        if(q <= F){return x;}
    }
}

Int_t PhysicsList::DiscreteRandomValue(std::vector<Double_t> prob)
{
    if(prob.size() < 2) return 0;
    
    Double_t Stat_Sum = 0.0;
    std::vector<Double_t> loc_vec;
    for (auto it = prob.begin() ; it != prob.end(); ++it)
    {
      Stat_Sum += (*it);
      loc_vec.push_back(Stat_Sum);
    }
    
    Int_t index = 0;
    Double_t ksi = gRandom->Rndm();
    for (auto it = loc_vec.begin() ; it != loc_vec.end(); ++it)
    { 
      if(ksi < ((*it)/Stat_Sum)) {index = std::distance(loc_vec.begin(), it); break;}
    }
    loc_vec.clear();

    return index;
}

TString PhysicsList::GetTGeoVolumeAtPoint(TVector3 point)
{
    RunManager::getInstance()->GetTrackPropagation()->LockMutex();
    
    TGeoNode* node = RunManager::getInstance()->GetTGeoManager()->FindNode(geo_size_factor * point.X(), geo_size_factor * point.Y(), geo_size_factor * point.Z());
    if(node == nullptr) {RunManager::getInstance()->GetTrackPropagation()->UnLockMutex(); return "";}

    TString name = node->GetVolume()->GetName();
    node = nullptr;

    RunManager::getInstance()->GetTrackPropagation()->UnLockMutex();

    return name;
}

void PhysicsList::CutTrackLine(Track* track)
{
    if(track->visible)
    {
        TGeoManager* geom = RunManager::getInstance()->GetTGeoManager();
    
        TVirtualGeoTrack* geo_track;

        RunManager::getInstance()->GetTrackPropagation()->LockMutex();
        int track_step = geom->AddTrack(geom->GetNtracks(), 1);
        RunManager::getInstance()->GetTrackPropagation()->UnLockMutex();

        geo_track = geom->GetTrack(track_step);
        geo_track->SetLineColor(track->geo_track->GetLineColor());
        geo_track->AddPoint(geo_size_factor * track->post_step_point.point.X(), geo_size_factor * track->post_step_point.point.Y(), geo_size_factor * track->post_step_point.point.Z(), geo_time_factor * track->post_step_point.time);
        
        track->geo_track = geo_track;

        geom = nullptr;
    }
}

void PhysicsList::SetTGeoTrack(Track* track)
{
    if(track->visible)
    {
        TGeoManager* geom = RunManager::getInstance()->GetTGeoManager();
    
        TVirtualGeoTrack* geo_track;
        int track_step = geom->AddTrack(geom->GetNtracks(), 1);
        geo_track = geom->GetTrack(track_step);

        if(track->particle_name == "e-") geo_track->SetLineColor(kRed);
        else if(track->particle_name == "h") geo_track->SetLineColor(kBlue);
        else if(track->particle_name == "ex") geo_track->SetLineColor(kMagenta);
        else if(track->particle_name == "photon") geo_track->SetLineColor(kGreen);
        else geo_track->SetLineColor(kBlack);

        geo_track->AddPoint(geo_size_factor * track->post_step_point.point.X(), geo_size_factor * track->post_step_point.point.Y(), geo_size_factor * track->post_step_point.point.Z(), geo_time_factor * track->post_step_point.time);

        track->geo_track = geo_track;

        geom = nullptr;
    }
}

Track* PhysicsList::CreateElectron(TrackPoint trp, TString type, TString volume_name, MaterialProperties* mp)
{
    Track* tr = new Track;

    tr->particle_name = "e-";
    tr->particle_type = type;
    tr->track_id = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    tr->mass = const_me;
    tr->charge = -const_e;
    tr->pre_step_point = trp;
    tr->initial_point = trp.point;
    tr->creation_time = trp.time;
    tr->creator_process = "";
    tr->killed = false;
    tr->should_be_killed = false;
    tr->visible = RunManager::getInstance()->GetVisible();
    tr->pre_step_point.velocity = CalculateVelocity(tr->mass, tr->pre_step_point.kinetic_energy);
    tr->pre_step_point.time_step = TIME_CUT_eh;
    
    Double_t Eg = 0;
    if(volume_name == "") volume_name = GetTGeoVolumeAtPoint(trp.point);
    if(volume_name != "") 
    {
        tr->pre_step_point.volume = volume_name;
        
        if(mp == nullptr) mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
        if(mp != nullptr)
        {
            tr->pre_step_point.matprop = mp;

            Eg = mp->GetBandGap();
        }
    }
    tr->pre_step_point.energy =  Eg + tr->pre_step_point.kinetic_energy;
    
    tr->post_step_point = tr->pre_step_point;

    UpdateMass(tr);
    UpdateVelocity(tr);

    return tr;
}

Track* PhysicsList::CreateHole(TrackPoint trp, TString type, TString volume_name, MaterialProperties* mp)
{
    Track* tr = new Track;

    tr->particle_name = "h";
    tr->particle_type = type;
    tr->track_id = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    tr->mass = const_me;
    tr->charge = const_e;
    tr->pre_step_point = trp;
    tr->initial_point = trp.point;
    tr->creation_time = trp.time;
    tr->creator_process = "";
    tr->killed = false;
    tr->should_be_killed = false;
    tr->visible = RunManager::getInstance()->GetVisible();
    tr->pre_step_point.velocity = CalculateVelocity(tr->mass, tr->pre_step_point.kinetic_energy);
    tr->pre_step_point.time_step = TIME_CUT_eh;

    tr->pre_step_point.kinetic_energy = -tr->pre_step_point.energy;

    if(volume_name == "") volume_name = GetTGeoVolumeAtPoint(trp.point);
    if(volume_name != "")
    {
        tr->pre_step_point.volume = volume_name;
        
        if(mp == nullptr) mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
        if(mp != nullptr)
        {
            tr->pre_step_point.matprop = mp;
            
            LinearInterpolation *dos = mp->GetDensityOfStates();
            if(dos != nullptr) 
            {
                if(tr->pre_step_point.energy < dos->GetXmin()) 
                {
                    tr->pre_step_point.kinetic_energy = 0; 
                }
            }
        }
    }
    
    tr->post_step_point = tr->pre_step_point;
    
    UpdateMass(tr);
    UpdateVelocity(tr);

    return tr;
}

Track* PhysicsList::CreateExciton(TrackPoint trp, TString type, TString volume_name, MaterialProperties* mp, Double_t mex)
{
    Track* tr = new Track;

    tr->particle_name = "ex";
    tr->particle_type = type;
    tr->track_id = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    tr->mass = mex;
    tr->charge = 0;
    tr->pre_step_point = trp;
    tr->initial_point = trp.point;
    tr->creation_time = trp.time;
    tr->creator_process = "";
    tr->killed = false;
    tr->should_be_killed = false;
    tr->visible = RunManager::getInstance()->GetVisible();
    tr->pre_step_point.velocity = CalculateVelocity(tr->mass, tr->pre_step_point.kinetic_energy);
    
    Double_t Eg = 0;
    if(volume_name == "") volume_name = GetTGeoVolumeAtPoint(trp.point);
    if(volume_name != "") 
    {
        tr->pre_step_point.volume = volume_name;
        
        if(mp == nullptr) mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
        if(mp != nullptr)
        {
            tr->pre_step_point.matprop = mp;

            Eg = mp->GetBandGap();
        }
    }
    tr->pre_step_point.energy =  Eg + tr->pre_step_point.kinetic_energy;
    
    tr->post_step_point = tr->pre_step_point;

    UpdateMass(tr);
    UpdateVelocity(tr);

    return tr;
}

Track* PhysicsList::CreatePhoton(TrackPoint trp, TString volume_name, MaterialProperties* mp)
{
    Track* tr = new Track;

    tr->particle_name = "photon";
    tr->particle_type = "";
    tr->track_id = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    tr->mass = 0;
    tr->charge = 0;
    tr->pre_step_point = trp;
    tr->initial_point = trp.point;
    tr->creation_time = trp.time;
    tr->creator_process = "";
    tr->killed = false;
    tr->should_be_killed = false;
    tr->visible = RunManager::getInstance()->GetVisible();
    tr->pre_step_point.velocity = CalculateVelocity(tr->mass, tr->pre_step_point.energy);
    tr->pre_step_point.kinetic_energy = tr->pre_step_point.energy;

    if(volume_name == "") volume_name = GetTGeoVolumeAtPoint(trp.point);
    if(volume_name != "") 
    {
        tr->pre_step_point.volume = volume_name;
        
        if(mp == nullptr) mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
        if(mp != nullptr) tr->pre_step_point.matprop = mp;
    }
    
    tr->post_step_point = tr->pre_step_point;

    UpdateVelocity(tr);

    return tr;
}

Track* PhysicsList::CreateCenter(TrackPoint trp, TString name, TString type, TString volume_name, MaterialProperties* mp)
{
    Track* tr = new Track;

    tr->particle_name = name;
    tr->particle_type = type;
    tr->track_id = RunManager::getInstance()->GetTrackPropagation()->SizeofTrackList();
    tr->mass = 0;
    tr->charge = 0;
    tr->pre_step_point = trp;
    tr->initial_point = trp.point;
    tr->creation_time = trp.time;
    tr->creator_process = "";
    tr->killed = false;
    tr->should_be_killed = false;
    tr->visible = RunManager::getInstance()->GetVisible();
    tr->pre_step_point.velocity = 0;
    tr->pre_step_point.kinetic_energy = 0;

    if(volume_name == "") volume_name = GetTGeoVolumeAtPoint(trp.point);
    if(volume_name != "") 
    {
        tr->pre_step_point.volume = volume_name;
        
        if(mp == nullptr) mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
        if(mp != nullptr) tr->pre_step_point.matprop = mp;
    }
    
    tr->post_step_point = tr->pre_step_point;

    return tr;
}

void PhysicsList::UpdateVelocity(Track* track)
{    
    if(track->post_step_point.matprop != nullptr)
    {
        LinearInterpolation *gv = track->post_step_point.matprop->GetGroupVelocity(track->particle_name);
        if(gv != nullptr)
        {
            if(track->particle_name == "h")
            {
                if((track->post_step_point.energy >= gv->GetXmin()) && (track->post_step_point.energy <= gv->GetXmax())) 
                {
                    track->post_step_point.velocity = gv->GetValue(track->post_step_point.energy);
                    return;
                }
            }
            else
            {
                if((track->post_step_point.kinetic_energy >= gv->GetXmin()) && (track->post_step_point.kinetic_energy <= gv->GetXmax())) 
                {
                    track->post_step_point.velocity = gv->GetValue(track->post_step_point.kinetic_energy);
                    return;
                }
            }
        }
    }

    track->post_step_point.velocity = CalculateVelocity(track->mass,  track->post_step_point.kinetic_energy);
}

void PhysicsList::UpdateMass(Track* track)
{
    if(track->post_step_point.matprop == nullptr) return;

    LinearInterpolation *effmass = track->post_step_point.matprop->GetEffectiveMass(track->particle_name);
    if(effmass == nullptr) return;

    if(track->particle_name == "h")
    {
        if((track->post_step_point.energy >= effmass->GetXmin()) && (track->post_step_point.energy <= effmass->GetXmax())) 
        track->mass = effmass->GetValue(track->post_step_point.energy);
    }
    else
    {
        if((track->post_step_point.kinetic_energy >= effmass->GetXmin()) && (track->post_step_point.kinetic_energy <= effmass->GetXmax())) 
        track->mass = effmass->GetValue(track->post_step_point.kinetic_energy); 
    }  
}

Double_t PhysicsList::CalculateMomentum(Double_t mass, Double_t kin_energy)
{
    return (kin_energy <= 0) ? 0 : sqrt(kin_energy * (kin_energy + 2.0 * mass * const_c * const_c)) / const_c;
}

Double_t PhysicsList::CalculateVelocity(Double_t mass, Double_t kin_energy)
{
    return (kin_energy <= 0) ? 0 : const_c * sqrt(kin_energy * (kin_energy + 2.0 * mass * pow(const_c, 2.0))) / (kin_energy + mass * pow(const_c, 2.0));
}

Double_t PhysicsList::GetDistanceBetweenTwoTracks(TVector3 r1, TVector3 r2, TVector3 r3, TVector3 r4)
{
    TVector3 dn = r4 - r2 - r3 + r1, dr0 = r2-r1;

    Double_t start_deriv = (dr0 + 0.0001*dn).Mag() - (dr0).Mag(),
             stop_deriv = (dr0 + dn).Mag() - (dr0 + 0.9999*dn).Mag();

    if(start_deriv * stop_deriv < 0)
    {
        Double_t dt = 0.01, t = 0;
        Double_t dD = -1;
        
        while ((dD < 0) && (t < 1))
        {
            t += dt;
            if(t > 1) t = 1;

            dD = (dr0 + t*dn).Mag() - (dr0 + (t-dt)*dn).Mag();
        }

        if(t < 1) t -= dt;

        return (dr0 + t*dn).Mag();
    }
    else return ((dr0 + dn).Mag() < (dr0).Mag()) ? (dr0 + dn).Mag() : (dr0).Mag();
}

void PhysicsList::ElectronHoleInteraction(std::vector<Track*> track_list)
{      
    size_t track_list_size = (size_t)track_list.size();

    TVector3 **dp = new TVector3*[track_list_size];
    Double_t **Uc = new Double_t*[track_list_size], 
             **Umax = new Double_t*[track_list_size],
              *Eps = new Double_t[track_list_size],
             **radprob = new Double_t*[track_list_size];

    for (size_t i1 = 0; i1 < track_list_size; ++i1)
    {
        dp[i1] = new TVector3[track_list_size];
        Uc[i1] = new Double_t[track_list_size];
        Umax[i1] = new Double_t[track_list_size];
        radprob[i1] = new Double_t[track_list_size];
    }

    #pragma omp parallel for
    for (size_t i1 = 0; i1 < track_list_size; ++i1)
    {
        LinearInterpolation *projection = nullptr;

        if(track_list[i1]->post_step_point.matprop->GetElectronHoleRadiationRecombinationTime() != nullptr) projection = track_list[i1]->post_step_point.matprop->GetElectronHoleRadiationRecombinationTime()->GetXProjection(track_list[i1]->post_step_point.kinetic_energy);
        
        radprob[i1][i1] = 0;
        Uc[i1][i1] = 0;
        Umax[i1][i1] = 0;
        dp[i1][i1] = TVector3(0,0,0);
        
        for (size_t i2 = i1+1; i2 < track_list_size; ++i2)
        {
            dp[i1][i2] = TVector3(0,0,0);
            Uc[i1][i2] = 0;
            Umax[i1][i2] = 0;
            radprob[i1][i2] = 0;
            dp[i2][i1] = TVector3(0,0,0);
            Uc[i2][i1] = 0;
            Umax[i2][i1] = 0;
            radprob[i2][i1] = 0;

            if((projection != nullptr) && (((track_list[i1]->particle_name == "e-") && (track_list[i2]->particle_name == "h")) || ((track_list[i1]->particle_name == "h") && (track_list[i2]->particle_name == "e-"))))
            {
                radprob[i1][i2] = 1.0 / projection->GetValue(track_list[i2]->post_step_point.kinetic_energy);
                radprob[i2][i1] = radprob[i1][i2];
            }
        }

        Eps[i1] = track_list[i1]->post_step_point.matprop->GetDielectricPermittivity();

        delete projection;
    }

    for (size_t i1 = 0; i1 < track_list_size; ++i1)
    {
        if(track_list[i1]->should_be_killed == false)
        {
            Double_t loc_radprob_tot = 0;
            std::vector<Double_t> loc_radprob_list;
            for (size_t i2 = 0; i2 < track_list_size; ++i2)
            {
                if(track_list[i2]->should_be_killed == false)
                {
                    loc_radprob_tot += radprob[i1][i2];
                    loc_radprob_list.push_back(radprob[i1][i2]);
                }
                else loc_radprob_list.push_back(0);
            }

            if(gRandom->Rndm() < 1.0 - TMath::Exp(-track_list[i1]->post_step_point.time_step * loc_radprob_tot))
            {            
                Int_t i2 = DiscreteRandomValue(loc_radprob_list);

                TrackPoint secondary_particle;
                secondary_particle.energy = abs(track_list[i1]->post_step_point.energy - track_list[i2]->post_step_point.energy);
                secondary_particle.momentum_direction = RandomTVector3();
                secondary_particle.point = (track_list[i1]->mass * track_list[i1]->post_step_point.point + track_list[i2]->mass * track_list[i2]->post_step_point.point) * pow(track_list[i1]->mass + track_list[i2]->mass, -1.0);
                secondary_particle.time = track_list[i1]->post_step_point.time + track_list[i1]->post_step_point.time_step;
                Track* ph_track = CreatePhoton(secondary_particle, track_list[i1]->post_step_point.volume, track_list[i1]->post_step_point.matprop);
                ph_track->creator_process = "Electron-hole radiation recombination";

                RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(ph_track);

                track_list[i1]->should_be_killed = true;
                track_list[i2]->should_be_killed = true;
            }

            loc_radprob_list.clear();
        }

        delete [] radprob[i1];
    }
    delete [] radprob;

    #pragma omp parallel for
    for (size_t i1 = 0; i1 < track_list_size; ++i1)
    {
        TVector3 ri1 = track_list[i1]->post_step_point.point, rf1 = track_list[i1]->post_step_point.point + track_list[i1]->post_step_point.momentum_direction * track_list[i1]->post_step_point.time_step * track_list[i1]->post_step_point.velocity;

        for (size_t i2 = i1+1; i2 < track_list_size; ++i2)
        {       
            Double_t eps_loc = pow(1.0/Eps[i1] + 1.0/Eps[i2], -1.0);

            TVector3 ri2 = track_list[i2]->post_step_point.point, rf2 = track_list[i2]->post_step_point.point + track_list[i2]->post_step_point.momentum_direction * track_list[i2]->post_step_point.time_step * track_list[i2]->post_step_point.velocity;
                   
            Double_t rmin = GetDistanceBetweenTwoTracks(ri1, ri2, rf1, rf2);

            Umax[i1][i2] = const_k / eps_loc * track_list[i1]->charge * track_list[i2]->charge / sqrt(rmin * rmin + unit_nm * unit_nm);
            Umax[i2][i1] = Umax[i1][i2];
            
            if(PeriodicCube == nullptr)
            {
                TVector3 F_loc = const_k / eps_loc * track_list[i1]->charge * track_list[i2]->charge * pow((ri1 - ri2).Mag2() + unit_nm * unit_nm, -1.0) * pow((ri1 - ri2).Mag() + 1E-10 * unit_nm, -1.0) * (ri1 - ri2);
                dp[i1][i2] = F_loc * track_list[i1]->post_step_point.time_step;
                dp[i2][i1] = -F_loc * track_list[i1]->post_step_point.time_step;

                Double_t efi = const_k / eps_loc * track_list[i1]->charge * track_list[i2]->charge * pow(sqrt((rf1 - rf2).Mag2() + unit_nm * unit_nm), -1.0);
                Uc[i1][i2] = efi;
                Uc[i2][i1] = efi;
            }
            else 
            {
                for (size_t j1 = -1; j1 < 2; j1++)
                {
                    for (size_t j2 = -1; j2 < 2; j2++)
                    {
                        for (size_t j3 = -1; j3 < 2; j3++)
                        {
                            TVector3 rr_f = rf1 - rf2 - TVector3(j1 * PeriodicCube->edges.X(), j2 * PeriodicCube->edges.Y(), j3 * PeriodicCube->edges.Z()),
                                     rr_i = ri1 - ri2 - TVector3(j1 * PeriodicCube->edges.X(), j2 * PeriodicCube->edges.Y(), j3 * PeriodicCube->edges.Z());

                            TVector3 F_loc = const_k / eps_loc * track_list[i1]->charge * track_list[i2]->charge * pow(rr_i.Mag2() + unit_nm * unit_nm, -1.0) * pow(rr_i.Mag() + 1E-10 * unit_nm, -1.0) * rr_i;
                            dp[i1][i2] = F_loc * track_list[i1]->post_step_point.time_step;
                            dp[i2][i1] = -F_loc * track_list[i1]->post_step_point.time_step;

                            Double_t efi = const_k / eps_loc * track_list[i1]->charge * track_list[i2]->charge * pow(sqrt(rr_f.Mag2() + unit_nm * unit_nm), -1.0);
                            Uc[i1][i2] = efi;
                            Uc[i2][i1] = efi;
                        }
                    }
                }
            }
        }
    }
    delete [] Eps;

    #pragma omp parallel for
    for (size_t i1 = 0; i1 < track_list_size; ++i1)
    {           
        TVector3 p = CalculateMomentum(track_list[i1]->mass, track_list[i1]->post_step_point.kinetic_energy) * track_list[i1]->post_step_point.momentum_direction;
        TVector3 loc_dp = TVector3(0,0,0);
        Double_t Uctot = 0;
        
        for (size_t i2 = 0; i2 < track_list_size; ++i2)
        {
            loc_dp += dp[i1][i2];
            Uctot += Uc[i1][i2];
        }
        
        loc_dp +=  track_list[i1]->charge * (Electric_Field_Strength + track_list[i1]->post_step_point.velocity * track_list[i1]->post_step_point.momentum_direction.Cross(Magnetic_Field_Induction)) * track_list[i1]->post_step_point.time_step;
        
        track_list[i1]->post_step_point.momentum_direction = pow((p+loc_dp).Mag(), -1.0) * (p+loc_dp);

        if(track_list[i1]->post_step_point.time == track_list[i1]->pre_step_point.time)
        {
            track_list[i1]->pre_step_point.Ecoulomb_energy = Uctot;
            track_list[i1]->post_step_point.Ecoulomb_energy = Uctot;
        }
        else track_list[i1]->post_step_point.Ecoulomb_energy = Uctot;
        
        if(track_list[i1]->post_step_point.kinetic_energy >= 0) 
        {
            Double_t dEk = sqrt((p+loc_dp).Mag2() + pow(track_list[i1]->mass * const_c, 2.0)) * const_c - track_list[i1]->mass * const_c * const_c - track_list[i1]->post_step_point.kinetic_energy;
            track_list[i1]->post_step_point.kinetic_energy += dEk;

            if(track_list[i1]->particle_name == "h") {track_list[i1]->post_step_point.energy -= dEk;}
            else {track_list[i1]->post_step_point.energy += dEk;}
        }
        else
        {
            if(track_list[i1]->particle_name == "h") {track_list[i1]->post_step_point.energy += track_list[i1]->post_step_point.Ecoulomb_energy - track_list[i1]->pre_step_point.Ecoulomb_energy;}
            else {track_list[i1]->post_step_point.energy -= track_list[i1]->post_step_point.Ecoulomb_energy - track_list[i1]->pre_step_point.Ecoulomb_energy;}
            
            track_list[i1]->post_step_point.kinetic_energy -= track_list[i1]->post_step_point.Ecoulomb_energy - track_list[i1]->pre_step_point.Ecoulomb_energy; 
        }

        UpdateMass(track_list[i1]);
        UpdateVelocity(track_list[i1]);

        delete [] dp[i1];
        delete [] Uc[i1];
    }
    delete [] dp;
    delete [] Uc;

    std::pair<Double_t,Int_t> *ex_list = new std::pair<Double_t,Int_t>[track_list_size];

    #pragma omp parallel for
    for (size_t i1 = 0; i1 < track_list_size; ++i1)
    {
        ex_list[i1].first = 10e10;
        ex_list[i1].second = 0;

        if(track_list[i1]->should_be_killed == false)
        {
            for (size_t i2 = 0; i2 < track_list_size; ++i2)
            {
                if((track_list[i2]->should_be_killed == false) && (((track_list[i1]->particle_name == "e-") && (track_list[i2]->particle_name == "h")) || ((track_list[i1]->particle_name == "h") && (track_list[i2]->particle_name == "e-"))))
                {
                    TVector3 p1 = CalculateMomentum(track_list[i1]->mass,  track_list[i1]->post_step_point.kinetic_energy) * track_list[i1]->post_step_point.momentum_direction,
                             p2 = CalculateMomentum(track_list[i2]->mass,  track_list[i2]->post_step_point.kinetic_energy) * track_list[i2]->post_step_point.momentum_direction;

                    Double_t Ktot = track_list[i1]->post_step_point.kinetic_energy + track_list[i2]->post_step_point.kinetic_energy,
                             Kc = const_c * sqrt((p1+p2).Mag2() + pow(track_list[i1]->mass + track_list[i2]->mass, 2.0) * const_c * const_c) - (track_list[i1]->mass + track_list[i2]->mass) * const_c * const_c;
                    
                    Double_t dE = Ktot + Umax[i1][i2] - Kc;

                    if(dE < ex_list[i1].first)
                    {
                        ex_list[i1].first = dE;
                        ex_list[i1].second = i2;
                    }
                } 
            }
        }

        delete [] Umax[i1];
    }
    delete [] Umax;

    for (size_t i3 = 0; i3 < track_list_size; ++i3)
    {
        if((track_list[i3]->should_be_killed == false) && (ex_list[i3].first < 0))
        {
            Int_t i2 = ex_list[i3].second;

            if((track_list[i2]->should_be_killed == false) && (ex_list[i2].first < 0))
            {
                Int_t i1 = ex_list[i2].second;
            
                TVector3 p1 = CalculateMomentum(track_list[i1]->mass,  track_list[i1]->post_step_point.kinetic_energy) * track_list[i1]->post_step_point.momentum_direction,
                         p2 = CalculateMomentum(track_list[i2]->mass,  track_list[i2]->post_step_point.kinetic_energy) * track_list[i2]->post_step_point.momentum_direction;

                Double_t Kc = const_c * sqrt((p1+p2).Mag2() + pow(track_list[i1]->mass + track_list[i2]->mass, 2.0) * const_c * const_c) - (track_list[i1]->mass + track_list[i2]->mass) * const_c * const_c;

                TrackPoint secondary_particle;
                secondary_particle.kinetic_energy = Kc;
                secondary_particle.momentum_direction = (p1 + p2) * pow((p1 + p2).Mag(), -1.0);
                secondary_particle.point = (track_list[i1]->mass * track_list[i1]->post_step_point.point + track_list[i2]->mass * track_list[i2]->post_step_point.point) * pow(track_list[i1]->mass + track_list[i2]->mass, -1.0);
                secondary_particle.time = track_list[i1]->post_step_point.time + track_list[i1]->post_step_point.time_step;

                Track* ex_track = CreateExciton(secondary_particle, "", track_list[i1]->post_step_point.volume, track_list[i1]->post_step_point.matprop, track_list[i1]->mass + track_list[i2]->mass);
                ex_track->creator_process = "Electron-hole coupling";

                RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(ex_track);

                ex_track = nullptr;

                track_list[i1]->should_be_killed = true;
                track_list[i2]->should_be_killed = true;
            }
        }
    }
    delete [] ex_list;
}

void PhysicsList::RadiationDecayOfExciton(Track* track)
{  
    TrackPoint secondary_particle;
    secondary_particle.energy = track->post_step_point.energy;
    secondary_particle.momentum_direction = RandomTVector3();
    secondary_particle.point = track->post_step_point.point + track->post_step_point.momentum_direction * track->post_step_point.time_step * track->post_step_point.velocity;
    secondary_particle.time = track->post_step_point.time + track->post_step_point.time_step;
    Track* ph_track = CreatePhoton(secondary_particle, track->post_step_point.volume, track->post_step_point.matprop);
    ph_track->creator_process = "Exciton radiation decay";

    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(ph_track);

    track->should_be_killed = true;
}

void PhysicsList::DecayOfCenter(Track* track, TString cname, TString ctype, Double_t cenergy, bool rad)
{
    if(rad)
    {
        TrackPoint secondary_particle;
        secondary_particle.energy = abs(track->post_step_point.energy - cenergy);
        secondary_particle.momentum_direction = RandomTVector3();
        secondary_particle.point = track->post_step_point.point;
        secondary_particle.time = track->post_step_point.time + track->post_step_point.time_step;
        Track* ph_track = CreatePhoton(secondary_particle, track->post_step_point.volume, track->post_step_point.matprop);
        ph_track->creator_process = ctype + cname + " radiation decay";

        RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(ph_track);   
    }
    
    if(ctype == "ground") track->should_be_killed = true;
    
    track->particle_type = ctype;
    track->post_step_point.energy = cenergy;
}

void PhysicsList::ExcitonDissociation(Track* track)
{  
    Double_t Eke = gRandom->Rndm() * (track->post_step_point.kinetic_energy + 90.0*unit_meV);
    
    Double_t birth_time = TIME_CUT_eh * TMath::Ceil((track->post_step_point.time + track->post_step_point.time_step) / TIME_CUT_eh);

    TrackPoint secondary_particle;
    secondary_particle.kinetic_energy = Eke;
    secondary_particle.momentum_direction = RandomTVector3();
    secondary_particle.point = track->post_step_point.point + track->post_step_point.momentum_direction * track->post_step_point.time_step * track->post_step_point.velocity;
    secondary_particle.time = birth_time;
    Track* e_track = CreateElectron(secondary_particle, "", track->post_step_point.volume, track->post_step_point.matprop);
    e_track->creator_process = "Exciton dissociation";

    TrackPoint secondary_particle2;
    secondary_particle2.energy = -(track->post_step_point.kinetic_energy + 90.0*unit_meV - Eke);
    secondary_particle2.momentum_direction = RandomTVector3();
    secondary_particle2.point = track->post_step_point.point + track->post_step_point.momentum_direction * track->post_step_point.time_step * track->post_step_point.velocity;
    secondary_particle2.time = birth_time;
    Track* h_track = CreateHole(secondary_particle2, "", track->post_step_point.volume, track->post_step_point.matprop);
    h_track->creator_process = "Exciton dissociation";
    h_track->track_id += 1;  

    TVector3 p1 = CalculateMomentum(e_track->mass,  e_track->post_step_point.kinetic_energy) * e_track->post_step_point.momentum_direction,
             p2 = CalculateMomentum(h_track->mass,  h_track->post_step_point.kinetic_energy) * h_track->post_step_point.momentum_direction;
    
    Double_t Ekc = const_c * sqrt((p1+p2).Mag2() + pow(e_track->mass + h_track->mass, 2.0) * const_c * const_c) - (e_track->mass + h_track->mass) * const_c * const_c;
    
    Double_t Ekeh = 0;
    if(e_track->post_step_point.kinetic_energy > 0) Ekeh += e_track->post_step_point.kinetic_energy;
    if(h_track->post_step_point.kinetic_energy > 0) Ekeh += h_track->post_step_point.kinetic_energy;

    Double_t disp_length = const_k / track->post_step_point.matprop->GetDielectricPermittivity() * const_e * const_e / (Ekeh - Ekc + unit_meV);
    if(disp_length < 2.0*unit_nm) disp_length = 2.0*unit_nm;
    TVector3 disp_vec = disp_length * RandomTVector3();

    e_track->pre_step_point.point += disp_vec;
    e_track->post_step_point.point += disp_vec;
    h_track->pre_step_point.point -= disp_vec;
    h_track->post_step_point.point -= disp_vec;

    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(e_track);
    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(h_track);

    track->should_be_killed = true;
}

void PhysicsList::GenerateVectorInElasticScattering(Double_t pe, Double_t hq, TVector3 InVec, TVector3& direct)
{
    Double_t cosa = (2.0 * pow(pe,2.0) - pow(hq,2.0)) / (2.0 * pow(pe,2.0));

    direct = TVector3(sqrt(1.0 - cosa*cosa), 0.0, cosa);

    Double_t fi = 2.0 * M_PI * gRandom->Rndm();

    TMatrixD mzc = TMatrixD(3,3);
    mzc[0] = {std::cos(fi), -std::sin(fi), 0.0};
    mzc[1] = {std::sin(fi),  std::cos(fi), 0.0};
    mzc[2] = {0.0         ,  0.0         , 1.0};

    Double_t sqz = sqrt(1.0 - InVec.Z() * InVec.Z());
    TMatrixD mrc = TMatrixD(3,3);
    mrc[0] = {InVec.X() * InVec.Z() / (sqz + 1E-8), -InVec.Y() / (sqz + 1E-8), InVec.X()};
    mrc[1] = {InVec.Y() * InVec.Z() / (sqz + 1E-8),  InVec.X() / (sqz + 1E-8), InVec.Y()};
    mrc[2] = {-sqz                                ,  0.0                     , InVec.Z()};

    direct = mrc * (mzc * direct);

    direct = direct * pow(direct.Mag(), -1.0);
}

void PhysicsList::ElasticScattering(Track* track, Int_t index)
{   
  //  Container_emfp* cont_emfp = track->post_step_point.matprop->GetListOfElasticAverageInteractionTime(track->particle_name)[index];
     
  //  Double_t pe = sqrt(2.0 * track->mass * track->post_step_point.kinetic_energy);

  //  TVector3 direct;

  //  Double_t hrq= cont_emfp->mls->GetXProjection(track->post_step_point.kinetic_energy)->GetValue(gRandom->Rndm());

  //  GenerateVectorInElasticScattering(pe, hrq, track->post_step_point.momentum_direction, direct);

    track->post_step_point.momentum_direction = RandomTVector3();
}

void PhysicsList::GenerateVectorInPhononScattering(Double_t sqrp0, Double_t sqrp1, Double_t sqrp2, TVector3 InVec, TVector3& direct1)
{
    Double_t a = (sqrp0 + sqrp1 - sqrp2) / (2.0 * sqrt(sqrp0 * sqrp1) + 1E-80);
    Double_t b = (sqrp0 - sqrp1 + sqrp2) / (2.0 * sqrt(sqrp0 * sqrp2) + 1E-80);

    if(a > 1) a = 1;
    if(a < -1) a = -1;
    if(b > 1) b = 1;
    if(b < -1) b = -1;

    direct1 = TVector3(sqrt(1.0 - a*a), 0.0, a);

    Double_t fi = 2.0 * M_PI * gRandom->Rndm();

    TMatrixD mzc = TMatrixD(3,3);
    mzc[0] = {std::cos(fi), -std::sin(fi), 0.0};
    mzc[1] = {std::sin(fi),  std::cos(fi), 0.0};
    mzc[2] = {0.0         ,  0.0         , 1.0};

    Double_t sqz = sqrt(1.0 - InVec.Z() * InVec.Z());
    TMatrixD mrc = TMatrixD(3,3);
    mrc[0] = {InVec.X() * InVec.Z() / (sqz + 1E-8), -InVec.Y() / (sqz + 1E-8), InVec.X()};
    mrc[1] = {InVec.Y() * InVec.Z() / (sqz + 1E-8),  InVec.X() / (sqz + 1E-8), InVec.Y()};
    mrc[2] = {-sqz                                ,  0.0                     , InVec.Z()};

    direct1 = mrc * (mzc * direct1);

    direct1 = direct1 * pow(direct1.Mag(), -1.0);
}

void PhysicsList::ThermalScattering(Track* track, Int_t index)
{    
    LinearInterpolation* dos = track->post_step_point.matprop->GetDensityOfStates();

    Container_ait* cont_ait = track->post_step_point.matprop->GetListOfThermalAverageInteractionTime(track->particle_name)[index];

    Double_t hw = 0, hrq = 0, post_energy = 0;

    LinearInterpolation *projection = cont_ait->mls->GetXProjection(track->post_step_point.kinetic_energy);

    Int_t while_iter = 0;
    while(true)
    {
        if(while_iter > 1000) 
        {
            std::cout << "PhysicsList: Warning - To much iterations for an thermal scattering event!\n"; 
            
            delete projection; 
            dos = nullptr;
            cont_ait = nullptr;
            
            return;
        }
    
        hrq = RandomValue(projection, projection->GetXmin(), projection->GetXmax());
        hw = cont_ait->displow->GetValue(hrq);
        post_energy = (track->particle_name == "h") ? (track->post_step_point.energy + hw) : (track->post_step_point.energy - hw);

        if(track->particle_name == "ex") 
        {
            if(track->post_step_point.kinetic_energy - hw > 0) break; 
        }
        else if(dos->GetValue(post_energy) > 0) break;

        while_iter++;
    }

    delete projection;

    Double_t sqrp0 = pow(CalculateMomentum(track->mass, track->post_step_point.kinetic_energy), 2.0);

    track->post_step_point.kinetic_energy -= hw;
    track->post_step_point.energy = post_energy;
    UpdateMass(track);
    UpdateVelocity(track);

    Double_t sqrp1 = pow(CalculateMomentum(track->mass, track->post_step_point.kinetic_energy), 2.0);

    TVector3 InVec = track->post_step_point.momentum_direction,          
             direct1 = track->post_step_point.momentum_direction;

    if((sqrp0 > 0) && (sqrp1 > 0)) 
    {
        GenerateVectorInPhononScattering(sqrp0, sqrp1, hrq*hrq, InVec, direct1);
        track->post_step_point.momentum_direction = direct1;
    }
    else track->post_step_point.momentum_direction = RandomTVector3();

    dos = nullptr;
    cont_ait = nullptr;
}

void PhysicsList::CaptureParticleByCenter(Track* track, TString cname, TString ctype, Double_t cenergy)
{
    TrackPoint secondary_particle;
    secondary_particle.energy = cenergy;
    secondary_particle.point = track->post_step_point.point + track->post_step_point.momentum_direction * track->post_step_point.time_step * track->post_step_point.velocity;
    secondary_particle.time = track->post_step_point.time + track->post_step_point.time_step;
    Track* c_track = CreateCenter(secondary_particle, cname, ctype, track->post_step_point.volume, track->post_step_point.matprop);
    c_track->creator_process = track->particle_type + track->particle_name + " capture by " + ctype + cname;

    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(c_track);

    track->should_be_killed = true;
}

bool PhysicsList::GenerateVectorInInelasticScattering(Double_t sqrp0, Double_t sqrp1, Double_t sqrp2, TVector3 InVec, TVector3& direct1, TVector3& direct2)
{
    Double_t a = (sqrp0 + sqrp1 - sqrp2) / (2.0 * sqrt(sqrp0 * sqrp1) + 1E-80);
    Double_t b = (sqrp0 - sqrp1 + sqrp2) / (2.0 * sqrt(sqrp0 * sqrp2) + 1E-80);

    if((a > 1) || (a < -1) || (b > 1) || (b < -1)) return false;

    direct1 = TVector3(sqrt(1.0 - a*a), 0.0, a);
    direct2 = TVector3(-sqrt(1.0 - b*b), 0.0, b);

    Double_t fi = 2.0 * M_PI * gRandom->Rndm();

    TMatrixD mzc = TMatrixD(3,3);
    mzc[0] = {std::cos(fi), -std::sin(fi), 0.0};
    mzc[1] = {std::sin(fi),  std::cos(fi), 0.0};
    mzc[2] = {0.0         ,  0.0         , 1.0};

    Double_t sqz = sqrt(1.0 - InVec.Z() * InVec.Z());

    TMatrixD mrc = TMatrixD(3,3);
    mrc[0] = {InVec.X() * InVec.Z() / (sqz + 1E-8), -InVec.Y() / (sqz + 1E-8), InVec.X()};
    mrc[1] = {InVec.Y() * InVec.Z() / (sqz + 1E-8),  InVec.X() / (sqz + 1E-8), InVec.Y()};
    mrc[2] = {-sqz                                ,  0.0                     , InVec.Z()};

    direct1 = mrc * (mzc * direct1);
    direct2 = mrc * (mzc * direct2);

    direct1 = direct1 * pow(direct1.Mag(), -1.0);
    direct2 = direct2 * pow(direct2.Mag(), -1.0); 

    return true;
}

Double_t PhysicsList::GenerateInelasticScatteringEnergy(Double_t Eg, Double_t Ek, TString key)
{
    if(key == "e-")
    {
        Double_t Fmax = 1.0;
        if(Ek < 2.0*Eg) Fmax = 25.0 * sqrt(5.0) * pow(Eg, 2.0) * sqrt(Ek - Eg) / pow(Ek + 3.0 * Eg, 2.5);

        while(true)
        {
            Double_t x = Eg + gRandom->Rndm() * (Ek-Eg);
            Double_t F = 25.0 * sqrt(5.0) * pow(Eg, 2.0) * sqrt(x - Eg) / pow(x + 3.0 * Eg, 2.5) / Fmax;
            Double_t q = gRandom->Rndm();
    
            if(q <= F){return x;}
        }
    }
    else if(key == "h")
    {
        return Eg + gRandom->Rndm() * (Ek-Eg);
    }
    else return Eg;
}

void PhysicsList::InelasticScattering(Track* track, Int_t index)
{       
    Double_t Eg = track->post_step_point.matprop->GetBandGap();

    LinearInterpolation* dos = track->post_step_point.matprop->GetDensityOfStates();

    Container_imfp* cont_imfp = track->post_step_point.matprop->GetListOfInelasticMeanFreePath(track->particle_name)[index];
    
    Double_t Ek0 = track->post_step_point.kinetic_energy,
             E0 = track->post_step_point.energy;

    Track *e_track, *h_track;
    TVector3 direct1, direct2;
    
    Int_t while_iter = 1;
    while(true)
    {
        if(while_iter > 1000) {std::cout << "PhysicsList: Warning - To much iterations for an inelastic scattering event!\n"; return;}

        Double_t Threshold_energy = cont_imfp->imfp->GetXmin();
        if((Threshold_energy - Eg) < 0) std::cout << "PhysicsList: Warning - Wrong Eg or threshold! \n";

        Double_t hw = (cont_imfp->els == nullptr) ? GenerateInelasticScatteringEnergy(Eg, track->post_step_point.kinetic_energy - (Threshold_energy - Eg), track->particle_name) + (Threshold_energy - Eg) : 
                                                    RandomValue(cont_imfp->els, cont_imfp->els->GetXmin(), (track->post_step_point.kinetic_energy < cont_imfp->els->GetXmax()) ? track->post_step_point.kinetic_energy : cont_imfp->els->GetXmax());
        Double_t Eh = 0;

        if(Threshold_energy <= Eg)
        {
            Double_t Eh_min = ((Eg - hw) > dos->GetXmin()) ? (Eg - hw) : dos->GetXmin(); 
            Eh = RandomValue(dos, Eh_min, 0);
        }
        else {Eh = -(Threshold_energy - Eg);}

        Double_t Ek1 = track->post_step_point.kinetic_energy - hw,
                 Eke2 = hw - Eg + Eh;

        track->post_step_point.kinetic_energy = Ek1 + Eke2;
        UpdateMass(track);
        UpdateVelocity(track);

        Double_t sqrp0 = pow(CalculateMomentum(track->mass, track->post_step_point.kinetic_energy), 2.0);

        track->post_step_point.kinetic_energy = Ek1;
        track->post_step_point.energy = Eg + Ek1;
        UpdateMass(track);
        UpdateVelocity(track);

        Double_t sqrp1 = pow(CalculateMomentum(track->mass, track->post_step_point.kinetic_energy), 2.0);

        TrackPoint secondary_particle;
        secondary_particle.kinetic_energy = Eke2;
        secondary_particle.momentum_direction = RandomTVector3();
        secondary_particle.point = track->post_step_point.point + track->post_step_point.momentum_direction * track->post_step_point.time_step * track->post_step_point.velocity;
        secondary_particle.time = track->post_step_point.time + track->post_step_point.time_step;
        Track* loc_e_track = CreateElectron(secondary_particle, "", track->post_step_point.volume, track->post_step_point.matprop);
        loc_e_track->creator_process = "Inelastic scattering";

        Double_t sqrp2 = pow(CalculateMomentum(loc_e_track->mass, loc_e_track->post_step_point.kinetic_energy), 2.0);

        TrackPoint secondary_particle2;
        secondary_particle2.energy = Eh;
        secondary_particle2.momentum_direction = RandomTVector3();
        secondary_particle2.point = track->post_step_point.point + track->post_step_point.momentum_direction * track->post_step_point.time_step * track->post_step_point.velocity;
        secondary_particle2.time = track->post_step_point.time + track->post_step_point.time_step;
        Track* loc_h_track = CreateHole(secondary_particle2, "", track->post_step_point.volume, track->post_step_point.matprop);
        loc_h_track->creator_process = "Inelastic scattering";
        loc_h_track->track_id += 1;

        TVector3 InVec = track->post_step_point.momentum_direction;

        if(GenerateVectorInInelasticScattering(sqrp0, sqrp1, sqrp2, InVec, direct1, direct2)) 
        {
            e_track = loc_e_track;
            h_track = loc_h_track;

            loc_e_track = nullptr;
            loc_h_track = nullptr; 

            break;
        }

        track->post_step_point.kinetic_energy = Ek0;
        track->post_step_point.energy = E0;
        UpdateMass(track);
        UpdateVelocity(track);

        delete loc_e_track;
        delete loc_h_track;

        while_iter++;
    }

    track->post_step_point.momentum_direction = direct1;
    e_track->post_step_point.momentum_direction = direct2;
    e_track->pre_step_point.momentum_direction = direct2;

    TVector3 p1 = CalculateMomentum(e_track->mass,  e_track->post_step_point.kinetic_energy) * e_track->post_step_point.momentum_direction,
             p2 = CalculateMomentum(h_track->mass,  h_track->post_step_point.kinetic_energy) * h_track->post_step_point.momentum_direction;
    Double_t Ekc = const_c * sqrt((p1+p2).Mag2() + pow(e_track->mass + h_track->mass, 2.0) * const_c * const_c) - (e_track->mass + h_track->mass) * const_c * const_c;

    Double_t Ekeh = 0;
    if(e_track->post_step_point.kinetic_energy > 0) Ekeh += e_track->post_step_point.kinetic_energy;
    if(h_track->post_step_point.kinetic_energy > 0) Ekeh += h_track->post_step_point.kinetic_energy;

    Double_t disp_length = const_k / track->post_step_point.matprop->GetDielectricPermittivity() * const_e * const_e / (Ekeh - Ekc + unit_meV);
    if(disp_length < 2.0*unit_nm) disp_length = 2.0*unit_nm;
    TVector3 disp_vec = disp_length * RandomTVector3();

    e_track->pre_step_point.point += disp_vec;
    e_track->post_step_point.point += disp_vec;
    h_track->pre_step_point.point -= disp_vec;
    h_track->post_step_point.point -= disp_vec;

    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(e_track);
    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(h_track);

    dos = nullptr;
    cont_imfp = nullptr;
    e_track = nullptr;
    h_track = nullptr; 
}

void PhysicsList::InelasticHoleScattering(Track* track)
{    
    Double_t Eg = track->post_step_point.matprop->GetBandGap();
       
    LinearInterpolation* dos = track->post_step_point.matprop->GetDensityOfStates();

    Double_t Eh0 = track->post_step_point.energy, 
             Eh1 = 0, Eh2 = 0, Eke = 0, Ekeh = 0;

    if(track->post_step_point.energy >= dos->GetXmin())
    {
        Double_t hw = GenerateInelasticScatteringEnergy(Eg, -track->post_step_point.energy, track->particle_name);
        Double_t Eh_min = ((Eg - hw) > dos->GetXmin()) ? (Eg - hw) : dos->GetXmin();

        Eh1 = track->post_step_point.energy + hw,
        Eh2 = RandomValue(dos, Eh_min, 0),
        Eke = hw - Eg + Eh2;   
        
        if(Eh2 < 0) Ekeh -= Eh2;
        if(Eke > 0) Ekeh += Eke;
    }
    else if(track->post_step_point.energy < dos->GetXmin())
    {
        std::vector<Container_imfp*> inmfp = track->post_step_point.matprop->GetListOfInelasticMeanFreePath("e-");
        if(inmfp.size() == 0) return;

        std::vector<Double_t> core_energies;
        for (auto it = inmfp.begin() ; it != inmfp.end(); ++it)
        {
            if(-((*it)->imfp->GetXmin() - Eg) > track->post_step_point.energy)
            {
                core_energies.push_back(-((*it)->imfp->GetXmin() - Eg));
            }
        }
        if(core_energies.size() == 0) return;

        Int_t while_iter = 1;
        while (true)
        {
            if(while_iter > 1000) {std::cout << "PhysicsList: Warning - To much iterations for an inelastic hole scattering event!\n"; return;}
            
            Int_t j1 = gRandom->Integer(core_energies.size());
            Int_t j2 = gRandom->Integer(core_energies.size());

            if(core_energies[j1] == 0) {Eh1 = RandomValue(dos, dos->GetXmin(), 0);} else {Eh1 = core_energies[j1];}
            if(core_energies[j2] == 0) {Eh2 = RandomValue(dos, dos->GetXmin(), 0);} else {Eh2 = core_energies[j2];}

            Eke = -Eh0 + Eh1 + Eh2 - Eg;

            if(Eke > 0) break;

            while_iter++;
        }

        if(Eke > 0) Ekeh += Eke;
    }

    track->post_step_point.kinetic_energy = -Eh1;
    track->post_step_point.energy = Eh1;
    UpdateMass(track);
    UpdateVelocity(track);

    TrackPoint secondary_particle;
    secondary_particle.kinetic_energy = Eke;
    secondary_particle.momentum_direction = RandomTVector3();
    secondary_particle.point = track->post_step_point.point;
    secondary_particle.time = track->post_step_point.time + track->post_step_point.time_step;
    Track *e_track = CreateElectron(secondary_particle, "", track->post_step_point.volume, track->post_step_point.matprop);
    e_track->creator_process = "Inelastic scattering";

    TrackPoint secondary_particle2;
    secondary_particle2.energy = Eh2;  
    secondary_particle2.momentum_direction = RandomTVector3();
    secondary_particle2.point = track->post_step_point.point;
    secondary_particle2.time = track->post_step_point.time + track->post_step_point.time_step;
    Track *h_track = CreateHole(secondary_particle2, "", track->post_step_point.volume, track->post_step_point.matprop);
    h_track->creator_process = "Inelastic scattering";
    h_track->track_id += 1;

    TVector3 p1 = CalculateMomentum(e_track->mass,  e_track->post_step_point.kinetic_energy) * e_track->post_step_point.momentum_direction,
             p2 = CalculateMomentum(h_track->mass,  h_track->post_step_point.kinetic_energy) * h_track->post_step_point.momentum_direction;
    Double_t Ekc = const_c * sqrt((p1+p2).Mag2() + pow(e_track->mass + h_track->mass, 2.0) * const_c * const_c) - (e_track->mass + h_track->mass) * const_c * const_c;

    Double_t disp_length = const_k / track->post_step_point.matprop->GetDielectricPermittivity() * const_e * const_e / (Ekeh - Ekc + unit_meV);
    if(disp_length < 2.0*unit_nm) disp_length = 2.0*unit_nm;
    TVector3 disp_vec = disp_length * RandomTVector3();

    e_track->pre_step_point.point += disp_vec;
    e_track->post_step_point.point += disp_vec;
    h_track->pre_step_point.point -= disp_vec;
    h_track->post_step_point.point -= disp_vec;

    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(e_track);
    RunManager::getInstance()->GetTrackPropagation()->AddNewTrack(h_track);

    dos = nullptr;
}

void PhysicsList::DrawTrack(Track* track) 
{       
    if(track->visible)
    {
        Int_t num_track_point = track->geo_track->GetNpoints();
        Double_t xx, yy, zz, tt;
        num_track_point = track->geo_track->GetPoint(num_track_point-1, xx, yy, zz, tt);
        TVector3 last_point(xx, yy, zz);
        
        if((last_point - geo_size_factor * track->post_step_point.point).Mag() > geo_size_factor * 1.0 * unit_nm)
        {
            track->geo_track->AddPoint(geo_size_factor * track->post_step_point.point.X(), geo_size_factor * track->post_step_point.point.Y(), geo_size_factor * track->post_step_point.point.Z(), geo_time_factor * track->post_step_point.time);
        }
    }
}

void PhysicsList::ParticleStep(std::vector<Track*> track_list)
{    
    EH_GLOBAL_TIME += TIME_CUT_eh;
    
    for (auto &it : track_list)
    {
        if((it->killed == false) && !(((it->particle_name == "e-") || (it->particle_name == "h")) && (it->post_step_point.time >= EH_GLOBAL_TIME)))
        {
            it->pre_step_point = it->post_step_point;

            it->post_step_point.time += it->post_step_point.time_step;
            
            if(it->post_step_point.velocity > 0)
            {
                it->post_step_point.point = it->post_step_point.point + it->post_step_point.momentum_direction * it->post_step_point.time_step * it->post_step_point.velocity;

                bool cut_flag = (it->particle_name != "photon") ? ApplyPeriodicity(it->post_step_point.point) : false;
                if(cut_flag) CutTrackLine(it);

                TString volume_name = GetTGeoVolumeAtPoint(it->post_step_point.point);

                if(volume_name != "") 
                {   
                    it->post_step_point.volume = volume_name;

                    MaterialProperties *mp = RunManager::getInstance()->GetDetectorConstruction()->GetMaterialProperties(volume_name);
                    if(mp != nullptr)
                    {
                        it->post_step_point.matprop = mp;
                    }
                }
                else
                {
                    it->post_step_point.volume = "";
                    it->post_step_point.matprop = nullptr;
                    it->killed = true; 

                    continue;
                }

                if(!cut_flag) DrawTrack(it);
            }

           // if((it->post_step_point.time > 1.0 * unit_ns) && (it->particle_name != "photon")) {it->killed = true;}
           // if((it->post_step_point.time > 100.0 * unit_ns) && (it->particle_name == "photon")) {it->killed = true;}
        }

        if(it->should_be_killed == true) {it->killed = true; it->should_be_killed = false;}
    }
}

void PhysicsList::WhichProcesse(Track* track)
{        
    if(track->post_step_point.matprop == nullptr) return;

    std::vector<Double_t> prob_list;
    std::vector<std::pair<TString,Int_t>> index_list;

    LinearInterpolation* dos = track->post_step_point.matprop->GetDensityOfStates();
    if(dos != nullptr)
    {
        std::vector<Container_imfp*> imfp = track->post_step_point.matprop->GetListOfInelasticMeanFreePath(track->particle_name);
        for (auto it = imfp.begin(); it != imfp.end(); ++it)
        {
            if((track->post_step_point.kinetic_energy >= (*it)->imfp->GetXmin()) && (track->post_step_point.kinetic_energy <= (*it)->imfp->GetXmax()))
            {
                Double_t lambda = (*it)->imfp->GetValue(track->post_step_point.kinetic_energy);
                prob_list.push_back(track->post_step_point.velocity/lambda);
                index_list.push_back(std::make_pair(TString("imfp"), std::distance(imfp.begin(), it)));
            }
        }
        imfp.clear();   

        if((track->particle_name == "h") && (-track->post_step_point.energy > track->post_step_point.matprop->GetBandGap())) InelasticHoleScattering(track); ///////////////////////

        if(!((track->particle_name == "h") && (track->post_step_point.energy < dos->GetXmin())))
        {
            std::vector<Container_ait*> termal_ait = track->post_step_point.matprop->GetListOfThermalAverageInteractionTime(track->particle_name);
            for (auto it = termal_ait.begin(); it != termal_ait.end(); ++it)
            {
                if((track->post_step_point.kinetic_energy >= (*it)->ait->GetXmin()) && (track->post_step_point.kinetic_energy <= (*it)->ait->GetXmax()))
                {
                    Double_t tau = (*it)->ait->GetValue(track->post_step_point.kinetic_energy);
                    prob_list.push_back(1.0/tau);
                    index_list.push_back(std::make_pair(TString("termal_ait"), std::distance(termal_ait.begin(), it)));
                }
            }
            termal_ait.clear();
        }
    
        dos = nullptr;
    }

    std::vector<Container_emfp*> el_ait = track->post_step_point.matprop->GetListOfElasticAverageInteractionTime(track->particle_name);
    for (auto it = el_ait.begin(); it != el_ait.end(); ++it)
    {
        if((track->post_step_point.kinetic_energy >= (*it)->ait->GetXmin()) && (track->post_step_point.kinetic_energy <= (*it)->ait->GetXmax()))
        {
            Double_t tau = (*it)->ait->GetValue(track->post_step_point.kinetic_energy);
            prob_list.push_back(1.0/tau);
            index_list.push_back(std::make_pair(TString("el_ait"), std::distance(el_ait.begin(), it)));
        }
    }
    el_ait.clear();

    std::vector<CenterProperties*> clist = track->post_step_point.matprop->GetCenterList();
    std::vector<std::tuple<TString,TString,Double_t>> center_prop_1;
    for (auto &it1 : clist)
    {
        std::vector<Container_ccr*> ccrates = it1->GetCaptureRates(track->particle_type + track->particle_name);
        for (auto &it2 : ccrates)
        {
            if((track->post_step_point.kinetic_energy >= it2->rate->GetXmin()) && (track->post_step_point.kinetic_energy <= it2->rate->GetXmax()))
            {
                Double_t crate = it2->rate->GetValue(track->post_step_point.kinetic_energy);
                prob_list.push_back(crate);
                index_list.push_back(std::make_pair(TString("carrier_capture"), (int)center_prop_1.size()));
                center_prop_1.push_back(std::make_tuple(it2->center_name, it2->center_type, it1->GetLevel(it2->center_type)->shell_energy));
            }
        }
        ccrates.clear();   
    }
    clist.clear();   

    CenterProperties* center = track->post_step_point.matprop->GetCenter(track->particle_name);
    std::vector<std::tuple<TString,Double_t,bool>> center_prop_2;
    if(center != nullptr)
    {
        shell* sh = center->GetLevel(track->particle_type);
        if(sh != nullptr)
        {
            for (auto &it : sh->radiation_times)
            {
                shell* loc_sh = center->GetLevel(it.first);
                if(loc_sh != nullptr)
                {
                    TString ctype = (loc_sh->shell_name == center->GetGroundState()) ? "ground" : loc_sh->shell_name;

                    prob_list.push_back(1.0/it.second);
                    index_list.push_back(std::make_pair(TString("center_decay"), (int)center_prop_2.size()));
                    center_prop_2.push_back(std::make_tuple(ctype, loc_sh->shell_energy, true));

                    loc_sh = nullptr;
                }
            }
            
            for (auto &it : sh->nonradiation_times)
            {
                shell* loc_sh = center->GetLevel(it.first);
                if(loc_sh != nullptr)
                {
                    TString ctype = (loc_sh->shell_name == center->GetGroundState()) ? "ground" : loc_sh->shell_name;

                    prob_list.push_back(1.0/it.second);
                    index_list.push_back(std::make_pair(TString("center_decay"), (int)center_prop_2.size()));
                    center_prop_2.push_back(std::make_tuple(ctype, loc_sh->shell_energy, false));

                    loc_sh = nullptr;
                }
            }
            sh = nullptr;
        }
        center = nullptr;
    }

    if(track->particle_name == "ex")
    {
        LinearInterpolation* decaytau = track->post_step_point.matprop->GetRadiationDecayTimeOfExciton();
        if(decaytau != nullptr)
        {
            if((track->post_step_point.kinetic_energy >= decaytau->GetXmin()) && (track->post_step_point.kinetic_energy <= decaytau->GetXmax()))
            {
                Double_t tau = decaytau->GetValue(track->post_step_point.kinetic_energy);
                prob_list.push_back(1.0/tau);
                index_list.push_back(std::make_pair(TString("ex_decay"), 0));
            }
            decaytau = nullptr;
        }

        Container_exadt* adt = track->post_step_point.matprop->GetExcitonDissociationTime();
        if(adt != nullptr)
        {
            if((track->post_step_point.kinetic_energy >= adt->ait->GetXmin()) && (track->post_step_point.kinetic_energy <= adt->ait->GetXmax()))
            {
                Double_t tau = adt->ait->GetValue(track->post_step_point.kinetic_energy);
                prob_list.push_back(1.0/tau);
                index_list.push_back(std::make_pair(TString("ex_dissociation"), 0));
            }
            adt = nullptr;
        }
    }

    if(prob_list.size() == 0) return;

    Double_t inv_total_tau = 0;
    for (auto it = prob_list.begin(); it != prob_list.end(); ++it)
    {
        inv_total_tau += (*it);
    }

    Double_t p = gRandom->Rndm();
    Double_t prob = 0;
    
    if((track->particle_name == "e-") || (track->particle_name == "h"))
    {
        prob = 1.0 - TMath::Exp(-track->post_step_point.time_step * inv_total_tau);
    }
    else
    {
        prob = 1.0;
        track->post_step_point.time_step = log(1.0 / (gRandom->Rndm() + 1E-10)) / (inv_total_tau + 1E-10); 
    }

    if(p < prob)
    {
        Int_t index = DiscreteRandomValue(prob_list);

        if(index_list[index].first == "imfp")
        {
            InelasticScattering(track, index_list[index].second);
        }
        else if(index_list[index].first == "termal_ait")
        {
            ThermalScattering(track, index_list[index].second);
        }
        else if(index_list[index].first == "el_ait")
        {
            ElasticScattering(track, index_list[index].second);
        }
        else if(index_list[index].first == "ex_decay")
        {
            RadiationDecayOfExciton(track);
        }
        else if(index_list[index].first == "ex_dissociation")
        {
            ExcitonDissociation(track);
        }
        else if(index_list[index].first == "carrier_capture")
        {
            CaptureParticleByCenter(track, std::get<0>(center_prop_1[index_list[index].second]), 
                                           std::get<1>(center_prop_1[index_list[index].second]), 
                                           std::get<2>(center_prop_1[index_list[index].second]));
        }
        else if(index_list[index].first == "center_decay")
        {           
            DecayOfCenter(track, track->particle_name, std::get<0>(center_prop_2[index_list[index].second]), 
                                                       std::get<1>(center_prop_2[index_list[index].second]), 
                                                       std::get<2>(center_prop_2[index_list[index].second]));
        }
        else return;
    }

    prob_list.clear();
    index_list.clear();
    center_prop_1.clear();
    center_prop_2.clear();
}

void PhysicsList::ApplySingleTrackProcesses(std::vector<Track*> track_list)
{
    std::vector<Track*> task_list;

    task_list.reserve(track_list.size());

    for (auto &it : track_list)
    {
        if((it->killed == false) && (it->should_be_killed == false) && !(((it->particle_name == "e-") || (it->particle_name == "h")) && (it->post_step_point.time >= EH_GLOBAL_TIME + TIME_CUT_eh))) task_list.push_back(it);
    }

    task_list.shrink_to_fit();

    if(task_list.size() > 0)
    {
        #pragma omp parallel for
        for (auto &it : task_list)
        {        
            WhichProcesse(it);
        }

        task_list.clear();
    }
}

void PhysicsList::ApplyMultiTrackProcesses(std::vector<Track*> track_list)
{
    std::vector<Track*> task_list;

    task_list.reserve(track_list.size());

    for (auto &it : track_list)
    {
        if((it->killed == false) && (it->should_be_killed == false))
        {
            if(it->post_step_point.matprop != nullptr)
            {        
                if((it->particle_name == "e-") && (it->post_step_point.time < EH_GLOBAL_TIME + TIME_CUT_eh))
                {
                    LinearInterpolation* dos = it->post_step_point.matprop->GetDensityOfStates();
                    if(dos != nullptr)
                    {
                        task_list.push_back(it);
                    }
                }
                else if((it->particle_name == "h") && (it->post_step_point.time < EH_GLOBAL_TIME + TIME_CUT_eh))
                {
                    LinearInterpolation* dos = it->post_step_point.matprop->GetDensityOfStates();
                    if(dos != nullptr)
                    {
                        if(it->post_step_point.energy > dos->GetXmin()) task_list.push_back(it);
                    }
                }
            }
        }
    }

    task_list.shrink_to_fit();
    
    if(task_list.size() > 0) ElectronHoleInteraction(task_list);
}