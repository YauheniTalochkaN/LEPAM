#include "MaterialProperties.hh"
#include "RunManager.hh"

MaterialProperties::MaterialProperties()
{
    
}

MaterialProperties::~MaterialProperties()
{
    if(Density_Of_States != nullptr) delete Density_Of_States;
    if(Radiation_Decay_Time_Of_Exciton != nullptr) delete Radiation_Decay_Time_Of_Exciton;
    if(Electron_Hole_Radiation_Recombination_Time != nullptr) delete Electron_Hole_Radiation_Recombination_Time;
    if(Exciton_Dissociation_Time != nullptr) delete Exciton_Dissociation_Time;

    if(Effective_Mass.size() > 0)
    {
        for (auto &it : Effective_Mass)
        {
            if(it != nullptr) delete it;
        }
        Effective_Mass.clear();
    }

    if(Group_Velocity.size() > 0)
    {
        for (auto &it : Group_Velocity)
        {
            if(it != nullptr) delete it;
        }
        Group_Velocity.clear();
    }

    if(Inelastic_Mean_Free_Paths.size() > 0)
    {
        for (auto &it : Inelastic_Mean_Free_Paths)
        {
            if(it != nullptr) delete it;
        }
        Inelastic_Mean_Free_Paths.clear();
    }

    if(Elastic_Mean_Free_Path.size() > 0)
    {
        for (auto &it : Elastic_Mean_Free_Path)
        {
            if(it != nullptr) delete it;
        }
        Elastic_Mean_Free_Path.clear();
    }

    if(Thermal_Average_Interaction_Time.size() > 0)
    {
        for (auto &it : Thermal_Average_Interaction_Time)
        {
            if(it != nullptr) delete it;
        }
        Thermal_Average_Interaction_Time.clear();
    }

    if(List_Of_Centers.size() > 0)
    {
        for (auto &it : List_Of_Centers)
        {
            if(it != nullptr) delete it;
        }
        List_Of_Centers.clear();
    }
}

void MaterialProperties::SetBandGap(Double_t Eg)
{
    if(Eg >= 0) Band_gap = Eg;
}

void MaterialProperties::SetDielectricPermittivity(Double_t eps)
{
    if(eps > 0) Dielectric_Permittivity = eps;
}

void MaterialProperties::SetEffectiveMass(Container* m)
{
    if((m->func != nullptr) && (m->full_particle_name != ""))
    Effective_Mass.push_back(m);
}

void MaterialProperties::SetUnitCellVolume(Double_t V)
{
    if(V > 0) Unit_Cell_Volume = V;
}

void MaterialProperties::SetDensityOfStates(LinearInterpolation* dos)
{
    if(dos != nullptr) Density_Of_States = dos;
}

void MaterialProperties::SetRadiationDecayTimeOfExciton(LinearInterpolation* tau)
{
    if(tau != nullptr) Radiation_Decay_Time_Of_Exciton = tau;
}

void MaterialProperties::SetElectronHoleRadiationRecombinationTime(LinearInterpolation2D* tau)
{
    if(tau != nullptr) Electron_Hole_Radiation_Recombination_Time = tau;
}

void MaterialProperties::SetGroupVelocity(Container* v)
{
    if((v->func != nullptr) && (v->full_particle_name != ""))
    Group_Velocity.push_back(v);
}

void MaterialProperties::SetInelasticMeanFreePath(Container_imfp* lambda)
{
    if((lambda->imfp != nullptr) /*&& (lambda->els != nullptr)*/ && (lambda->full_particle_name != ""))
    Inelastic_Mean_Free_Paths.push_back(lambda);
}

void MaterialProperties::SetElasticMeanFreePath(Container_emfp* lambda)
{
    if((lambda->ait != nullptr) /*&& (lambda->mls != nullptr)*/ && (lambda->full_particle_name != ""))
    Elastic_Mean_Free_Path.push_back(lambda);
}

void MaterialProperties::SetThermalAverageInteractionTime(Container_ait* tau)
{
    if((tau->ait != nullptr) && (tau->displow != nullptr) && (tau->mls != nullptr) && (tau->full_particle_name != ""))
    Thermal_Average_Interaction_Time.push_back(tau);
}

void MaterialProperties::SetExcitonDissociationTime(Container_exadt* tau)
{
    if((tau->ait != nullptr) /*&& (tau->des != nullptr)*/) Exciton_Dissociation_Time = tau;
}

void MaterialProperties::SetCenterProperties(CenterProperties* obj)
{
    if(obj != nullptr)
    List_Of_Centers.push_back(obj);
}

Double_t MaterialProperties::GetBandGap()
{
    return Band_gap;
}

Double_t MaterialProperties::GetDielectricPermittivity()
{
    return Dielectric_Permittivity;
}

LinearInterpolation* MaterialProperties::GetEffectiveMass(TString particle_name)
{
    for (auto &it : Effective_Mass)
    {
        if(it->full_particle_name == particle_name) return it->func;
    }

    return nullptr;
}

Double_t MaterialProperties::GetUnitCellVolume()
{
    return Unit_Cell_Volume;
}

LinearInterpolation* MaterialProperties::GetDensityOfStates()
{
    return Density_Of_States;
}

LinearInterpolation* MaterialProperties::GetRadiationDecayTimeOfExciton()
{
    return Radiation_Decay_Time_Of_Exciton;
}

LinearInterpolation2D* MaterialProperties::GetElectronHoleRadiationRecombinationTime()
{
    return Electron_Hole_Radiation_Recombination_Time;
}

LinearInterpolation* MaterialProperties::GetGroupVelocity(TString particle_name)
{
    for (auto &it : Group_Velocity)
    {
        if(it->full_particle_name == particle_name) return it->func;
    }

    return nullptr;
}

std::vector<Container_imfp*> MaterialProperties::GetListOfInelasticMeanFreePath(TString particle_name)
{
    std::vector<Container_imfp*> list_IMFP;

    for (auto &it : Inelastic_Mean_Free_Paths)
    {
        if(it->full_particle_name == particle_name) list_IMFP.push_back(it);
    }

    return list_IMFP;
}

std::vector<Container_emfp*> MaterialProperties::GetListOfElasticAverageInteractionTime(TString particle_name)
{
    std::vector<Container_emfp*> list_EMFP;

    for (auto &it : Elastic_Mean_Free_Path)
    {
        if(it->full_particle_name == particle_name) list_EMFP.push_back(it);
    }

    return list_EMFP;
}

std::vector<Container_ait*> MaterialProperties::GetListOfThermalAverageInteractionTime(TString particle_name)
{
    std::vector<Container_ait*> list_AIT;
    
    for (auto &it : Thermal_Average_Interaction_Time)
    {
        if(it->full_particle_name == particle_name) list_AIT.push_back(it);
    }

    return list_AIT;
}

Container_exadt* MaterialProperties::GetExcitonDissociationTime()
{
    return Exciton_Dissociation_Time;
}

std::vector<CenterProperties*> MaterialProperties::GetCenterList()
{    
    return List_Of_Centers;
}

CenterProperties* MaterialProperties::GetCenter(TString cname)
{
    for (auto &it : List_Of_Centers)
    {
        if(it->GetName() == cname) return it;
    }

    return nullptr;
}
