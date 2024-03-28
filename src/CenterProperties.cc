#include "CenterProperties.hh"

CenterProperties::CenterProperties(TString nam) : Name(nam)
{
    
}

CenterProperties::~CenterProperties()
{
    if(Level_List.size() > 0)
    {
       for (auto &it : Level_List)
       {
          delete it;
       }
       Level_List.clear();
    }
    
    if(Capture_Rates.size() > 0)
    {
       for (auto &it : Capture_Rates)
       {
          delete it;
       }
      Capture_Rates.clear();
    }
}

void CenterProperties::SetName(TString nam)
{
   Name = nam;
}

void CenterProperties::SetGroundState(TString gr)
{
   ground = gr;
}

void CenterProperties::SetLevel(shell* lv)
{
   if(lv != nullptr) Level_List.push_back(lv);
}

void CenterProperties::SetCaptureRate(Container_ccr* crate)
{
   if(crate != nullptr) Capture_Rates.push_back(crate);
}

TString CenterProperties::GetName()
{
   return Name;
}

TString CenterProperties::GetGroundState()
{
   return ground;
}

shell* CenterProperties::GetLevel(TString nam)
{
   for (auto &it : Level_List)
   {
       if(it->shell_name == nam) return it;
   }

   return nullptr;
}

std::vector<Container_ccr*> CenterProperties::GetCaptureRates(TString nam)
{
   std::vector<Container_ccr*> list_ccr;
   
   for (auto &it : Capture_Rates)
   {
       if(it->full_particle_name == nam) list_ccr.push_back(it);
   }

   return list_ccr;
}