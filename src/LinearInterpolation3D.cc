#include "LinearInterpolation3D.hh"

LinearInterpolation3D* LinearInterpolation3D::Initialization(TString innam, Double_t unit1, Double_t unit2, Double_t unit3, Double_t unit4)
{
    FILE *file = fopen(innam, "r");
    if (!file)
    {
        std::cout << "\n" << "LinearInterpolation3D: Wrong name of the input file of data for interpolation! \n" << "File name: " << innam << "\n\n";
        return nullptr;
    }
    else
    {
        fseek(file, 0, SEEK_END);
        if(ftell(file) == 0) 
        {
            std::cout << "LinearInterpolation3D: Empty file - " << innam << "\n"; 
            fclose(file); 
            return nullptr;
        }
        {
            fseek(file, 0, SEEK_SET);
            return new LinearInterpolation3D(file, unit1, unit2, unit3, unit4);
        }
    }
}

LinearInterpolation3D::LinearInterpolation3D(FILE* file, Double_t unit1, Double_t unit2, Double_t unit3, Double_t unit4)
{     
    Double_t xx = 0, ksi = 0, yy = 0;
    int prenit = 0;
    std::vector<Double_t> list_zz, list_pp;
    std::vector<projectionX*> pjZP;
    
    char STR[256];
    while (STR == fgets(STR, 255, file))
    {
        double v1, v2;
        int nit = sscanf(STR, "%lf %lf", &v1, &v2);
        if(nit > 1)
        {
            if(prenit == 1) yy = ksi * unit2;

            list_zz.push_back(v1 * unit3);
            list_pp.push_back(v2 * unit4);
        }
        else if(nit > 0)
        {
            if(list_zz.size() > 0)
            {
                pjZP.push_back(new projectionX(yy, new LinearInterpolation(list_zz, list_pp)));
                list_zz.clear();
                list_pp.clear();
            }
            
            if(prenit == 1) 
            {
                xx = ksi * unit1;

                if(pjZP.size() > 0) 
                {
                    LinInter2DList.push_back(new sliceYZ(xx, new LinearInterpolation2D(pjZP)));
                    pjZP.clear();
                }
            }

            ksi = v1;
        }
        prenit = nit;
    }

    if(list_zz.size() > 0)
    {
        pjZP.push_back(new projectionX(yy, new LinearInterpolation(list_zz, list_pp)));
        list_zz.clear();
        list_pp.clear();
    }
    
    if(pjZP.size() > 0) 
    {
        LinInter2DList.push_back(new sliceYZ(xx, new LinearInterpolation2D(pjZP)));
        pjZP.clear();
    }

    fclose(file); 
}

LinearInterpolation3D::~LinearInterpolation3D()
{
    if(LinInter2DList.size() == 0) return;
    
    for (auto &it : LinInter2DList)
    {
        delete it;
    }

    LinInter2DList.clear();
}

Double_t LinearInterpolation3D::GetValue(Double_t xx, Double_t yy, Double_t zz)
{    
    if(LinInter2DList.size() == 0) {std::cout << "LinearInterpolation3D: Empty  arrays!\n"; return -1;}
    if(LinInter2DList.size() == 1) return LinInter2DList[0]->surface->GetValue(yy, zz);

    if(xx >= GetXmax()) return LinInter2DList[LinInter2DList.size()-1]->surface->GetValue(yy, zz);
    if(xx <= GetXmin()) return LinInter2DList[0]->surface->GetValue(yy, zz);
    
    int k = 0;
    for (int i = 0; i < (int)(LinInter2DList.size()-1); i++)
    {
        if(xx < LinInter2DList[i+1]->X) {k = i; break;}
    }
    
    return LinInter2DList[k]->surface->GetValue(yy, zz) + (xx - LinInter2DList[k]->X) * ((LinInter2DList[k+1]->surface->GetValue(yy, zz) - LinInter2DList[k]->surface->GetValue(yy, zz))/(LinInter2DList[k+1]->X - LinInter2DList[k]->X));
}

Double_t LinearInterpolation3D::GetXmin()
{    
    return LinInter2DList[0]->X;
}

Double_t LinearInterpolation3D::GetXmax()
{
    return LinInter2DList[LinInter2DList.size()-1]->X;
}
/*
Double_t LinearInterpolation3D::GetYmin()
{
    Double_t yy_min[LinInter1DList.size()];
    int k=0;
    for (auto it = LinInter1DList.begin() ; it != LinInter1DList.end(); ++it)
    {
        yy_min[k] = (*it)->slice->GetXmin();
        k++;
    }
    
    return TMath::MinElement(LinInter1DList.size(), yy_min);
}

Double_t LinearInterpolation3D::GetYmax()
{
    Double_t yy_max[LinInter1DList.size()];
    int k=0;
    for (auto it = LinInter1DList.begin() ; it != LinInter1DList.end(); ++it)
    {
        yy_max[k] = (*it)->slice->GetXmax();
        k++;
    }
    
    return TMath::MaxElement(LinInter1DList.size(), yy_max);
}

Double_t LinearInterpolation3D::GetZmin()
{
    Double_t zz_min[LinInter1DList.size()];
    int k=0;
    for (auto it = LinInter1DList.begin() ; it != LinInter1DList.end(); ++it)
    {
        zz_min[k] = (*it)->slice->GetYmin();
        k++;
    }
    
    return TMath::MinElement(LinInter1DList.size(), zz_min);
}

Double_t LinearInterpolation3D::GetZmax()
{
    Double_t zz_max[LinInter1DList.size()];
    int k=0;
    for (auto it = LinInter1DList.begin() ; it != LinInter1DList.end(); ++it)
    {
        zz_max[k] = (*it)->slice->GetYmax();
        k++;
    }
    
    return TMath::MaxElement(LinInter1DList.size(), zz_max);
}*/
