#include "LinearInterpolation2D.hh"

LinearInterpolation2D* LinearInterpolation2D::Initialization(TString innam, Double_t unit1, Double_t unit2, Double_t unit3)
{
    FILE *file = fopen(innam, "r");
    if (!file)
    {
        std::cout << "\n" << "LinearInterpolation2D: Wrong name of the input file of data for interpolation! \n" << "File name: " << innam << "\n\n";
        return nullptr;
    }
    else
    {
        fseek(file, 0, SEEK_END);
        if(ftell(file) == 0) 
        {
            std::cout << "LinearInterpolation2D: Empty file - " << innam << "\n"; 
            fclose(file); 
            return nullptr;
        }
        {
            fseek(file, 0, SEEK_SET);
            return new LinearInterpolation2D(file, unit1, unit2, unit3);
        }
    }
}

LinearInterpolation2D::LinearInterpolation2D(FILE* file, Double_t unit1, Double_t unit2, Double_t unit3)
{     
    Double_t xx = 0;
    std::vector<Double_t> list_yy, list_zz;
    
    char STR[256];
    while (STR == fgets(STR, 255, file))
    {
        double v1, v2;
        int nit = sscanf(STR, "%lf %lf", &v1, &v2);
        if(nit > 1)
        {
            list_yy.push_back(v1 * unit2);
            list_zz.push_back(v2 * unit3);
        }
        else if(nit > 0)
        {
            if(list_yy.size() > 0)
            {
                LinInter1DList.push_back(new projectionX(xx, new LinearInterpolation(list_yy, list_zz)));
                list_yy.clear();
                list_zz.clear();
            }
            
            xx = v1 * unit1;
        }
    }

    if(list_yy.size() > 0)
    {
        LinInter1DList.push_back(new projectionX(xx, new LinearInterpolation(list_yy, list_zz)));
        list_yy.clear();
        list_zz.clear();
    }

    fclose(file); 
}

LinearInterpolation2D::LinearInterpolation2D(std::vector<projectionX*> data)
{
    for (auto &it : data)
    {
        LinInter1DList.push_back(it);
    }
}

LinearInterpolation2D::~LinearInterpolation2D()
{
    if(LinInter1DList.size() == 0) return;
    
    for (auto &it : LinInter1DList)
    {
        delete it;
    }

    LinInter1DList.clear();
}

Double_t LinearInterpolation2D::GetValue(Double_t xx, Double_t yy)
{    
    if(LinInter1DList.size() == 0) {std::cout << "LinearInterpolation2D: Empty  arrays!\n"; return -1;}
    if(LinInter1DList.size() == 1) return LinInter1DList[0]->slice->GetValue(yy);

    if(xx >= GetXmax()) return LinInter1DList[LinInter1DList.size()-1]->slice->GetValue(yy);
    if(xx <= GetXmin()) return LinInter1DList[0]->slice->GetValue(yy);
    
    int k = 0;
    for (int i = 0; i < (int)(LinInter1DList.size()-1); i++)
    {
        if(xx < LinInter1DList[i+1]->X) {k = i; break;}
    }
    
    return LinInter1DList[k]->slice->GetValue(yy) + (xx - LinInter1DList[k]->X) * ((LinInter1DList[k+1]->slice->GetValue(yy) - LinInter1DList[k]->slice->GetValue(yy))/(LinInter1DList[k+1]->X - LinInter1DList[k]->X));
}

LinearInterpolation* LinearInterpolation2D::GetXProjection(Double_t xx)
{
    if(LinInter1DList.size() == 0) {std::cout << "LinearInterpolation2D: Empty  arrays!\n"; return nullptr;}
    if(LinInter1DList.size() == 1) return LinInter1DList[0]->slice;

    if(xx >= GetXmax()) 
    {
        Int_t nbins = LinInter1DList[LinInter1DList.size()-1]->slice->GetNBins();
        Double_t *array_yy = LinInter1DList[LinInter1DList.size()-1]->slice->GetXMass(), *array_zz = LinInter1DList[LinInter1DList.size()-1]->slice->GetYMass();
        return new LinearInterpolation(array_yy, array_zz, nbins);
    }
    
    if(xx <= GetXmin()) 
    {   
        Int_t nbins = LinInter1DList[0]->slice->GetNBins();
        Double_t *array_yy = LinInter1DList[0]->slice->GetXMass(), *array_zz = LinInter1DList[0]->slice->GetYMass();
        return new LinearInterpolation(array_yy, array_zz, nbins);
    }
    
    int k = 0;
    for (int i = 0; i < (int)(LinInter1DList.size()-1); i++)
    {
        if(xx < LinInter1DList[i+1]->X) {k = i; break;}
    }
    
    Double_t y_min = 0.5 * (LinInter1DList[k]->slice->GetXmin() + LinInter1DList[k+1]->slice->GetXmin());
    Double_t y_max = 0.5 * (LinInter1DList[k]->slice->GetXmax() + LinInter1DList[k+1]->slice->GetXmax());
    Int_t nbins = (LinInter1DList[k]->slice->GetNBins() > LinInter1DList[k+1]->slice->GetNBins()) ? LinInter1DList[k]->slice->GetNBins() : LinInter1DList[k+1]->slice->GetNBins();

    Double_t array_yy[nbins], array_zz[nbins];
    for (int i = 0; i < nbins; i++)
    {
        array_yy[i] = y_min + (y_max-y_min) * i / (double)(nbins-1);
        array_zz[i] = LinInter1DList[k]->slice->GetValue(array_yy[i]) + (xx - LinInter1DList[k]->X) * ((LinInter1DList[k+1]->slice->GetValue(array_yy[i]) - LinInter1DList[k]->slice->GetValue(array_yy[i]))/(LinInter1DList[k+1]->X - LinInter1DList[k]->X));
    }   

    return new LinearInterpolation(array_yy, array_zz, nbins);
}

LinearInterpolation* LinearInterpolation2D::GetYProjection(Double_t yy)
{
    Int_t nbins = LinInter1DList.size();
    Double_t array_xx[nbins], array_zz[nbins]; 

    for (int i = 0; i < nbins; i++)
    {
        array_xx[i] = LinInter1DList[i]->X;
        array_zz[i] = LinInter1DList[i]->slice->GetValue(yy);
    } 

    return new LinearInterpolation(array_xx, array_zz, nbins); 
}

Double_t LinearInterpolation2D::GetXmin()
{    
    return LinInter1DList[0]->X;
}

Double_t LinearInterpolation2D::GetXmax()
{
    return LinInter1DList[LinInter1DList.size()-1]->X;
}

Double_t LinearInterpolation2D::GetYmin()
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

Double_t LinearInterpolation2D::GetYmax()
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

Double_t LinearInterpolation2D::GetZmin()
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

Double_t LinearInterpolation2D::GetZmax()
{
    Double_t zz_max[LinInter1DList.size()];
    int k=0;
    for (auto it = LinInter1DList.begin() ; it != LinInter1DList.end(); ++it)
    {
        zz_max[k] = (*it)->slice->GetYmax();
        k++;
    }
    
    return TMath::MaxElement(LinInter1DList.size(), zz_max);
}
